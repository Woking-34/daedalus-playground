// ---------------------------------------------------------------------------------//
// Copyright (c) 2013, Regents of the University of Pennsylvania                    //
// All rights reserved.                                                             //
//                                                                                  //
// Redistribution and use in source and binary forms, with or without               //
// modification, are permitted provided that the following conditions are met:      //
//     * Redistributions of source code must retain the above copyright             //
//       notice, this list of conditions and the following disclaimer.              //
//     * Redistributions in binary form must reproduce the above copyright          //
//       notice, this list of conditions and the following disclaimer in the        //
//       documentation and/or other materials provided with the distribution.       //
//     * Neither the name of the <organization> nor the                             //
//       names of its contributors may be used to endorse or promote products       //
//       derived from this software without specific prior written permission.      //
//                                                                                  //
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  //
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    //
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           //
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY               //
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       //
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     //
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      //
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       //
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    //
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     //
//                                                                                  //
// Contact Tiantian Liu (ltt1598@gmail.com) if you have any questions.              //
//----------------------------------------------------------------------------------//

#pragma warning( disable : 4996)

#include "simulation.h"
#include "timer_wrapper.h"

Simulation::Simulation()

{
}

Simulation::~Simulation()
{
    clearConstraints();
}

void Simulation::Reset()
{    
    m_inertia_y.resize(m_mesh->m_system_dimension);
    m_external_force.resize(m_mesh->m_system_dimension);

    setupConstraints();
    SetReprefactorFlag();

    m_selected_attachment_constraint = NULL;
}

void Simulation::Update()
{
    // update inertia term
    calculateInertiaY();

    // update external force
    calculateExternalForce();

    // update cloth
    switch (m_integration_method)
    {
    case INTEGRATION_EXPLICIT_EULER:
        integrateExplicitEuler();
        break;
    case INTEGRATION_EXPLICIT_SYMPLECTIC:
        integrateExplicitSymplectic();
        break;
    case INTEGRATION_IMPLICIT_EULER_BARAFF_WITKIN:
        integrateImplicitBW();
        break;
    case INTEGRATION_GRADIENT_DESCENT:
    case INTEGRATION_NEWTON_DESCENT:
    case INTEGRATION_NWETON_DESCENT_PCG:
    case INTEGRATION_LOCAL_GLOBAL:
        integrateOptimizationMethod();
        break;
    }

    // Add collision detection here using pos_next;
    VectorX penetration = collisionDetection(m_mesh->m_current_positions);
    m_mesh->m_current_positions -= penetration;

    // update velocity and damp
    dampVelocity();
}

void Simulation::DrawConstraints(const VBO& vbos)
{
    for (std::vector<Constraint*>::iterator it = m_constraints.begin(); it != m_constraints.end(); ++it)
    {
        (*it)->Draw(vbos);
    }
}

ScalarType Simulation::TryToSelectAttachmentConstraint(const EigenVector3& p0, const EigenVector3& dir)
{
    ScalarType ray_point_dist;
    ScalarType min_dist = 100.0;
    AttachmentConstraint* best_candidate = NULL;

    bool current_state_on = false;
    for (std::vector<Constraint*>::iterator c = m_constraints.begin(); c != m_constraints.end(); ++c)
    {
        AttachmentConstraint* ac;
        if (ac = dynamic_cast<AttachmentConstraint*>(*c)) // is attachment constraint
        {
            ray_point_dist = ((ac->GetFixedPoint()-p0).cross(dir)).norm();
            if (ray_point_dist < min_dist)
            {
                min_dist = ray_point_dist;
                best_candidate = ac;
            }
        }
    }
    // exit if no one fits
    if (min_dist > DEFAULT_SELECTION_RADIUS)
    {
        UnselectAttachmentConstraint();

        return -1;
    }
    else
    {
        SelectAtttachmentConstraint(best_candidate);
        EigenVector3 fixed_point_temp = m_mesh->m_current_positions.block_vector(m_selected_attachment_constraint->GetConstrainedVertexIndex());

        return (fixed_point_temp-p0).dot(dir); // this is m_cached_projection_plane_distance
    }
}

bool Simulation::TryToToggleAttachmentConstraint(const EigenVector3& p0, const EigenVector3& dir)
{
    EigenVector3 p1;

    ScalarType ray_point_dist;
    ScalarType min_dist = 100.0;
    unsigned int best_candidate = 0;
    // first pass: choose nearest point
    for (unsigned int i = 0; i != m_mesh->m_vertices_number; i++)
    {
        p1 = m_mesh->m_current_positions.block_vector(i);
        
        ray_point_dist = ((p1-p0).cross(dir)).norm();
        if (ray_point_dist < min_dist)
        {
            min_dist = ray_point_dist;
            best_candidate = i;
        }
    }
    for (std::vector<Constraint*>::iterator c = m_constraints.begin(); c != m_constraints.end(); ++c)
    {
        AttachmentConstraint* ac;
        if (ac = dynamic_cast<AttachmentConstraint*>(*c)) // is attachment constraint
        {
            ray_point_dist = ((ac->GetFixedPoint()-p0).cross(dir)).norm();
            if (ray_point_dist < min_dist)
            {
                min_dist = ray_point_dist;
                best_candidate = ac->GetConstrainedVertexIndex();
            }
        }
    }
    // exit if no one fits
    if (min_dist > DEFAULT_SELECTION_RADIUS)
    {
        return false;
    }
    // second pass: toggle that point's fixed position constraint
    bool current_state_on = false;
    for (std::vector<Constraint*>::iterator c = m_constraints.begin(); c != m_constraints.end(); ++c)
    {
        AttachmentConstraint* ac;
        if (ac = dynamic_cast<AttachmentConstraint*>(*c)) // is attachment constraint
        {
            if (ac->GetConstrainedVertexIndex() == best_candidate)
            {
                current_state_on = true;
                m_constraints.erase(c);
                break;
            }
        }
    }
    if (!current_state_on)
    {
        AddAttachmentConstraint(best_candidate);
    }

    return true;
}

void Simulation::SelectAtttachmentConstraint(AttachmentConstraint* ac)
{
    m_selected_attachment_constraint = ac;
    m_selected_attachment_constraint->Select();
}

void Simulation::UnselectAttachmentConstraint()
{
    if (m_selected_attachment_constraint)
    {
        m_selected_attachment_constraint->UnSelect();
    }
    m_selected_attachment_constraint = NULL;
}

void Simulation::AddAttachmentConstraint(unsigned int vertex_index)
{
    AttachmentConstraint* ac = new AttachmentConstraint(&m_stiffness_attachment, vertex_index, m_mesh->m_current_positions.block_vector(vertex_index));
    m_constraints.push_back(ac);
}

void Simulation::MoveSelectedAttachmentConstraintTo(const EigenVector3& target)
{
    if (m_selected_attachment_constraint)
        m_selected_attachment_constraint->SetFixedPoint(target);
}

void Simulation::clearConstraints()
{
    for (unsigned int i = 0; i < m_constraints.size(); ++i)
    {
        delete m_constraints[i];
    }
    m_constraints.clear();
}

void Simulation::setupConstraints()
{
    clearConstraints();

    switch(m_mesh->m_mesh_type)
    {
    case MESH_TYPE_CLOTH:
        // procedurally generate constraints including to attachment constraints
        {
            // generating attachment constraints.
            AddAttachmentConstraint(0);
            AddAttachmentConstraint(m_mesh->m_dim[1]*(m_mesh->m_dim[0]-1));

            // generate stretch constraints. assign a stretch constraint for each edge.
            EigenVector3 p1, p2;
            for(std::vector<Edge>::iterator e = m_mesh->m_edge_list.begin(); e != m_mesh->m_edge_list.end(); ++e)
            {
                p1 = m_mesh->m_current_positions.block_vector(e->m_v1);
                p2 = m_mesh->m_current_positions.block_vector(e->m_v2);
                SpringConstraint* c = new SpringConstraint(&m_stiffness_stretch, e->m_v1, e->m_v2, (p1-p2).norm());
                m_constraints.push_back(c);
            }

            // generate bending constraints. naive
            unsigned int i, k;
            for(i = 0; i < m_mesh->m_dim[0]; ++i)
            {
                for(k = 0; k < m_mesh->m_dim[1]; ++k)
                {
                    unsigned int index_self = m_mesh->m_dim[1] * i + k;
                    p1 = m_mesh->m_current_positions.block_vector(index_self);
                    if (i+2 < m_mesh->m_dim[0])
                    {
                        unsigned int index_row_1 = m_mesh->m_dim[1] * (i + 2) + k;
                        p2 = m_mesh->m_current_positions.block_vector(index_row_1);
                        SpringConstraint* c = new SpringConstraint(&m_stiffness_bending, index_self, index_row_1, (p1-p2).norm());
                        m_constraints.push_back(c);
                    }
                    if (k+2 < m_mesh->m_dim[1])
                    {
                        unsigned int index_column_1 = m_mesh->m_dim[1] * i + k + 2;
                        p2 = m_mesh->m_current_positions.block_vector(index_column_1);
                        SpringConstraint* c = new SpringConstraint(&m_stiffness_bending, index_self, index_column_1, (p1-p2).norm());
                        m_constraints.push_back(c);
                    }
                }
            }
        }
        break;
    case MESH_TYPE_TET:
        {
            // generate stretch constraints. assign a stretch constraint for each edge.
            EigenVector3 p1, p2;
            for(std::vector<Edge>::iterator e = m_mesh->m_edge_list.begin(); e != m_mesh->m_edge_list.end(); ++e)
            {
                p1 = m_mesh->m_current_positions.block_vector(e->m_v1);
                p2 = m_mesh->m_current_positions.block_vector(e->m_v2);
                SpringConstraint *c = new SpringConstraint(&m_stiffness_stretch, e->m_v1, e->m_v2, (p1-p2).norm());
                m_constraints.push_back(c);
            }
        }
        break;
    }
}

void Simulation::dampVelocity()
{
    if (std::abs(m_damping_coefficient) < EPSILON)
        return;

    // post-processing damping
    EigenVector3 pos_mc(0.0, 0.0, 0.0), vel_mc(0.0, 0.0, 0.0);
    unsigned int i, size;
    ScalarType denominator(0.0), mass(0.0);
    size = m_mesh->m_vertices_number;
    for(i = 0; i < size; ++i)
    {
        mass = m_mesh->m_mass_matrix.coeff(i*3, i*3);

        pos_mc += mass * m_mesh->m_current_positions.block_vector(i);
        vel_mc += mass * m_mesh->m_current_velocities.block_vector(i);
        denominator += mass;
    }
    assert(denominator != 0.0);
    pos_mc /= denominator;
    vel_mc /= denominator;

    EigenVector3 angular_momentum(0.0, 0.0, 0.0), r(0.0, 0.0, 0.0);
    EigenMatrix3 inertia, r_mat;
    inertia.setZero(); r_mat.setZero();

    for(i = 0; i < size; ++i)
    {
        mass = m_mesh->m_mass_matrix.coeff(i*3, i*3);

        r = m_mesh->m_current_positions.block_vector(i) - pos_mc;
        angular_momentum += r.cross(mass * m_mesh->m_current_velocities.block_vector(i));

        //r_mat = EigenMatrix3(0.0,  r.z, -r.y,
        //                    -r.z, 0.0,  r.x,
        //                    r.y, -r.x, 0.0);

        r_mat.coeffRef(0, 1) = r[2];
        r_mat.coeffRef(0, 2) = -r[1];
        r_mat.coeffRef(1, 0) = -r[2];
        r_mat.coeffRef(1, 2) = r[0];
        r_mat.coeffRef(2, 0) = r[1];
        r_mat.coeffRef(2, 1) = -r[0];

        inertia += r_mat * r_mat.transpose() * mass;
    }
    EigenVector3 angular_vel = inertia.inverse() * angular_momentum;

    EigenVector3 delta_v(0.0, 0.0, 0.0);
    for(i = 0; i < size; ++i)
    {
        r = m_mesh->m_current_positions.block_vector(i) - pos_mc;
        delta_v = vel_mc + angular_vel.cross(r) - m_mesh->m_current_velocities.block_vector(i);     
        m_mesh->m_current_velocities.block_vector(i) += m_damping_coefficient * delta_v;
    }
}

void Simulation::calculateInertiaY()
{
    m_inertia_y = m_mesh->m_current_positions + m_mesh->m_current_velocities * m_h;
}

void Simulation::calculateExternalForce()
{
    m_external_force.resize(m_mesh->m_system_dimension);
    m_external_force.setZero();

    // gravity
    for (unsigned int i = 0; i < m_mesh->m_vertices_number; ++i)
    {
        m_external_force[3*i+1] += -m_gravity_constant;
    }

    m_external_force = m_mesh->m_mass_matrix * m_external_force;
}

VectorX Simulation::collisionDetection(const VectorX x)
{
    // Naive implementation of collision detection
    VectorX penetration(m_mesh->m_system_dimension);
    penetration.setZero();
    EigenVector3 normal;
    ScalarType dist;

    for (unsigned int i = 0; i < m_mesh->m_vertices_number; ++i)
    {
        EigenVector3 xi = x.block_vector(i);

        if (m_scene->StaticIntersectionTest(xi, normal, dist))
        {
            penetration.block_vector(i) += (dist) * normal;
        }
    }

    return penetration;
}

void Simulation::integrateExplicitEuler()
{
    // q_{n+1} - 2q_n + q_{n-1} = h^2 * M^{-1} * force(q_{n-1})

    // inertia term 2q_n - q_{n-1} is calculated in calculateInertiaY function

    // calculate force(q_{n-1})
    VectorX position_previous = m_mesh->m_current_positions - m_mesh->m_current_velocities*m_h;
    VectorX force_previous;
    computeForces(position_previous, force_previous);

    // update q_{n+1}
    ScalarType h_square = m_h*m_h;
    VectorX pos_next = m_inertia_y + h_square*m_mesh->m_inv_mass_matrix*force_previous;
    updatePosAndVel(pos_next);
}

void Simulation::integrateExplicitSymplectic()
{
    // q_{n+1} - 2q_n + q_{n-1} = h^2 * M^{-1} * force(q_{n})

    // inertia term 2q_n - q_{n-1} is calculated in calculateInertiaY function

    // calculate force(q_{n})
    VectorX force_current;
    computeForces(m_mesh->m_current_positions, force_current);

    // update q_{n+1}
    ScalarType h_square = m_h*m_h;
    VectorX pos_next = m_inertia_y + h_square*m_mesh->m_inv_mass_matrix*force_current;
    updatePosAndVel(pos_next);
}

void Simulation::integrateImplicitBW()
{
    // (M+h^2*hessian(q_n))q_{n+1} =  M*(2q_n - q_{n-1}) + h^2*force(q_{n}) + h^2*hessian(q_n)*q_n

    // inertia term 2q_n - q_{n-1} is calculated in calculateInertiaY function

    // calculate force(q_{n})
    VectorX force_current;
    computeForces(m_mesh->m_current_positions, force_current);

    // calculate hessian(q_{n})
    SparseMatrix stiffness_matrix;
    computeStiffnessMatrix(m_mesh->m_current_positions, stiffness_matrix);

    ScalarType h_square = m_h*m_h;
    // formulate A*x = b
    VectorX b = m_mesh->m_mass_matrix*m_inertia_y + h_square * (force_current + stiffness_matrix*m_mesh->m_current_positions);
    SparseMatrix A = (m_mesh->m_mass_matrix + h_square * stiffness_matrix);
    // solve A*x = b using direct solver
    Eigen::SimplicialLDLT<SparseMatrix, Eigen::Upper> ldltSolver;
    factorizeDirectSolverLDLT(A, ldltSolver, "Baraff-Witkin method");
    VectorX pos_next = ldltSolver.solve(b);

    // update q_{n+1}
    updatePosAndVel(pos_next);
}

void Simulation::integrateOptimizationMethod()
{
    // take a initial guess
    VectorX pos_next = m_inertia_y;

    // while loop until converge or exceeds maximum iterations
    bool converge = false;

    for (unsigned int iteration_num = 0; !converge && iteration_num < m_iterations_per_frame; ++iteration_num)
    {
        switch (m_integration_method)
        {
        case INTEGRATION_GRADIENT_DESCENT:
            converge = integrateGradientDescentOneIteration(pos_next);
            break;
        case INTEGRATION_NEWTON_DESCENT:
            converge = integrateNewtonDescentOneIteration(pos_next);
            break;
        case INTEGRATION_NWETON_DESCENT_PCG:
            converge = integrateNewtonDescentWithPCGOneIteration(pos_next);
            break;
        case INTEGRATION_LOCAL_GLOBAL:
            converge = integrateLocalGlobalOneIteration(pos_next);
            break;
        }
    }

    // update q_{n+1}
    updatePosAndVel(pos_next);
}

bool Simulation::integrateGradientDescentOneIteration(VectorX& x)
{
    // evaluate gradient direction
    VectorX gradient;
    evaluateGradient(x, gradient);

    if (gradient.squaredNorm() < EPSILON)
        return true;

    // assign descent direction
    VectorX descent_dir = -gradient;

    // line search
    ScalarType step_size = lineSearch(x, gradient, descent_dir);

    // update x
    x = x + descent_dir * step_size;

    // report convergence
    if (step_size < EPSILON)
        return true;
    else
        return false;
}

bool Simulation::integrateNewtonDescentOneIteration(VectorX& x)
{
    // evaluate gradient direction
    VectorX gradient;
    evaluateGradient(x, gradient);

    if (gradient.squaredNorm() < EPSILON)
        return true;

    // evaluate hessian matrix
    SparseMatrix hessian;
    evaluateHessian(x, hessian);

    Eigen::SimplicialLDLT<SparseMatrix, Eigen::Upper> ldltSolver;
    factorizeDirectSolverLDLT(hessian, ldltSolver, "Newton Descent");
    VectorX descent_dir = -ldltSolver.solve(gradient);

    // line search
     ScalarType step_size = lineSearch(x, gradient, descent_dir);

    // update x
    x = x + descent_dir * step_size;
    // report convergence
    if (step_size < EPSILON)
        return true;
    else
        return false;
}

bool Simulation::integrateNewtonDescentWithPCGOneIteration(VectorX& x)
{
    // evaluate gradient direction
    VectorX gradient;
    evaluateGradient(x, gradient);

    if (gradient.squaredNorm() < EPSILON)
        return true;

    // evaluate hessian matrix
    SparseMatrix hessian;
    evaluateHessian(x, hessian);

    // assign descent direction using pcg solver
    VectorX descent_dir = gradient;
    Eigen::ConjugateGradient<SparseMatrix> cg;
    cg.compute(hessian);
    cg.setTolerance(LARGER_EPSILON);
    int cg_iterations = 0;
    int max_cg_iterations = 2000;
    do 
    {
        descent_dir = cg.solveWithGuess(gradient, descent_dir);
        cg_iterations += cg.iterations();
    }while(cg.info() != Eigen::Success && cg_iterations < max_cg_iterations);
    if (cg_iterations >= max_cg_iterations)
    {
        std::cout << "Warning: CG Solver does not converge in " << max_cg_iterations << " iterations. " << std::endl;
    }
    descent_dir = -descent_dir;

    // line search
    ScalarType step_size = lineSearch(x, gradient, descent_dir);

    // update x
    x = x + descent_dir * step_size;

    // report convergence
    if (step_size < EPSILON)
        return true;
    else
        return false;
}

bool Simulation::integrateLocalGlobalOneIteration(VectorX& x)
{
    prefactorize(PREFACTOR_M_PLUS_H2L);

    VectorX d;
    evaluateDVector(x, d);

    VectorX b = m_mesh->m_mass_matrix * m_inertia_y + m_h*m_h*(m_J_matrix*d+m_external_force);

    x = m_prefactored_solver[PREFACTOR_M_PLUS_H2L].solve(b);

    return false;
}

ScalarType Simulation::evaluateObjectiveFunction(const VectorX& x)
{
    ScalarType potential_term = 0.0;

    for (std::vector<Constraint*>::iterator it = m_constraints.begin(); it != m_constraints.end(); ++it)
    {
        potential_term += (*it)->EvaluatePotentialEnergy(x);
    }

    // external force
    potential_term -= x.transpose()*m_external_force;

    ScalarType inertia_term = 0.5 * (x-m_inertia_y).transpose() * m_mesh->m_mass_matrix * (x-m_inertia_y);
    ScalarType h_square = m_h*m_h;

    return inertia_term + potential_term * h_square;    
}

void Simulation::evaluateGradient(const VectorX& x, VectorX& gradient)
{
    gradient.resize(m_mesh->m_system_dimension);
    gradient.setZero();

    // springs
    for (std::vector<Constraint*>::iterator it = m_constraints.begin(); it != m_constraints.end(); ++it)
    {
        (*it)->EvaluateGradient(x, gradient);
    }

    // external forces
    gradient -= m_external_force;

    ScalarType h_square = m_h*m_h;
    gradient = m_mesh->m_mass_matrix * (x - m_inertia_y) + h_square*gradient;
}

void Simulation::evaluateHessian(const VectorX& x, SparseMatrix& hessian_matrix)
{
    hessian_matrix.resize(m_mesh->m_system_dimension, m_mesh->m_system_dimension);
    std::vector<SparseMatrixTriplet> h_triplets;
    h_triplets.clear();

    for (std::vector<Constraint*>::iterator it = m_constraints.begin(); it != m_constraints.end(); ++it)
    {
        (*it)->EvaluateHessian(x, h_triplets);
    }

    hessian_matrix.setFromTriplets(h_triplets.begin(), h_triplets.end());
    ScalarType h_square = m_h*m_h;
    hessian_matrix = m_mesh->m_mass_matrix + h_square*hessian_matrix;
}

// for local global method only
void Simulation::evaluateDVector(const VectorX& x, VectorX& d)
{
    d.resize(m_constraints.size()*3);
    d.setZero();

    for (unsigned int index = 0; index < m_constraints.size(); ++index)
    {
        m_constraints[index]->EvaluateDVector(index, x, d);
    }
}
void Simulation::evaluateJMatrix(SparseMatrix& J)
{
    J.resize(m_mesh->m_system_dimension, m_constraints.size()*3);
    std::vector<SparseMatrixTriplet> J_triplets;
    J_triplets.clear();

    for (unsigned int index = 0; index < m_constraints.size(); ++index)
    {
        m_constraints[index]->EvaluateJMatrix(index, J_triplets);
    }
    J.setFromTriplets(J_triplets.begin(), J_triplets.end());
}


ScalarType Simulation::lineSearch(const VectorX& x, const VectorX& gradient_dir, const VectorX& descent_dir)
{
    if (m_enable_line_search)
    {
        VectorX x_plus_tdx(m_mesh->m_system_dimension);
        ScalarType t = 1.0/m_ls_beta;
        ScalarType lhs, rhs;

        ScalarType currentObjectiveValue = evaluateObjectiveFunction(x);

        do
        {
            t *= m_ls_beta;
            x_plus_tdx = x + t*descent_dir;
        
            lhs = evaluateObjectiveFunction(x_plus_tdx);
            rhs = currentObjectiveValue + m_ls_alpha * t * (gradient_dir.transpose() * descent_dir)(0);

        } while (lhs >= rhs && t > EPSILON);

        if (t < EPSILON)
        {
            t = 0.0;
        }
        else
        {
            m_ls_step_size = t;
        }

        return t;
    }
    else
    {
        return m_ls_step_size;
    }
}

#pragma region implicit/explicit euler
void Simulation::computeForces(const VectorX& x, VectorX& force)
{
    VectorX gradient;

    gradient.resize(m_mesh->m_system_dimension);
    gradient.setZero();

    // springs
    for (std::vector<Constraint*>::iterator it = m_constraints.begin(); it != m_constraints.end(); ++it)
    {
        (*it)->EvaluateGradient(x, gradient);
    }

    // external forces
    gradient -= m_external_force;

    force = -gradient;
}

void Simulation::computeStiffnessMatrix(const VectorX& x, SparseMatrix& stiffness_matrix)
{
    evaluateHessian(x, stiffness_matrix);
}

#pragma endregion

#pragma region matrices and prefactorization
void Simulation::setWeightedLaplacianMatrix()
{
    m_weighted_laplacian.resize(m_mesh->m_system_dimension, m_mesh->m_system_dimension);
    std::vector<SparseMatrixTriplet> l_triplets;
    l_triplets.clear();

    for (std::vector<Constraint*>::iterator it = m_constraints.begin(); it != m_constraints.end(); ++it)
    {
        (*it)->EvaluateWeightedLaplacian(l_triplets);
    }

    m_weighted_laplacian.setFromTriplets(l_triplets.begin(), l_triplets.end());
}

void Simulation::setJMatrix()
{
    m_J_matrix;

    evaluateJMatrix(m_J_matrix);
}

void Simulation::prefactorize(PrefactorType type)
{
    if (m_prefatorization_flag[type] == false)
    {
        SparseMatrix A;
        ScalarType h2 = m_h*m_h;

        // choose matrix
        switch (type)
        {
        case PREFACTOR_L:
            setWeightedLaplacianMatrix();
            setJMatrix();
            A = m_weighted_laplacian;
            break;
        case PREFACTOR_M_PLUS_H2L:
            setWeightedLaplacianMatrix();
            setJMatrix();
            A = m_mesh->m_mass_matrix + h2 * m_weighted_laplacian;
            break;
        }

        factorizeDirectSolverLLT(A, m_prefactored_solver[type]);
        m_prefatorization_flag[type] = true;
    }
}

#pragma endregion

#pragma region utilities
void Simulation::updatePosAndVel(const VectorX& new_pos)
{
    m_mesh->m_current_velocities = (new_pos - m_mesh->m_current_positions)/m_h;
    m_mesh->m_current_positions = new_pos;
}

void Simulation::factorizeDirectSolverLLT(const SparseMatrix& A, Eigen::SimplicialLLT<SparseMatrix, Eigen::Upper>& lltSolver, char* warning_msg)
{
    SparseMatrix A_prime = A;
    lltSolver.analyzePattern(A_prime);
    lltSolver.factorize(A_prime);
    ScalarType Regularization = 0.00001;
    bool success = true;
    while (lltSolver.info() != Eigen::Success)
    {
        Regularization *= 10;
        A_prime = A_prime + Regularization*m_mesh->m_identity_matrix;
        lltSolver.factorize(A_prime);
        success = false;
    }
    if (!success)
        std::cout << "Warning: " << warning_msg <<  " adding "<< Regularization <<" identites.(llt solver)" << std::endl;
}

void Simulation::factorizeDirectSolverLDLT(const SparseMatrix& A, Eigen::SimplicialLDLT<SparseMatrix, Eigen::Upper>& ldltSolver, char* warning_msg)
{
    SparseMatrix A_prime = A;
    ldltSolver.analyzePattern(A_prime);
    ldltSolver.factorize(A_prime);
    ScalarType Regularization = 0.00001;
    bool success = true;
    while (ldltSolver.info() != Eigen::Success)
    {
        Regularization *= 10;
        A_prime = A_prime + Regularization*m_mesh->m_identity_matrix;
        ldltSolver.factorize(A_prime);
        success = false;
    }
    if (!success)
        std::cout << "Warning: " << warning_msg <<  " adding "<< Regularization <<" identites.(ldlt solver)" << std::endl;
}

void Simulation::generateRandomVector(const unsigned int size, VectorX& x)
{
    x.resize(size);

    for (unsigned int i = 0; i < size; i++)
    {
        x(i) = ((ScalarType)(rand())/(ScalarType)(RAND_MAX+1)-0.5)*2;
    }
}

#pragma endregion
