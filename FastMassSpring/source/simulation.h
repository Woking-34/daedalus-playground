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

#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <vector>

#include "global_headers.h"
#include "anttweakbar_wrapper.h"
#include "mesh.h"
#include "constraint.h"
#include "scene.h"

class Mesh;
class AntTweakBarWrapper;

typedef enum
{
    PREFACTOR_L,
    PREFACTOR_M_PLUS_H2L,
    PREFACTOR_TOTAL_NUM

} PrefactorType;

typedef enum
{
    INTEGRATION_LOCAL_GLOBAL,
    INTEGRATION_EXPLICIT_EULER,
    INTEGRATION_EXPLICIT_SYMPLECTIC,
    INTEGRATION_IMPLICIT_EULER_BARAFF_WITKIN,
    INTEGRATION_GRADIENT_DESCENT,
    INTEGRATION_NEWTON_DESCENT,
    INTEGRATION_NWETON_DESCENT_PCG,
    INTEGRATION_TOTAL_NUM

} IntegrationMethod;

class Simulation
{
    friend class AntTweakBarWrapper;

public:
    Simulation();
    virtual ~Simulation();

    void Reset();
    void Update();
    void DrawConstraints(const VBO& vbos);

    // select/unselect/move attachment constratins
    ScalarType TryToSelectAttachmentConstraint(const EigenVector3& p0, const EigenVector3& dir); // return ray_projection_plane_distance if hit; return -1 otherwise.
    bool TryToToggleAttachmentConstraint(const EigenVector3& p0, const EigenVector3& dir); // true if hit some vertex/constraint
    void SelectAtttachmentConstraint(AttachmentConstraint* ac);
    void UnselectAttachmentConstraint();
    void AddAttachmentConstraint(unsigned int vertex_index); // add one attachment constraint at vertex_index
    void MoveSelectedAttachmentConstraintTo(const EigenVector3& target); // move selected attachement constraint to target

    // inline functions
    inline void SetReprefactorFlag() 
    {
        for (unsigned int i = 0; i < PREFACTOR_TOTAL_NUM; ++i)
        {
            m_prefatorization_flag[i] = false;
        }
    }
    inline void SetMesh(Mesh* mesh) {m_mesh = mesh;}
    inline void SetScene(Scene* scene) {m_scene = scene;}
    
protected:

    // simulation constants
    ScalarType m_h; // time_step

    // simulation constants
    ScalarType m_gravity_constant;
    ScalarType m_stiffness_attachment;
    ScalarType m_stiffness_stretch;
    ScalarType m_stiffness_bending;
    ScalarType m_damping_coefficient;

    // integration method
    IntegrationMethod m_integration_method;

    // key simulation components: mesh and scene
    Mesh *m_mesh;
    Scene *m_scene;
    // key simulation components: constraints
    std::vector<Constraint*> m_constraints;
    AttachmentConstraint* m_selected_attachment_constraint;

    // inertia term
    VectorX m_inertia_y;

    // external force (gravity, wind, etc...)
    VectorX m_external_force;

    // for optimization method, number of iterations
    unsigned int m_iterations_per_frame;

    // line search for gradient descent and newton's method
    bool m_enable_line_search;
    ScalarType m_ls_alpha;
    ScalarType m_ls_beta;
    ScalarType m_ls_step_size;

    // for prefactorization
    SparseMatrix m_weighted_laplacian;
    SparseMatrix m_J_matrix;
    bool m_prefatorization_flag[PREFACTOR_TOTAL_NUM];
    Eigen::SimplicialLLT<SparseMatrix, Eigen::Upper> m_prefactored_solver[PREFACTOR_TOTAL_NUM];

private:

    // main update sub-routines
    void clearConstraints(); // cleanup all constraints
    void setupConstraints(); // initialize constraints
    void dampVelocity(); // damp velocity at the end of each iteration.
    void calculateInertiaY(); // calculate the inertia term: y = current_pos + current_vel*h
    void calculateExternalForce(); // wind force is propotional to the area of triangles projected on the tangential plane
    VectorX collisionDetection(const VectorX x); // detect collision and return a vector of penetration

    void integrateExplicitEuler();
    void integrateExplicitSymplectic();
    void integrateImplicitBW();
    void integrateOptimizationMethod();

    // all those "OneIteration" functions will be called in a loop
    // x is initially passed as the initial guess of the next postion (i.e. inertia term): x = y = current_pos + current_vel*h
    // x will be changed during these subroutines in EVERY iteration
    // the final value of x will be the next_pos that we used to update all vertices.
    bool integrateGradientDescentOneIteration(VectorX& x);
    bool integrateNewtonDescentOneIteration(VectorX& x);
    bool integrateNewtonDescentWithPCGOneIteration(VectorX& x);
    bool integrateLocalGlobalOneIteration(VectorX& x); // our method

    // for explicit/implicit integration only
    void computeForces(const VectorX& x, VectorX& force);
    void computeStiffnessMatrix(const VectorX& x, SparseMatrix& stiffness_matrix);

    // evaluate energy
    ScalarType evaluateObjectiveFunction(const VectorX& x);
    // evaluate gradient
    void evaluateGradient(const VectorX& x, VectorX& gradient);
    // evaluate Hessian Matrix
    void evaluateHessian(const VectorX& x, SparseMatrix& hessian_matrix); 

    // for local global method only
    void evaluateDVector(const VectorX& x, VectorX& d); // d-vector will be evaluate every iteration
    void evaluateJMatrix(SparseMatrix& J); // J matrix is only dependent on connectivity and stiffness, so it is evaluated only once.

    // line search (for newton's method and gradient descent)
    ScalarType lineSearch(const VectorX& x, const VectorX& gradient_dir, const VectorX& descent_dir);

    // matrices and prefactorizations
    void setWeightedLaplacianMatrix();
    void setJMatrix();
    void prefactorize(PrefactorType type);

    // utility functions
    void updatePosAndVel(const VectorX& new_pos); // current_vel = (next_pos-current_pos)/h; current_pos = next_pos; 
    void factorizeDirectSolverLLT(const SparseMatrix& A, Eigen::SimplicialLLT<SparseMatrix, Eigen::Upper>& lltSolver, char* warning_msg = ""); // factorize matrix A using LLT decomposition
    void factorizeDirectSolverLDLT(const SparseMatrix& A, Eigen::SimplicialLDLT<SparseMatrix, Eigen::Upper>& ldltSolver, char* warning_msg = ""); // factorize matrix A using LDLT decomposition
    void generateRandomVector(const unsigned int size, VectorX& x); // generate random vector varing from [-1 1].
};

#endif