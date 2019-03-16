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

#include "constraint.h"

//----------Constraint Class----------//
Constraint::Constraint(ScalarType *stiffness) : 
    m_stiffness(stiffness)
{
}

Constraint::Constraint(const Constraint& other) : 
    m_stiffness(other.m_stiffness)
{
}

Constraint::~Constraint()
{
}

//----------AttachmentConstraint Class----------//
AttachmentConstraint::AttachmentConstraint(ScalarType *stiffness) : 
    Constraint(stiffness)
{
    m_selected = false;
}

AttachmentConstraint::AttachmentConstraint(ScalarType *stiffness, unsigned int p0, const EigenVector3& fixedpoint) : 
    Constraint(stiffness),
    m_p0(p0),
    m_fixd_point(fixedpoint)
{
    m_selected = false;
}

AttachmentConstraint::AttachmentConstraint(const AttachmentConstraint& other) : 
    Constraint(other),
    m_p0(other.m_p0),
    m_fixd_point(other.m_fixd_point),
    m_selected(other.m_selected)
{
    
}

AttachmentConstraint::~AttachmentConstraint()
{
}

// 0.5*k*(current_length)^2
ScalarType AttachmentConstraint::EvaluatePotentialEnergy(const VectorX& x)
{
    ScalarType e_i = 0.5*(*(m_stiffness))*(x.block_vector(m_p0) - m_fixd_point).squaredNorm();

    return e_i;
}

// attachment spring gradient: k*(current_length)*current_direction
void AttachmentConstraint::EvaluateGradient(const VectorX& x, VectorX& gradient)
{
    EigenVector3 g_i = (*(m_stiffness))*(x.block_vector(m_p0) - m_fixd_point);
    gradient.block_vector(m_p0) += g_i;
}

void AttachmentConstraint::EvaluateHessian(const VectorX& x, std::vector<SparseMatrixTriplet>& hessian_triplets)
{
    ScalarType ks = *(m_stiffness);
    hessian_triplets.push_back(SparseMatrixTriplet(3*m_p0, 3*m_p0, ks));
    hessian_triplets.push_back(SparseMatrixTriplet(3*m_p0+1, 3*m_p0+1, ks));
    hessian_triplets.push_back(SparseMatrixTriplet(3*m_p0+2, 3*m_p0+2, ks));
}

void AttachmentConstraint::EvaluateWeightedLaplacian(std::vector<SparseMatrixTriplet>& laplacian_triplets)
{
    ScalarType ks = *(m_stiffness);
    laplacian_triplets.push_back(SparseMatrixTriplet(3*m_p0, 3*m_p0, ks));
    laplacian_triplets.push_back(SparseMatrixTriplet(3*m_p0+1, 3*m_p0+1, ks));
    laplacian_triplets.push_back(SparseMatrixTriplet(3*m_p0+2, 3*m_p0+2, ks));
}

void AttachmentConstraint::EvaluateWeightedDiagonal(std::vector<SparseMatrixTriplet>& diagonal_triplets)
{
    ScalarType ks = *(m_stiffness);
    diagonal_triplets.push_back(SparseMatrixTriplet(3*m_p0, 3*m_p0, ks));
    diagonal_triplets.push_back(SparseMatrixTriplet(3*m_p0+1, 3*m_p0+1, ks));
    diagonal_triplets.push_back(SparseMatrixTriplet(3*m_p0+2, 3*m_p0+2, ks));
}

void AttachmentConstraint::EvaluateDVector(unsigned int index, const VectorX& x, VectorX& d)
{
    d.block_vector(index) = m_fixd_point;
}

void AttachmentConstraint::EvaluateJMatrix(unsigned int index, std::vector<SparseMatrixTriplet>& J_triplets)
{
    ScalarType ks = *(m_stiffness);
    J_triplets.push_back(SparseMatrixTriplet(3*m_p0, 3*index, ks));
    J_triplets.push_back(SparseMatrixTriplet(3*m_p0+1, 3*index+1, ks));
    J_triplets.push_back(SparseMatrixTriplet(3*m_p0+2, 3*index+2, ks));
}

void AttachmentConstraint::Draw(const VBO& vbos)
{
    m_attachment_constraint_body.move_to(Eigen2GLM(m_fixd_point));
    if (m_selected)
        m_attachment_constraint_body.change_color(glm::vec3(0.8, 0.8, 0.2));
    else
        m_attachment_constraint_body.change_color(glm::vec3(0.8, 0.2, 0.2));
        
    m_attachment_constraint_body.Draw(vbos);
}

//----------SpringConstraint Class----------//
SpringConstraint::SpringConstraint(ScalarType *stiffness) : 
    Constraint(stiffness)
{
}

SpringConstraint::SpringConstraint(ScalarType *stiffness, unsigned int p1, unsigned int p2, ScalarType length) : 
    Constraint(stiffness),
    m_p1(p1),
    m_p2(p2),
    m_rest_length(length)
{
}

SpringConstraint::SpringConstraint(const SpringConstraint& other) : 
    Constraint(other),
    m_p1(other.m_p1),
    m_p2(other.m_p2),
    m_rest_length(other.m_rest_length)
{
}

SpringConstraint::~SpringConstraint()
{
}

// 0.5*k*(current_length - rest_length)^2
ScalarType SpringConstraint::EvaluatePotentialEnergy(const VectorX& x)
{
    EigenVector3 x_ij = x.block_vector(m_p1) - x.block_vector(m_p2);
    ScalarType length_difference_ij = x_ij.norm()-m_rest_length;
    ScalarType e_ij = 0.5*(*(m_stiffness))*length_difference_ij*length_difference_ij;

    return e_ij;
}

// sping gradient: k*(current_length-rest_length)*current_direction;
void SpringConstraint::EvaluateGradient(const VectorX& x, VectorX& gradient)
{
    EigenVector3 x_ij = x.block_vector(m_p1) - x.block_vector(m_p2);
    EigenVector3 g_ij = (*(m_stiffness))*(x_ij.norm()-m_rest_length)*x_ij.normalized();
    gradient.block_vector(m_p1) += g_ij;
    gradient.block_vector(m_p2) -= g_ij;
}

void SpringConstraint::EvaluateHessian(const VectorX& x, std::vector<SparseMatrixTriplet>& hessian_triplets)
{
    EigenVector3 x_ij = x.block_vector(m_p1) - x.block_vector(m_p2);
    ScalarType l_ij = x_ij.norm();
    ScalarType l0 = m_rest_length;
    ScalarType ks = *(m_stiffness);

    EigenMatrix3 k = ks * (EigenMatrix3::Identity() - l0/l_ij*(EigenMatrix3::Identity() - (x_ij*x_ij.transpose())/(l_ij*l_ij)));

    for (int row = 0; row < 3; row ++)
    {
        for (int col = 0; col < 3; col ++)
        {
            ScalarType val = k(row,col);
            //Update the global hessian matrix
            hessian_triplets.push_back(SparseMatrixTriplet(3*m_p1+row, 3*m_p1+col, val));
            hessian_triplets.push_back(SparseMatrixTriplet(3*m_p1+row, 3*m_p2+col, -val));
            hessian_triplets.push_back(SparseMatrixTriplet(3*m_p2+row, 3*m_p1+col, -val));
            hessian_triplets.push_back(SparseMatrixTriplet(3*m_p2+row, 3*m_p2+col, val));
        }
    }
}

void SpringConstraint::EvaluateWeightedLaplacian(std::vector<SparseMatrixTriplet>& laplacian_triplets)
{
    ScalarType ks = *(m_stiffness);
    // block 1 1
    laplacian_triplets.push_back(SparseMatrixTriplet(3*m_p1+0, 3*m_p1+0, ks));
    laplacian_triplets.push_back(SparseMatrixTriplet(3*m_p1+1, 3*m_p1+1, ks));
    laplacian_triplets.push_back(SparseMatrixTriplet(3*m_p1+2, 3*m_p1+2, ks));
    // block 1 2
    laplacian_triplets.push_back(SparseMatrixTriplet(3*m_p1+0, 3*m_p2+0, -ks));
    laplacian_triplets.push_back(SparseMatrixTriplet(3*m_p1+1, 3*m_p2+1, -ks));
    laplacian_triplets.push_back(SparseMatrixTriplet(3*m_p1+2, 3*m_p2+2, -ks));
    // block 2 1
    laplacian_triplets.push_back(SparseMatrixTriplet(3*m_p2+0, 3*m_p1+0, -ks));
    laplacian_triplets.push_back(SparseMatrixTriplet(3*m_p2+1, 3*m_p1+1, -ks));
    laplacian_triplets.push_back(SparseMatrixTriplet(3*m_p2+2, 3*m_p1+2, -ks));
    // block 2 2
    laplacian_triplets.push_back(SparseMatrixTriplet(3*m_p2+0, 3*m_p2+0, ks));
    laplacian_triplets.push_back(SparseMatrixTriplet(3*m_p2+1, 3*m_p2+1, ks));
    laplacian_triplets.push_back(SparseMatrixTriplet(3*m_p2+2, 3*m_p2+2, ks));
}

void SpringConstraint::EvaluateWeightedDiagonal(std::vector<SparseMatrixTriplet>& diagonal_triplets)
{
    ScalarType ks = *(m_stiffness);
    // block 1 1
    diagonal_triplets.push_back(SparseMatrixTriplet(3*m_p1+0, 3*m_p1+0, ks));
    diagonal_triplets.push_back(SparseMatrixTriplet(3*m_p1+1, 3*m_p1+1, ks));
    diagonal_triplets.push_back(SparseMatrixTriplet(3*m_p1+2, 3*m_p1+2, ks));
    // block 2 2
    diagonal_triplets.push_back(SparseMatrixTriplet(3*m_p2+0, 3*m_p2+0, ks));
    diagonal_triplets.push_back(SparseMatrixTriplet(3*m_p2+1, 3*m_p2+1, ks));
    diagonal_triplets.push_back(SparseMatrixTriplet(3*m_p2+2, 3*m_p2+2, ks));
}

void SpringConstraint::EvaluateDVector(unsigned int index, const VectorX& x, VectorX& d)
{
    EigenVector3 x_ij = x.block_vector(m_p1) - x.block_vector(m_p2);
    EigenVector3 di = x_ij.normalized() * m_rest_length;

    d.block_vector(index) = di;
}

void SpringConstraint::EvaluateJMatrix(unsigned int index, std::vector<SparseMatrixTriplet>& J_triplets)
{
    ScalarType ks = *(m_stiffness);
    
    // block 1 1
    J_triplets.push_back(SparseMatrixTriplet(3*m_p1+0, 3*index+0, ks));
    J_triplets.push_back(SparseMatrixTriplet(3*m_p1+1, 3*index+1, ks));
    J_triplets.push_back(SparseMatrixTriplet(3*m_p1+2, 3*index+2, ks));
    // block 2 2
    J_triplets.push_back(SparseMatrixTriplet(3*m_p2+0, 3*index+0, -ks));
    J_triplets.push_back(SparseMatrixTriplet(3*m_p2+1, 3*index+1, -ks));
    J_triplets.push_back(SparseMatrixTriplet(3*m_p2+2, 3*index+2, -ks));
}
