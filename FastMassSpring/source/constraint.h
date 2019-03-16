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

#ifndef _CONSTRAINT_H_
#define _CONSTRAINT_H_

#include <vector>
#include <iostream>

#include "global_headers.h"
#include "math_headers.h"
#include "opengl_headers.h"
#include "primitive.h"

class Constraint
{
public:
    Constraint(ScalarType *stiffness);
    Constraint(const Constraint& other);
    virtual ~Constraint();

    virtual ScalarType  EvaluatePotentialEnergy(const VectorX& x) {std::cout << "Warning: reach <Constraint> base class virtual function." << std::endl; return 0;}
    virtual void  EvaluateGradient(const VectorX& x, VectorX& gradient) {std::cout << "Warning: reach <Constraint> base class virtual function." << std::endl;}
    virtual void  EvaluateHessian(const VectorX& x, std::vector<SparseMatrixTriplet>& hessian_triplets) {std::cout << "Warning: reach <Constraint> base class virtual function." << std::endl;}
    virtual void  EvaluateWeightedLaplacian(std::vector<SparseMatrixTriplet>& laplacian_triplets) {std::cout << "Warning: reach <Constraint> base class virtual function." << std::endl;}
    virtual void  EvaluateWeightedDiagonal(std::vector<SparseMatrixTriplet>& diagonal_triplets) {std::cout << "Warning: reach <Constraint> base class virtual function." << std::endl;}

    // for local global method
    virtual void EvaluateDVector(unsigned int index, const VectorX& x, VectorX& d){std::cout << "Warning: reach <Constraint> base class virtual function." << std::endl;}
    virtual void EvaluateJMatrix(unsigned int index, std::vector<SparseMatrixTriplet>& J_triplets) {std::cout << "Warning: reach <Constraint> base class virtual function." << std::endl;}

    inline const ScalarType& Stiffness() {return (*m_stiffness);}

protected:
    ScalarType *m_stiffness;

// for visualization and selection
public:
    virtual void Draw(const VBO& vbos) { /*do nothing*/ }
    //virtual ScalarType RayConstraintIntersection() {return false;}
};

class AttachmentConstraint : public Constraint
{
public:
    AttachmentConstraint(ScalarType *stiffness);
    AttachmentConstraint(ScalarType *stiffness, unsigned int p0, const EigenVector3& fixedpoint);
    AttachmentConstraint(const AttachmentConstraint& other);
    virtual ~AttachmentConstraint();

    virtual ScalarType  EvaluatePotentialEnergy(const VectorX& x);
    virtual void  EvaluateGradient(const VectorX& x, VectorX& gradient);
    virtual void  EvaluateHessian(const VectorX& x, std::vector<SparseMatrixTriplet>& hessian_triplets);
    virtual void  EvaluateWeightedLaplacian(std::vector<SparseMatrixTriplet>& laplacian_triplets);
    virtual void  EvaluateWeightedDiagonal(std::vector<SparseMatrixTriplet>& diagonal_triplets);

    // for local global method
    virtual void EvaluateDVector(unsigned int index, const VectorX& x, VectorX& d);
    virtual void EvaluateJMatrix(unsigned int index, std::vector<SparseMatrixTriplet>& J_triplets);

protected:
    unsigned int m_p0;
    EigenVector3 m_fixd_point;

// for visualization and selection
public:
    virtual void Draw(const VBO& vbos);
    inline void Select() {m_selected = true;}
    inline void UnSelect() {m_selected = false;}
    inline EigenVector3 GetFixedPoint() {return m_fixd_point;}
    inline void SetFixedPoint(const EigenVector3& target) {m_fixd_point = target;}
    inline unsigned int GetConstrainedVertexIndex() {return m_p0;}

private: 
    bool m_selected;
    Sphere m_attachment_constraint_body;
};

class SpringConstraint : public Constraint
{
public:
    SpringConstraint(ScalarType *stiffness);
    SpringConstraint(ScalarType *stiffness, unsigned int p1, unsigned int p2, ScalarType length);
    SpringConstraint(const SpringConstraint& other);
    virtual ~SpringConstraint();

    virtual ScalarType  EvaluatePotentialEnergy(const VectorX& x);
    virtual void  EvaluateGradient(const VectorX& x, VectorX& gradient);
    virtual void  EvaluateHessian(const VectorX& x, std::vector<SparseMatrixTriplet>& hessian_triplets);
    virtual void  EvaluateWeightedLaplacian(std::vector<SparseMatrixTriplet>& laplacian_triplets);
    virtual void  EvaluateWeightedDiagonal(std::vector<SparseMatrixTriplet>& diagonal_triplets);

    // for local global method
    virtual void EvaluateDVector(unsigned int index, const VectorX& x, VectorX& d);
    virtual void EvaluateJMatrix(unsigned int index, std::vector<SparseMatrixTriplet>& J_triplets);

protected:
    unsigned int m_p1, m_p2;
    // rest length
    ScalarType m_rest_length;
};

#endif