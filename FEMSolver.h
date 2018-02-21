#pragma once

#include "globalincludes.h"
#include "mesh/TetraMesh.h"
#include "mesh/Tetrahedron.h"
#include "integrator/ForwardEuler.h"



// values are for rubber;
template<class T, int dim>
double TetraMesh<T,dim>::k = 0.05;
template<class T, int dim>
double TetraMesh<T,dim>::nu = 0.49;

constexpr float cTimeStep = 0.001f;


template<class T, int dim>
class FEMSolver {
private:
    //ForwardEuler mIntegrator;
    TetraMesh<T,dim> *mTetraMesh;
    int mSteps;
    double mu;
    double lambda;
    ForwardEuler mIntegrator;

    void calculateMaterialConstants();    // calculates mu and lambda values for material
    void precomputeTetraConstants();      // precompute tetrahedron constant values
    void computeDs(Eigen::Matrix<T,dim,dim>& Ds,
                    const Tetrahedron<T,dim>& t);       // assembles Ds matrix
    void computeF(Eigen::Matrix<T,dim,dim>& F,
                    const Eigen::Matrix<T,dim,dim>& Ds,
                    const Tetrahedron<T,dim>& t);       // computes F matrix
    void computeR(Eigen::Matrix<T,dim,dim>& R,
                    const Eigen::Matrix<T,dim,dim>& F); // computes R matrix from F using SVD
    void computeFInvTran(Eigen::Matrix<T,dim,dim> F);   // computes det(F) * (F^-1)^T

public:
    FEMSolver(int steps);

    void initializeMesh();
    void cookMyJello();
};

template<class T, int dim>
FEMSolver<T,dim>::FEMSolver(int steps) : mSteps(steps), mIntegrator("explicit") {}

template<class T, int dim>
void FEMSolver<T,dim>::initializeMesh() {

    // Initialize mTetraMesh here
    TetraMesh<T,dim> *tetraMesh = new TetraMesh<T,dim>("objects/cube.1");
    mTetraMesh = tetraMesh;
    mTetraMesh->generateTetras();
}

template<class T, int dim>
void FEMSolver<T,dim>::cookMyJello() {

    // calculate deformation constants
    calculateMaterialConstants();
    // precompute tetrahedron constant values
    precomputeTetraConstants();

    // deformation gradient matrix
    Eigen::Matrix<T,dim,dim> F = Eigen::Matrix<T,dim,dim>::Zero(dim,dim);
    // deformed tetrahedron matrix
    Eigen::Matrix<T,dim,dim> Ds = Eigen::Matrix<T,dim,dim>::Zero(dim,dim);
    // SVD rotation matrix
    Eigen::Matrix<T,dim,dim> R = Eigen::Matrix<T,dim,dim>::Zero(dim,dim);
    // Piola stress tensor
    Eigen::Matrix<T,dim,dim> P = Eigen::Matrix<T,dim,dim>::Zero(dim,dim);
    // momentum conservation
    Eigen::Matrix<T,dim,1> force = Eigen::Matrix<T,dim,1>::Zero(dim);

    // <<<<< Time Loop BEGIN
    for(int i = 0; i < mSteps; ++i) {
        mTetraMesh->mParticles.zeroForces();
        // <<<<< force update BEGIN
        for(Tetrahedron<T,dim> t : mTetraMesh->mTetras){
            force = Eigen::Matrix<T,dim,1>::Zero(dim);
            computeDs(Ds, t);
            computeF(F, Ds, t);
            computeR(R, F);
            P = mu * (F - R) + lambda * (F.determinant() - 1) * F.determinant() * F.transpose().inverse();
            P *= t.mVolDmInv;
            for(int j = 1; j < dim + 1; ++j){
                mTetraMesh->mParticles.forces[t.mPIndices[j]] += P.col(j);
                force += P.col(j);
            }
            mTetraMesh->mParticles.forces[t.mPIndices[0]] += -force;
        }
        // <<<<< force update END

        // <<<<< Integration BEGIN
        for(Tetrahedron<T,dim> t : mTetraMesh->mTetras){

            State<T, dim> currState;
            State<T, dim> newState;

            currState.mComponents[POS] = mTetraMesh->mParticles.positions[t.mPIndices[dim]];
            currState.mComponents[VEL] = mTetraMesh->mParticles.velocities[t.mPIndices[dim]];
            currState.mComponents[FOR] = mTetraMesh->mParticles.forces[t.mPIndices[dim]] + Eigen::Matrix<T,dim,1>(0, -9.8f, 0);
            currState.mMass = mTetraMesh->mParticles.masses[t.mPIndices[dim]];

            mIntegrator.integrate(cTimeStep, 0, currState, newState);

            mTetraMesh->mParticles.positions[t.mPIndices[dim]] = newState.mComponents[POS];
            mTetraMesh->mParticles.velocities[t.mPIndices[dim]] = newState.mComponents[VEL];
        }
        // <<<<< Integration END

        // collision check. Loop through particles of mTetraMesh
            // Update velocity if needed.
            //TODO make a new function for velocity update?
        //TODO do we have dt???? -> Added a new variable called cTimeStep = 0.1

        mTetraMesh->outputFrame(i);
    }
    // <<<<< Time Loop END
}

template<class T, int dim>
void FEMSolver<T,dim>::calculateMaterialConstants(){
    mu = TetraMesh<T,dim>::k / (1. + TetraMesh<T,dim>::nu);
    lambda = TetraMesh<T,dim>::k * TetraMesh<T,dim>::nu / ((1. + TetraMesh<T,dim>::nu)*(1. - 2*TetraMesh<T,dim>::nu));
}

template<class T, int dim>
void FEMSolver<T,dim>::precomputeTetraConstants(){
    std::vector<Eigen::Matrix<T,dim,1>> positions;
    for(Tetrahedron<T,dim> t : mTetraMesh->mTetras){
        positions.clear();
        for(int i = 0; i < dim + 1; ++i){
            positions.push_back(mTetraMesh->mParticles.positions[t.mPIndices[i]]);
        }
        t.precompute(positions);
    }
}

template<class T, int dim>
void FEMSolver<T,dim>::computeDs(Eigen::Matrix<T,dim,dim>& Ds, const Tetrahedron<T,dim>& t){
    std::vector<Eigen::Matrix<T,dim,1>> x;
    for(int i = 0; i < dim + 1; ++i){
        x.push_back(mTetraMesh->mParticles.positions[t.mPIndices[i]]);
    }
    for(int i = 1; i < dim + 1; ++i){
        Ds.col(i) = x[i] - x[0];
    }
}

template<class T, int dim>
void FEMSolver<T,dim>::computeF(Eigen::Matrix<T,dim,dim>& F,
                    const Eigen::Matrix<T,dim,dim>& Ds,
                    const Tetrahedron<T,dim>& t){
    F = Ds * t.mDmInv;
}

template<class T, int dim>
void FEMSolver<T,dim>::computeR(Eigen::Matrix<T,dim,dim>& R,
                    const Eigen::Matrix<T,dim,dim>& F){
    Eigen::JacobiSVD<Eigen::Matrix<T,dim,dim>> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::Matrix<T,dim,dim> U = svd.matrixU();
    Eigen::Matrix<T,dim,dim> V = svd.matrixV(); // make sure this does need to be transposed

    if(U.determinant() < 0){
        U.col(dim - 1) = -1 * U.col(dim - 1);
    }
    if(V.determinant() < 0){
        V.col(dim - 1) = -1 * V.col(dim - 1);
    }

    R = U * V.transpose();
}

template<class T, int dim>
void FEMSolve<T,dim>::computeFInvTran(Eigen::Matrix<T,dim,dim> F){
    det = F.determinant();
    switch(dim){
        case 2:
            
            break;
        case 3:
            break;
        default: std::cout << "error: dimension must be 2 or 3" << std::endl;
    }
}