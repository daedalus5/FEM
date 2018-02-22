#pragma once

#include "globalincludes.h"
#include "mesh/TetraMesh.h"
#include "mesh/Tetrahedron.h"
#include "integrator/ForwardEuler.h"
#include "scene/squareplane.h"



// values are for rubber;
template<class T, int dim>
double TetraMesh<T,dim>::k = 100.0;
template<class T, int dim>
double TetraMesh<T,dim>::nu = 0.1;

constexpr float cTimeStep = 0.0001f;


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
    void computeJFinvT(Eigen::Matrix<T,dim,dim>& JFinvT,
                    const Eigen::Matrix<T,dim,dim>& F); // computes det(F) * (F^-1)^T

public:
    FEMSolver(int steps);
    ~FEMSolver();

    void initializeMesh();
    void cookMyJello();
};

template<class T, int dim>
FEMSolver<T,dim>::FEMSolver(int steps) : mTetraMesh(nullptr), mSteps(steps), mu(0.0), lambda(0.0), mIntegrator("explicit") {}

template<class T, int dim>
FEMSolver<T,dim>::~FEMSolver(){
    //delete mTetraMesh;
}

template<class T, int dim>
void FEMSolver<T,dim>::initializeMesh() {

    // Initialize mTetraMesh here
    mTetraMesh = new TetraMesh<T,dim>("objects/cube.1");
    mTetraMesh->generateTetras();
}

template<class T, int dim>
void FEMSolver<T,dim>::cookMyJello() {

    // Create a ground plane
    SquarePlane<T, dim> ground = SquarePlane<T, dim>();

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
    // Force matrix
    Eigen::Matrix<T,dim,dim> G = Eigen::Matrix<T,dim,dim>::Zero(dim,dim);
    // momentum conservation
    Eigen::Matrix<T,dim,1> force = Eigen::Matrix<T,dim,1>::Zero(dim);
    // det(F) * (F^-1)^T term
    Eigen::Matrix<T,dim,dim> JFinvT = Eigen::Matrix<T,dim,dim>::Zero(dim,dim);

    // <<<<< Time Loop BEGIN
    for(int i = 0; i < mSteps; ++i) {
        mTetraMesh->mParticles.zeroForces();
        // <<<<< force update BEGIN
        for(Tetrahedron<T,dim> &t : mTetraMesh->mTetras){
            force = Eigen::Matrix<T,dim,1>::Zero(dim);
            computeDs(Ds, t);
            computeF(F, Ds, t);
            computeR(R, F);
            computeJFinvT(JFinvT, F);
            P = mu * (F - R) + lambda * (F.determinant() - 1) * JFinvT;
            G = P * t.mVolDmInvT;
            // if(i == 600){
            //     std::cout << P << std::endl;
            // }
            for(int j = 1; j < dim + 1; ++j){
                mTetraMesh->mParticles.forces[t.mPIndices[j]] += G.col(j - 1);
                force += G.col(j - 1);
            }
            mTetraMesh->mParticles.forces[t.mPIndices[0]] += -force;
            if(i > 500){
                //std::cout << force[1] << std::endl;
            }
        }
        // <<<<< force update END

        // <<<<< Integration BEGIN

        int size = mTetraMesh->mParticles.positions.size();
        std::vector<Eigen::Matrix<T, dim, 1>> past_pos = mTetraMesh->mParticles.positions;

        // if(i == 600){
        //     std::cout << mTetraMesh->mParticles.forces[5] << std::endl;
        // }

        for(int j = 0; j < size; ++j) {

            State<T, dim> currState;
            State<T, dim> newState;

            currState.mComponents[POS] = mTetraMesh->mParticles.positions[j];
            currState.mComponents[VEL] = mTetraMesh->mParticles.velocities[j];
            currState.mMass = mTetraMesh->mParticles.masses[j];
            currState.mComponents[FOR] = mTetraMesh->mParticles.forces[j] / currState.mMass + Eigen::Matrix<T,dim,1>(0, -9.8f, 0);

            mIntegrator.integrate(cTimeStep, 0, currState, newState);

            mTetraMesh->mParticles.positions[j] = newState.mComponents[POS];
            mTetraMesh->mParticles.velocities[j] = newState.mComponents[VEL];

        }

        // <<<<< Integration END

        // collision check. Loop through particles of mTetraMesh
        Eigen::Matrix<T, dim, 1> temp_pos;
        for (int j = 0; j < size; ++j) {
            if (ground.checkCollisions(mTetraMesh->mParticles.positions[j], temp_pos)) {
                // Update velocity if needed.
                mTetraMesh->mParticles.velocities[j] = (temp_pos - past_pos[j]) / cTimeStep;
                mTetraMesh->mParticles.positions[j] = temp_pos;
            }
        }

        mTetraMesh->outputFrame(i);
    }
    // <<<<< Time Loop END
}

template<class T, int dim>
void FEMSolver<T,dim>::calculateMaterialConstants(){
    mu = TetraMesh<T,dim>::k / (2 * (1. + TetraMesh<T,dim>::nu));
    lambda = (TetraMesh<T,dim>::k * TetraMesh<T,dim>::nu) / ((1. + TetraMesh<T,dim>::nu)*(1. - 2*TetraMesh<T,dim>::nu));
}

template<class T, int dim>
void FEMSolver<T,dim>::precomputeTetraConstants(){
    std::vector<Eigen::Matrix<T,dim,1>> positions;
    for(Tetrahedron<T,dim> &t : mTetraMesh->mTetras){
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
        Ds.col(i - 1) = x[i] - x[0];
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
void FEMSolver<T,dim>::computeJFinvT(Eigen::Matrix<T,dim,dim>& JFinvT, const Eigen::Matrix<T,dim,dim>& F){
    switch(dim){
        case 2:
            JFinvT(0,0) = F(1,1);
            JFinvT(0,1) = -F(0,1);
            JFinvT(1,0) = -F(1,0);
            JFinvT(1,1) = F(0,0);
            break;
        case 3:
            JFinvT(0,0) = F(1,1)*F(2,2) - F(1,2)*F(2,1);
            JFinvT(0,1) = F(0,2)*F(2,1) - F(0,1)*F(2,2);
            JFinvT(0,2) = F(0,1)*F(1,2) - F(0,2)*F(1,1);
            JFinvT(1,0) = F(1,2)*F(2,0) - F(1,0)*F(2,2);
            JFinvT(1,1) = F(0,0)*F(2,2) - F(0,2)*F(2,0);
            JFinvT(1,2) = F(0,2)*F(1,0) - F(0,0)*F(1,2);
            JFinvT(2,0) = F(1,0)*F(2,1) - F(1,1)*F(2,0);
            JFinvT(2,1) = F(0,1)*F(2,0) - F(0,0)*F(2,1);
            JFinvT(2,2) = F(0,0)*F(1,1) - F(0,1)*F(1,0);
            break;
        default: std::cout << "error: dimension must be 2 or 3" << std::endl;
    }
    JFinvT.transpose();
}