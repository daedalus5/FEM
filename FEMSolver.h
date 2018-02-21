#pragma once

#include "mesh/TetraMesh.h"
#include "mesh/Tetrahedron.h"

using T = float;
const int dim = 3;

// values are for rubber;
template<class T, int dim>
double TetraMesh<T,dim>::k = 0.05;
template<class T, int dim>
double TetraMesh<T,dim>::nu = 0.49;


template<class T, int dim>
class FEMSolver {
private:
    //ForwardEuler mIntegrator;
    TetraMesh<T,dim> *mTetraMesh;
    int mSteps;
    double mu;
    double lambda;

    void calculateMaterialConstants();    // calculates mu and lambda values for material
    void precomputeTetraConstants();      // precompute tetrahedron constant values
    void computeDs(Eigen::Matrix<T,dim,dim>& Ds,
                    const Tetrahedron<T,dim>& t);       // assembles Ds matrix
    void computeF(Eigen::Matrix<T,dim,dim>& F,
                    const Eigen::Matrix<T,dim,dim>& Ds,
                    const Tetrahedron<T,dim>& t);       // computes F matrix
    void computeR(Eigen::Matrix<T,dim,dim>& R,
                    const Eigen::Matrix<T,dim,dim>& F); // computes R matrix from F using SVD

public:
    FEMSolver(int steps);

    void initializeMesh();
    void cookMyJello();
};

template<class T, int dim>
FEMSolver<T,dim>::FEMSolver(int steps) : mSteps(steps) {}

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

    for(int i = 0; i < mSteps; ++i){
        mTetraMesh->mParticles.zeroForces();
        for(Tetrahedron<T,dim> t : mTetraMesh->mTetras){
            computeDs(Ds, t);
            computeF(F, Ds, t);
            computeR(R, F);
        }
    }

      // NOTE: mTetraMesh->outputFrame(frame#);
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
    for(int i = 0; i < dim; ++i){
        x.push_back(mTetraMesh->mParticles.positions[t.mPIndices[i]]);
    }
    for(int i = 0; i < dim; ++i){
        Ds.col(i) = x[i] - x.back();
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
