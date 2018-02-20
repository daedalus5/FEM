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

    void calculateConstants();      // calculates mu and lambda values for material

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
    TetraMesh<T,dim> *tetraMesh = new TetraMesh<T,dim>("hello");
    mTetraMesh = tetraMesh;

}

template<class T, int dim>
void FEMSolver<T,dim>::cookMyJello() {

    // <<<<<<< BEGIN TESTING >>>>>>>>
    Eigen::Matrix<float, 3, 1> x1 = {2, 0, 0};
    Eigen::Matrix<float, 3, 1> x2 = {0, 2, 0};
    Eigen::Matrix<float, 3, 1> x3 = {0, 0, 2};
    Eigen::Matrix<float, 3, 1> x4 = {0, 0, 0};
    std::vector<Eigen::Matrix<float, 3, 1>> x = {x1, x2, x3, x4};

    std::vector<int> indices = {0, 1, 2, 3};
    Tetrahedron<float, 3> test(indices);

    test.precompute(x);
    test.print_info();
    // <<<<<<< END TESTING >>>>>>>>

    // calculate deformation constants
    calculateConstants();

    // precompute tetrahedron constant values
    std::vector<Eigen::Matrix<T,dim,1>> positions;
    for(Tetrahedron<T,dim> t : mTetraMesh->mTetras){
        positions.clear();
        for(int i = 0; i < dim + 1; ++i){
            positions.push_back(mTetraMesh->mParticles.positions[t.mPIndices[i]]);
        }
        t.precompute(positions);
    }

    for(int i = 0; i < mSteps; ++i){
        mTetraMesh->mParticles.zeroForces();
        for(Tetrahedron<T,dim> t : mTetraMesh->mTetras){

        }
    }
    // time loop
        // zero forces
        // compute forces
    // end loop

    // time loop
        // update particles
        // write particle data
    // end loop
}

template<class T, int dim>
void FEMSolver<T,dim>::calculateConstants(){
    mu = TetraMesh<T,dim>::k / (1. + TetraMesh<T,dim>::nu);
    lambda = TetraMesh<T,dim>::k * TetraMesh<T,dim>::nu / ((1. + TetraMesh<T,dim>::nu)*(1. - 2*TetraMesh<T,dim>::nu));
}