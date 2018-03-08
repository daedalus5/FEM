
#include "BackwardEuler.h"

BackwardEuler::BackwardEuler(std::string name) : BaseIntegrator(name) {}

BackwardEuler::~BackwardEuler() {}

void BackwardEuler::integrate(float timeStep, int params, const State<T, dim> &currentState, State<T, dim> &newState) {


    float mass = currentState.mMass;
    Eigen::Matrix<T, dim, 1> velocity = currentState.mComponents[VEL];
    Eigen::Matrix<T, dim, 1> force = currentState.mComponents[FOR];

    //Eigen::Matrix<T, dim, 1> B = mass * velocity / timeStep + force;
    Eigen::Matrix<T, dim, dim> A = Eigen::Matrix<T,dim,dim>::Zero(dim,dim);

    float temp = mass / (timeStep * timeStep);

    A(0, 0) = temp;
    A(1, 1) = temp;
    A(2, 2) = temp;

    Eigen::Matrix<T, dim, dim> J; // TO be filled in

    Eigen::Matrix<T, dim, 1> dx; // To be computed from Ax = B

    newState.mComponents[VEL] = dx / timeStep;

    newState.mComponents[POS] = currentState.mComponents[POS] + dx;

    newState.mMass = currentState.mMass;

}