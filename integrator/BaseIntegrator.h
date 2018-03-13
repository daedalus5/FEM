#pragma once

#include "../globalincludes.h"
#include <string>
#include <vector>

const static int POS = 0;
const static int VEL = 1;
const static int ACC = 2;
const static int FOR = 3;
const static int MASS = 4;

template <class T, int dim>
struct State {
    int mStateId;
    T mMass;
    std::vector<Eigen::Matrix<T, dim, 1>> mComponents;
    std::vector<Eigen::Matrix<T, dim, 1>> mComponentDot;

    State() {
        mComponents.resize(4);
        mComponentDot.resize(3);
    }
};

template<class T, int dim>
class BaseIntegrator {

protected:
    std::string mName;

public:
    BaseIntegrator(std::string name);

    ~BaseIntegrator();

    virtual void integrate(float timeStep, int params, const State<T, dim> &currentState, State<T, dim> &newState) = 0;

    const std::string& name() const;

};


template<class T, int dim>
BaseIntegrator<T, dim>::BaseIntegrator(std::string name) : mName(name) {}

template<class T, int dim>
BaseIntegrator<T, dim>::~BaseIntegrator() {}

template<class T, int dim>
const std::string& BaseIntegrator<T, dim>::name() const{
    return mName;
}
