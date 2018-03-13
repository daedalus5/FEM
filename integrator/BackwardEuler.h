#pragma once

#include "BaseIntegrator.h"

class BackwardEuler : public BaseIntegrator {

public:
    BackwardEuler(std::string name);

    ~BackwardEuler();

    virtual void integrate(double timeStep, int params, const State<T, dim> &currentState, State<T, dim> &newState);

};