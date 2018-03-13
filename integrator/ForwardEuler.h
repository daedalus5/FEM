#pragma once

#include "BaseIntegrator.h"

class ForwardEuler : public BaseIntegrator {

public:
    ForwardEuler(std::string name);

    ~ForwardEuler();

    virtual void integrate(double timeStep, int params, const State<T, dim> &currentState, State<T, dim> &newState);

};
