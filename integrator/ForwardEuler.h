#pragma once

#include "BaseIntegrator.h"

class ForwardEuler : public BaseIntegrator {

public:
    ForwardEuler(std::string name);

    ~ForwardEuler();

    virtual void integrate(float timeStep, int params, const State &currentState, State &newState);

};
