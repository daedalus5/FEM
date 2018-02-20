#pragma once

#include "../globalincludes.h"
#include <string>
#include <vector>


const static int POS = 0;
const static int VEL = 1;
const static int ACC = 2;
const static int FOR = 3;
const static int MASS = 4;

struct State {
    int mStateId;
    float mMass;
    std::vector<Eigen::Vector3f> mComponents;
    std::vector<Eigen::Vector3f> mComponentDot;

    State() {
        mComponents.reserve(4);
        mComponentDot.reserve(3);
    }
};

class BaseIntegrator {

protected:
    std::string mName;

public:
    BaseIntegrator(std::string name);

    ~BaseIntegrator();

    virtual void integrate(float timeStep, int params, const State &currentState, State &newState) = 0;

    const std::string& name() const;

};




