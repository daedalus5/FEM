
#include "ForwardEuler.h"

ForwardEuler::ForwardEuler(std::string name) : BaseIntegrator(name) {}


ForwardEuler::~ForwardEuler() {}


void ForwardEuler::integrate(float timeStep, int params, const State &currentState, State &newState) {

    if(currentState.mComponents.size() == 0) {
        return;
    }

    newState.mComponentDot[POS] = currentState.mComponents[VEL];
    newState.mComponentDot[VEL] = currentState.mComponents[ACC];
    newState.mComponentDot[ACC] = currentState.mComponents[FOR] / currentState.mMass;

    newState.mComponents[POS] = currentState.mComponents[POS] + newState.mComponentDot[POS] * timeStep;
    newState.mComponents[VEL] = currentState.mComponents[POS] + newState.mComponentDot[POS] * timeStep;
    newState.mComponents[ACC] = currentState.mComponents[POS] + newState.mComponentDot[POS] * timeStep;

    newState.mMass = currentState.mMass;

}