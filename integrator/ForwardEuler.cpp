
#include "ForwardEuler.h"

template<class T, int dim>
ForwardEuler<T, dim>::ForwardEuler(std::string name) : BaseIntegrator<T, dim>(name) {}

template<class T, int dim>
ForwardEuler<T, dim>::~ForwardEuler() {}

template<class T, int dim>
void ForwardEuler<T, dim>::integrate(float timeStep, int params, const State<T, dim> &currentState, State<T, dim> &newState) {

    if(currentState.mComponents.size() == 0) {
        return;
    }

    // Derivative of position is the current velocity of particle
    newState.mComponentDot[POS] = currentState.mComponents[VEL];

    // Derivative of velocity is = F/m (from Newtons second law)
    newState.mComponentDot[VEL] = currentState.mComponents[FOR];

    // Calculate the new positions and velocity
    newState.mComponents[POS] = currentState.mComponents[POS] + newState.mComponentDot[POS] * timeStep;
    newState.mComponents[VEL] = currentState.mComponents[VEL] + newState.mComponentDot[VEL] * timeStep;

    newState.mMass = currentState.mMass;
}
