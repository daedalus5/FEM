#include "BaseIntegrator.h"

template<class T, int dim>
BaseIntegrator<T, dim>::BaseIntegrator(std::string name) : mName(name) {}

template<class T, int dim>
BaseIntegrator<T, dim>::~BaseIntegrator() {}

template<class T, int dim>
const std::string& BaseIntegrator<T, dim>::name() const{
    return mName;
}
