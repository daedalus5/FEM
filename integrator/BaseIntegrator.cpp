#include "BaseIntegrator.h"

BaseIntegrator::BaseIntegrator(std::string name) : mName(name) {}

BaseIntegrator::~BaseIntegrator() {}

const std::string& BaseIntegrator::name() const{
    return mName;
}
