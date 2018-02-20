
#include "Spring.h"

Spring::Spring() : mDampConstant(0.f) , mSpringConstant(0.f), mRestLength(0.f), mCurrLength(0.f) {}

Spring::Spring(float kConstant, float bDampConstant, float restLen) : mDampConstant(bDampConstant), mSpringConstant(kConstant) , mRestLength(restLen), mCurrLength(0.f) {}

Spring::Spring(const Spring& spring) : mDampConstant(spring.getDampConstant()), mSpringConstant(spring.getSpringConstant()) , mRestLength(spring.getRestLength()), mCurrLength(spring.getCurrLength()) {}

float Spring::getRestLength() const {
    return mRestLength;
}

float Spring::getSpringConstant() const {
    return mSpringConstant;
}

float Spring::getDampConstant() const {
    return mDampConstant;
}

const std::vector<Eigen::Vector3f>& Spring::getEndPoints() {
    return mEndPoints;
}

const std::vector<Eigen::Vector3f>& Spring::getPointVelocities() {
    return mPtVelocities;
}

float Spring::getCurrLength() const {
    return mCurrLength;
}

const Eigen::Vector3f& Spring::getSpringForce() const {
    return mSpringForce;
}

const Eigen::Vector3f& Spring::getDampingForce() const {
    return mDampForce;
}

void Spring::setSpringConstant(float springConstant) {
    mSpringConstant = springConstant;
}

void Spring::setDampConstant(float dampConstant) {
    mDampConstant = dampConstant;
}

void Spring::setRestLength(float restLength) {
    mRestLength = restLength;
}

void Spring::setEndPoints(const Eigen::Vector3f &endPoint1, const Eigen::Vector3f &endPoint2) {
    Eigen::Vector3f point1 = Eigen::Vector3f(endPoint1);
    Eigen::Vector3f point2 = Eigen::Vector3f(endPoint2);
    mEndPoints.push_back(point1);
    mEndPoints.push_back(point2);
    mCurrLength = (point1 - point2).norm();
}

void Spring::setEndPoints(const std::vector<Eigen::Vector3f> &endPoints) {
    Eigen::Vector3f point1 = Eigen::Vector3f(endPoints[0]);
    Eigen::Vector3f point2 = Eigen::Vector3f(endPoints[1]);
    mEndPoints.push_back(point1);
    mEndPoints.push_back(point2);
    mCurrLength = (point1 - point2).norm();
}

void Spring::setPointVelocities(const Eigen::Vector3f &velocity1, const Eigen::Vector3f &velocity2) {
    Eigen::Vector3f vel1 = Eigen::Vector3f(velocity1);
    Eigen::Vector3f vel2 = Eigen::Vector3f(velocity2);
    mPtVelocities.push_back(vel1);
    mPtVelocities.push_back(vel2);
}

void Spring::setPointVelocities(const std::vector<Eigen::Vector3f> &velocities) {
    Eigen::Vector3f vel1 = Eigen::Vector3f(velocities[0]);
    Eigen::Vector3f vel2 = Eigen::Vector3f(velocities[1]);
    mPtVelocities.push_back(vel1);
    mPtVelocities.push_back(vel2);
}

void Spring::recompute() {
    float scalar =  -mSpringConstant * (mCurrLength / mRestLength - 1.f);
    Eigen::Vector3f dir = mEndPoints[0] - mEndPoints[1];
    dir.normalize();
    mSpringForce = -scalar * dir;

    Eigen::Vector3f n = mEndPoints[1] - mEndPoints[0];
    n.normalize();

    mDampForce = -mDampConstant * (n * n.transpose()) * (mPtVelocities[1] - mPtVelocities[0]);

}


