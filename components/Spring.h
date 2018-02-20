
#pragma once

#include <vector>
#include "Dense"

class Spring {

private:

    float mDampConstant;

    float mSpringConstant;

    float mRestLength;

    float mCurrLength;

    std::vector<Eigen::Vector3f> mEndPoints;

    std::vector<Eigen::Vector3f> mPtVelocities;

    Eigen::Vector3f mSpringForce;

    Eigen::Vector3f mDampForce;

public:

    Spring();

    Spring(float kConstant, float bDampConstant, float restLen);

    Spring(const Spring& spring);

    float getRestLength() const;

    float getSpringConstant() const;

    float getDampConstant() const;

    const std::vector<Eigen::Vector3f>& getEndPoints();

    const std::vector<Eigen::Vector3f>& getPointVelocities();

    float getCurrLength() const;

    const Eigen::Vector3f& getSpringForce() const;

    const Eigen::Vector3f& getDampingForce() const;

    void setSpringConstant(float springConstant);

    void setDampConstant(float dampConstant);

    void setRestLength(float restLength);

    void setEndPoints(const Eigen::Vector3f &endPoint1, const Eigen::Vector3f &endPoint2);

    void setEndPoints(const std::vector<Eigen::Vector3f> &endPoints);

    void setPointVelocities(const Eigen::Vector3f &velocity1, const Eigen::Vector3f &velocity2);

    void setPointVelocities(const std::vector<Eigen::Vector3f> &velocities);

    void recompute();

};

