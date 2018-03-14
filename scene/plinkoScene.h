#pragma once

#include "shape.h"
#include "squareplane.h"
#include "sphere.h"

template<class T, int dim>
class PlinkoScene : public Scene<T,dim>{
public:
    PlinkoScene();
    virtual ~PlinkoScene();
};

template<class T, int dim>
PlinkoScene<T, dim>::~PlinkoScene() {}

template<class T, int dim>
PlinkoScene<T, dim>::PlinkoScene() {

    // Create a ground plane
    Shape<T, dim>* ground = new SquarePlane<T, dim>();
    Eigen::Matrix<T,dim,1> gCenter = Eigen::Matrix<T,dim,1>(0.f, -5.f, 0.f);
    ground->setCenter(gCenter);
    this->shapes.push_back(ground);

    // Create sphere A
    Shape<T, dim>* sphereA = new Sphere<T, dim>();
    Eigen::Matrix<T,dim,1> aCenter = Eigen::Matrix<T,dim,1>(-0.2, -1.5, -0.5);
    sphereA->setCenter(aCenter);
    this->shapes.push_back(sphereA);

    // Create sphere B
    Shape<T, dim>* sphereB = new Sphere<T, dim>();
    Eigen::Matrix<T,dim,1> bCenter = Eigen::Matrix<T,dim,1>(3.0, -3.0f, -0.5);
    sphereB->setCenter(bCenter);
    this->shapes.push_back(sphereB);

}
