#pragma once

#include "shape.h"
#include "squareplane.h"
#include "sphere.h"

template<class T, int dim>
class BulldozeScene : public Scene<T,dim>{
public:
    BulldozeScene();
    virtual ~BulldozeScene();
};

template<class T, int dim>
BulldozeScene<T, dim>::~BulldozeScene() {
    //this->~Shape<T, dim>();
}

template<class T, int dim>
BulldozeScene<T, dim>::BulldozeScene() {

    // Create a ground plane
    Shape<T, dim>* ground = new SquarePlane<T, dim>();
    Eigen::Matrix<T,dim,1> gCenter = Eigen::Matrix<T,dim,1>(0.f, -1.5f, 0.f);
    ground->setCenter(gCenter);
    this->shapes.push_back(ground);

    // Create sphere moving
    Shape<T, dim>* sphereA = new Sphere<T, dim>("bulldoze_output");
    Eigen::Matrix<T,dim,1> aCenter = Eigen::Matrix<T,dim,1>(-1.5, -0.5f, -1.5);
    sphereA->setCenter(aCenter);
    Eigen::Matrix<T,dim,1> aVel = Eigen::Matrix<T,dim,1>(0.008, 0.0, 0.008);
    sphereA->setVelocity(aVel);
    this->shapes.push_back(sphereA);

}
