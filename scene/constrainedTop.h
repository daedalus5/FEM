#pragma once

#include "shape.h"
#include "squareplane.h"
#include "sphere.h"

template<class T, int dim>
class ConstrainedTop : public Scene<T,dim>{
public:
    ConstrainedTop();
    virtual ~ConstrainedTop();
};

template<class T, int dim>
ConstrainedTop<T, dim>::~ConstrainedTop() {}

template<class T, int dim>
ConstrainedTop<T, dim>::ConstrainedTop() {

    // Create a plate plane
    Shape<T, dim>* plane = new SquarePlane<T, dim>();
    Eigen::Matrix<T,dim,1> plate = Eigen::Matrix<T,dim,1>(0.f, 1.0f, 0.f);
    plane->setCenter(plate);
    this->shapes.push_back(plane);
}