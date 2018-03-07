#pragma once

#include "shape.h"
#include "squareplane.h"
#include "sphere.h"

template<class T, int dim>
class DefaultScene : public Scene<T,dim>{
public:
    DefaultScene();
    virtual ~DefaultScene();
};

template<class T, int dim>
DefaultScene<T, dim>::~DefaultScene() {
    //this->~Shape<T, dim>();
}

template<class T, int dim>
DefaultScene<T, dim>::DefaultScene() {

    // Create a ground plane
    Shape<T, dim>* ground = new SquarePlane<T, dim>();
    this->shapes.push_back(ground);

}
