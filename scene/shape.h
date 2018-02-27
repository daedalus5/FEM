//base copied from CIS 561

#pragma once

#include "../mesh/Particles.h"

template<class T, int dim>
class Shape
{
public:
    Shape() {}

    virtual ~Shape(){}
    virtual bool checkCollisions(Eigen::Matrix<T, dim, 1> &pos, Eigen::Matrix<T, dim, 1> &out_pos) const = 0;
    virtual void setCenter(Eigen::Matrix<T, dim, 1> &n_cen);

protected:
    Eigen::Matrix<T, dim, 1> center; //Default initialize to zero
    
};


template<class T, int dim>
void Shape<T, dim>::setCenter(Eigen::Matrix<T, dim, 1> &n_cen) {
    center = n_cen;
}

