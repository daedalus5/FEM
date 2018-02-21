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
};
