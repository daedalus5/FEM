// base copied from CIS 561

#pragma  once
#include "shape.h"

//A SquarePlane is assumed to have a center of <0,0,0> and is horizontal.
template<class T, int dim>
class SquarePlane : public Shape<T, dim>
{
    public:
        virtual bool checkCollisions(Eigen::Matrix<T, dim, 1> &pos, Eigen::Matrix<T, dim, 1> &out_pos) const override;

    private:
        float length_half = 50.f;
        Eigen::Matrix<T, dim, 1> center; //Default initialize to zero
};
