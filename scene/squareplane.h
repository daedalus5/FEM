// base copied from CIS 561

#pragma  once
#include <shape.h>

//A SquarePlane is assumed to have a radius of 1 and a center of <0,0,0>.
//These attributes can be altered by applying a transformation matrix to the SquarePlane.
class SquarePlane : public Shape
{
public:
    virtual bool Intersect(const Ray &ray, Intersection *isect) const;
};
