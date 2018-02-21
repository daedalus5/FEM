//base copied from CIS 561

#pragma once


//Geometry is an abstract class since it contains a pure virtual function (i.e. a virtual function that is set to 0)
class Shape
{
public:
    Shape(): transform()
    {}

    virtual ~Shape(){}
    // Find the intersection of the ray with this Shape
    virtual bool Intersect(const Ray &ray, Intersection *isect) const = 0;
    // A helper function for computing the world-space normal, tangent, and bitangent at local-space point P

    Transform transform;
};
