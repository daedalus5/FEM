// base copied from CIS 561

#pragma  once
#include "shape.h"

//A SquarePlane is assumed to have a center of <0,0,0> and is horizontal.
template<class T, int dim>
class Sphere : public Shape<T, dim>
{
    public:
        Sphere() : Shape<T, dim>() {};
        Sphere(std::string file);
        bool checkCollisions(const Eigen::Matrix<T, dim, 1> &pos, Eigen::Matrix<T, dim, 1> &out_pos) const override;

    private:
        float radius = 1.0f;
};


template<class T, int dim>
Sphere<T, dim>::Sphere(std::string file) : Shape<T, dim>(file) {
    //center[1] = -1.0;
    this->center[0] = 0.3;
    this->center[1] = -1.0;
    this->center[2] = 0.5;
}

template<class T, int dim>
bool Sphere<T, dim>::checkCollisions(const Eigen::Matrix<T, dim, 1> &pos, Eigen::Matrix<T, dim, 1> &out_pos) const
{

    //if (dim == 3) {
        Eigen::Matrix<T, dim, 1> dist = pos - this->center;
        if (dist.norm() < radius) {
        //if (((pos[0] - this->center[0]) * (pos[0] - this->center[0])) + 
        //    ((pos[1] - this->center[1]) * (pos[1] - this->center[1])) + 
        //    ((pos[2] - this->center[2]) * (pos[2] - this->center[2])) - 
        //    (radius * radius) 
        //    < 0) {

            // Collision in 3D
            Eigen::Matrix<T, dim, 1> normal = pos - this->center;
            normal = normal /  normal.norm();
            normal = normal * radius;

            out_pos =  normal + this->center;
            
            //std::cout << this->center[0] << " " << this->center[1] << " " << this->center[2] << std::endl;

            return true;
        }
    //}
/*
    else if (dim == 2) {
        if (((pos[0] - this->center[0]) * (pos[0] - this->center[0])) + 
            ((pos[1] - this->center[1]) * (pos[1] - this->center[1])) + 
            (radius * radius) 
            < 0) {

            // Coliision in 2D
            Eigen::Matrix<T, dim, 1> normal = pos - this->center;
            normal = normal /  normal.norm();
            normal = normal * radius;

            out_pos =  normal + this->center;
            
            return true;
        }
    }
*/
    return false;
}

