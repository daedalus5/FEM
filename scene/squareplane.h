// base copied from CIS 561

#pragma  once
#include "shape.h"

//A SquarePlane is assumed to have a center of <0,0,0> and is horizontal.
template<class T, int dim>
class SquarePlane : public Shape<T, dim>
{
    public:
        SquarePlane() : Shape<T, dim>() {};
        SquarePlane(std::string file);
        bool checkCollisions(const Eigen::Matrix<T, dim, 1> &pos, Eigen::Matrix<T, dim, 1> &out_pos) const override;

    private:
        float length_half = 200.f;
};


template<class T, int dim>
SquarePlane<T, dim>::SquarePlane(std::string file) : Shape<T, dim>(file) {
    this->center[1] = -1.0;
}

template<class T, int dim>
bool SquarePlane<T, dim>::checkCollisions(const Eigen::Matrix<T, dim, 1> &pos, Eigen::Matrix<T, dim, 1> &out_pos) const
{
    //for (unsigned int i = 0; i < particles.positions.size(); ++i) {

        // If particle is under plane
        if (pos[1] < this->center[1] - 0.00001) {

            // Check if x, z is within plane
            if (std::abs(pos[0] - this->center[0]) < length_half
                || (dim == 3 && std::abs(pos[2] - this->center[2]) < length_half)) {
                //TODO is there a better way to see if this is 2D or 3D?

                // Collision
                // Move particle out of plane
                //particles.positions[i][1] = center[1];
                // Update velocity?

                // fill out_pos with new position
                // Some function outside this will calculate velocity
                out_pos[0] = pos[0];
                out_pos[1] = this->center[1];
                if (dim == 3) {
                    out_pos[2] = pos[2];
                }
                return true;

            }

        }
    //}
    return false;
}

