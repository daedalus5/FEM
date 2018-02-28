#pragma once

#include "shape.h"
#include "squareplane.h"
#include "sphere.h"

template<class T, int dim>
class Scene
{
    public:
        Scene();
        virtual ~Scene();
        bool checkCollisions(const Eigen::Matrix<T, dim, 1> pos, Eigen::Matrix<T, dim,1> &out_pos) const;
        void outputFrame(int currFrame);
        void updatePosition(T dt);
        
        std::vector<Shape<T, dim>*> shapes;
};

template<class T, int dim>
Scene<T, dim>::Scene() {

    // Create a ground plane
    // Shape<T, dim>* ground = new SquarePlane<T, dim>();
    // shapes.push_back(ground);

    // Create sphere
    // Shape<T, dim>* sphere = new Sphere<T, dim>();
    // shapes.push_back(sphere);

}

template<class T, int dim>
Scene<T, dim>::~Scene() {
    for (unsigned int i = 0; i < shapes.size(); ++i) {
        delete shapes[i];
    }
}

template<class T, int dim>
bool Scene<T, dim>::checkCollisions(const Eigen::Matrix<T, dim, 1> pos, Eigen::Matrix<T, dim, 1> &out_pos) const {

    bool collide = false;
    Eigen::Matrix<T, dim, 1> temp_pos;
    for (unsigned int i = 0; i < dim; ++i) {
        temp_pos[i] = T(pos[i]);
    }
    for (unsigned int i = 0; i < shapes.size(); ++i) {
        if (shapes[i]->checkCollisions(pos, out_pos)) {
            
            collide = true;
            for (unsigned int j = 0; j < dim; ++j) {
                temp_pos[j] = T(out_pos[j]);
            }
        }
    }
    return collide;
}

template<class T, int dim>
void Scene<T, dim>::outputFrame(int currFrame) {
    for (unsigned int i = 0; i < shapes.size(); ++i) {
        shapes[i]->outputFrame(currFrame);
    }
}

template<class T, int dim>
void Scene<T, dim>::updatePosition(T dt) {
    for (unsigned int i = 0; i < shapes.size(); ++i) {
        shapes[i]->updatePosition(dt);
    }
}
