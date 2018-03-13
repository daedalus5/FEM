#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <Partio.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <math.h>
#include <iterator>
#include <fstream>

template<class T, int dim>
class Particles{
public:
	std::vector<Eigen::Matrix<T, dim, 1>> positions;
	std::vector<Eigen::Matrix<T, dim, 1>> velocities;
	std::vector<Eigen::Matrix<T, dim, 1>> forces;
	std::vector<Eigen::Matrix<T, dim, 1>> drags;
	std::vector<T> masses;
    std::vector<int> tets;

	Particles();
	~Particles();

    void zeroForces();
    void addParticle(Eigen::Matrix<T, dim, 1> pos);


};

template<class T, int dim>
Particles<T,dim>::Particles() : positions(), velocities(), forces(), drags(), masses() {}

template<class T, int dim>
Particles<T,dim>::~Particles() {}

template<class T, int dim>
void Particles<T,dim>::zeroForces(){
    for(unsigned int i = 0; i < forces.size(); ++i){
        forces[i] = Eigen::Matrix<T,dim,1>::Zero(dim);
    }
}

template<class T, int dim>
void Particles<T,dim>::addParticle(Eigen::Matrix<T, dim, 1> pos) {
    positions.push_back(Eigen::Matrix<T,dim,1>(pos));
    velocities.push_back(Eigen::Matrix<T,dim,1>(0.0,0.0,0.0));
    masses.push_back(0.0);
    forces.push_back(Eigen::Matrix<T,dim,1>(0.0,0.0,0.0));
    drags.push_back(Eigen::Matrix<T,dim,1>(0.0,0.0,0.0));
}
