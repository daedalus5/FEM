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
	void write_data_bgeo(const std::string& s, int frame);


    void addParticle(Eigen::Matrix<T, dim, 1> pos);


};

template<class T, int dim>
Particles<T,dim>::Particles() : positions(), velocities(), forces(), drags(), masses() {}

template<class T, int dim>
Particles<T,dim>::~Particles() {}

template<class T, int dim>
void Particles<T,dim>::zeroForces(){
    for(unsigned int i = 0; i < forces.size(); ++i){
        for(int j = 0; j < dim; ++j){
            forces[i][j] = 0;
        }
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

template<class T, int dim>
void Particles<T,dim>::write_data_bgeo(const std::string& s, int frame){
	std::string file = s + std::to_string(frame) + ".bgeo";
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, vH, mH;
    mH = parts->addAttribute("m", Partio::VECTOR, 1);
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    vH = parts->addAttribute("v", Partio::VECTOR, 3);

    for(unsigned int i = 0; i < positions.size(); ++i){
        int idx = parts->addParticle();
        float* m = parts->dataWrite<float>(mH, idx);
        float* p = parts->dataWrite<float>(posH, idx);
        float* v = parts->dataWrite<float>(vH, idx);
        m[0] = (T)masses[i];
        for (int k = 0; k < 3; ++k){
        	p[k] = (T)positions[i][k];
        	v[k] = (T)velocities[i][k];
        }
    }
    Partio::write(file.c_str(), *parts);
    parts->release();
}
