#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <Partio.h>
#include <Core>
#include <Dense>
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

	Particles();
	~Particles();

	void write_data_bgeo(const std::string& s);
};

template<class T, int dim>
Particles<T,dim>::Particles() : positions(), velocities(), forces(), drags(), masses() {}

template<class T, int dim>
Particles<T,dim>::~Particles() {}
