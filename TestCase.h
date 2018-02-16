#ifndef TEST_CASE_H
#define TEST_CASE_H

#include <iostream>
#include <vector>
#include <string>
#include <Partio.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <math.h>
#include <iterator>
#include <fstream>
#include <utility>

#include "Spring.h"

template<class T, int dim>
class TestCase{
public:
	Particles<T,dim> particles;
	std::vector<Spring<T,dim>> springs;

	TestCase();
	~TestCase();

	// parses a line containing all relevant particle & spring info
	void initialize(std::string& s);
	void print_info(std::ofstream& file);
private:
	void computeForces();
};

template<class T, int dim>
TestCase<T,dim>::TestCase(){}

template<class T, int dim>
TestCase<T,dim>::~TestCase(){}

template<class T, int dim>
void TestCase<T,dim>::initialize(std::string& s){
	// number of particles
	const int N = 4;

	std::string delimiter = " ";
	std::size_t pos = 0;
	std::string token;
	int n = 0;

	Eigen::Matrix<T,dim,1> x = {0, 0, 0};
	Eigen::Matrix<T,dim,1> v = {0, 0, 0};
	Eigen::Matrix<T,dim,1> f = {0, 0, 0};
	Eigen::Matrix<T,dim,1> fd = {0, 0, 0};
	T k = 0.0;
	T L0 = 0.0;
	T gamma = 0.0;

	for(int i = 0; i < 3; ++i){
		pos = s.find(delimiter);
		token = s.substr(0, pos);
		switch(n){
			case 0 : k = std::stof(token);
				 	 break;
			case 1 : gamma = std::stof(token);
					 break;
			case 2 : L0 = std::stof(token);
					 break;
			default : std::cout << "error" << std::endl;  
					 break;
		}
		s.erase(0, pos + delimiter.length());
		n++;
	}

	for(int i = 0; i < N; ++i){
		for(int j = 0; j < dim; ++j){
			pos = s.find(delimiter);
			token = s.substr(0, pos);
			x[j] = std::stof(token);
			s.erase(0, pos + delimiter.length());
		}
		particles.positions.push_back(x);
	}

	for(int i = 0; i < N; ++i){
		for(int j = 0; j < dim; ++j){
			pos = s.find(delimiter);
			token = s.substr(0, pos);
			v[j] = std::stof(token);
			s.erase(0, pos + delimiter.length());
		}
		particles.velocities.push_back(v);
	}

	for(int i = 0; i < N; ++i){
		for(int j = 0; j < dim; ++j){
			particles.forces.push_back(f);
			particles.drags.push_back(fd);
			particles.masses.push_back((T)1.0);
		}
	}

	for(int i = 0; i < N - 1; ++i){
		Spring<T,dim> spring(i, i + 1, k, L0, gamma);
		springs.push_back(spring);
	}
	Spring<T,dim> springA(3, 0, k, L0, gamma);
	springs.push_back(springA);
	Spring<T,dim> springB(1, 3, k, L0, gamma);
	springs.push_back(springB);

	computeForces();
}

template<class T, int dim>
void TestCase<T,dim>::print_info(std::ofstream& file){
	for(unsigned int i = 0; i < particles.positions.size(); ++i){
		for(int j = 0; j < dim; ++j){
			file << particles.forces[i][j] << " ";
		}
	}
	for(unsigned int i = 0; i < particles.positions.size(); ++i){
		for(int j = 0; j < dim; ++j){
			file << particles.drags[i][j] << " ";
		}
	}
	file << std::endl;
}

template<class T, int dim>
void TestCase<T,dim>::computeForces(){
	Eigen::Matrix<T,dim,1> x1 = {0, 0, 0};
	Eigen::Matrix<T,dim,1> x2 = {0, 0, 0};
	Eigen::Matrix<T,dim,1> diffX = {0, 0, 0};
	Eigen::Matrix<T,dim,1> n = {0, 0, 0};
	Eigen::Matrix<T,dim,1> v1 = {0, 0, 0};
	Eigen::Matrix<T,dim,1> v2 = {0, 0, 0};
	Eigen::Matrix<T,dim,1> diffV = {0, 0, 0};
	Eigen::Matrix<T,dim,1> f1 = {0, 0, 0};
	Eigen::Matrix<T,dim,1> f2 = {0, 0, 0};
	Eigen::Matrix<T,dim,1> fd1 = {0, 0, 0};
	Eigen::Matrix<T,dim,1> fd2 = {0, 0, 0};
	T k = 0.0;
	T L0 = 0.0;
	T gamma = 0.0;
	T length = 0.0;

	int particle1;
	int particle2;

	for(unsigned int i = 0; i < springs.size(); ++i){
		k = springs[i].k;
		L0 = springs[i].L0;
		gamma = springs[i].gamma;

		particle1 = springs[i].endPoints.first;
		particle2 = springs[i].endPoints.second;

		x1 = particles.positions[particle1];
		x2 = particles.positions[particle2];
		v1 = particles.velocities[particle1];
		v2 = particles.velocities[particle2];

		diffX = x1 - x2;
		length = diffX.norm();
		n = diffX / length;

		f1 = -k * ((length / L0) - 1) * n;
		f2 = -f1;
		particles.forces[particle1] += f1;
		particles.forces[particle2] += f2;

		diffV = v1 - v2;

		fd1 = -gamma * (diffV.dot(n)) * n;	
		fd2 = -fd1;
		particles.drags[particle1] += fd1;
		particles.drags[particle2] += fd2;
	}
}

#endif