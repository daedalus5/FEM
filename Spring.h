#ifndef SPRING_H
#define SPRING_H

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

template<class T, int dim>
class Spring{
public:
	std::pair<int, int> endPoints; // the points to which the spring is attached
	T k;                           // spring constant
	T L0;
	T gamma;

	Spring(int ptA, int ptB, T k, T L0, T gamma);
	~Spring();

	void print_info() const;
};

template<class T, int dim>
Spring<T,dim>::Spring(int ptA, int ptB, T k, T L0, T gamma) : endPoints(std::pair<int, int> (ptA, ptB)), k(k), L0(L0), gamma(gamma) {}

template<class T, int dim>
Spring<T,dim>::~Spring(){}

template<class T, int dim>
void Spring<T,dim>::print_info() const{
	std::cout << "Point A: " << endPoints.first << std::endl;
	std::cout << "Point B: " << endPoints.second << std::endl;
	std::cout << "k: " << k << std::endl;
	std::cout << "L0: " << L0 << std::endl;
	std::cout << "Gamma: " << gamma << std::endl;
}

#endif