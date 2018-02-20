#pragma once

#include <iostream>
#include <vector>
#include <string>

#include <Core>
#include <Dense>
#include <math.h>

class Tetrahedron{
public:
	std::vector<int> mPIndices;     // tetrahedron vertices

	Tetrahedron(const std::vector<int>& indices);
	~Tetrahedron();

	void print_info() const;
};

Tetrahedron::Tetrahedron(const std::vector<int>& indices) : mPIndices(indices) {}

Tetrahedron::~Tetrahedron(){}

void Tetrahedron::print_info() const{
}
