#pragma once

#include <iostream>
#include <vector>
#include <string>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <math.h>

template<class T, int dim>
class Tetrahedron{
public:
	std::vector<int> mPIndices;         // tetrahedron vertices
    Eigen::Matrix<T,dim,dim> mDm;
    Eigen::Matrix<T,dim,dim> mDmInv;    // rest configuration Dm inverse
    T volume;                           // volume of tetrahedron in rest configuration
    Eigen::Matrix<T,dim,dim> mVolDmInvT;// volume * Dm inverse
    T mass;                             // mass of the tetrahedron

	Tetrahedron(const std::vector<int>& indices);
	~Tetrahedron();

    // precompute populates Dm inverse matrix and tetrahedron volume in rest configuration
    void precompute();
	void print_info() const;           // for debugging
};

template<class T, int dim>
Tetrahedron<T,dim>::Tetrahedron(const std::vector<int>& indices) :
                                mPIndices(indices),
                                mDm(Eigen::Matrix<T,dim,dim>::Zero(dim,dim)),
                                mDmInv(Eigen::Matrix<T,dim,dim>::Zero(dim, dim)),
                                volume(0.0f),
                                mVolDmInvT(Eigen::Matrix<T,dim,dim>::Zero(dim, dim)),
                                mass(0.0f) {}

template<class T, int dim>
Tetrahedron<T,dim>::~Tetrahedron(){}

template<class T, int dim>
void Tetrahedron<T,dim>::precompute(){
    switch(dim){
        case 2: volume = std::abs(mDm.determinant() / 2.f); break;
        case 3: volume = std::abs(mDm.determinant()) * 0.16666666666667f; break;
        default: std::cout << "error: dimension must be 2 or 3" << std::endl;
    }
    if(volume != 0){    // check if Dm is nonsingular
        mDmInv = mDm.inverse();
    }
    else{
        std::cout << "bad tetrahedron" << std::endl;
    }
    mVolDmInvT = volume * (mDmInv.transpose());
    mass = 1000.f * volume; // density is 1000
}

template<class T, int dim>
void Tetrahedron<T,dim>::print_info() const{
    std::cout << mDmInv << std::endl;
    std::cout << volume << std::endl;
}
