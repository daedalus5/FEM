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
	std::vector<int> mPIndices;        // tetrahedron vertices
    Eigen::Matrix<T,dim,dim> mDmInv;   // rest configuration Dm inverse
    T volume;                          // volume of tetrahedron in rest configuration
    Eigen::Matrix<T,dim,dim> mVolDmInvT;// volume * Dm inverse

	Tetrahedron(const std::vector<int>& indices);
	~Tetrahedron();

    // precompute populates Dm inverse matrix and tetrahedron volume in rest configuration
    void precompute(const std::vector<Eigen::Matrix<T,dim,1>>& x);
	void print_info() const;           // for debugging
};

template<class T, int dim>
Tetrahedron<T,dim>::Tetrahedron(const std::vector<int>& indices) :
                                mPIndices(indices),
                                mDmInv(Eigen::Matrix<T,dim,dim>::Zero(dim, dim)),
                                volume(0.0),
                                mVolDmInvT(Eigen::Matrix<T,dim,dim>::Zero(dim, dim)) {}

template<class T, int dim>
Tetrahedron<T,dim>::~Tetrahedron(){}

template<class T, int dim>
void Tetrahedron<T,dim>::precompute(const std::vector<Eigen::Matrix<T,dim,1>>& x){
    Eigen::Matrix<T,dim,dim> Dm = Eigen::Matrix<T,dim,dim>::Zero(dim,dim);
    for(int i = 1; i < dim + 1; ++i){
        Dm.col(i - 1) = x[i] - x[0];
    }
    switch(dim){
        case 2: volume = std::abs(Dm.determinant() / 2.); break;
        case 3: volume = std::abs(Dm.determinant() / 6.); break;
        default: std::cout << "error: dimension must be 2 or 3" << std::endl;
    }
    if(volume != 0){    // check if Dm is nonsingular
        mDmInv = Dm.inverse();
    }
    else{
        std::cout << "bad tetrahedron" << std::endl;
    }
    mVolDmInvT = -1 * volume * mDmInv.transpose();
}

template<class T, int dim>
void Tetrahedron<T,dim>::print_info() const{
    std::cout << mDmInv << std::endl;
    std::cout << volume << std::endl;
}
