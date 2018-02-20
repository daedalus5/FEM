#pragma once

#include <iostream>
#include <string>

#include "Mesh.h"
#include "Particles.h"
#include "Tetrahedron.h"

template<class T, int dim>
class TetraMesh : public Mesh<T,dim>{
public:
    TetraMesh(std::string s);
    virtual ~TetraMesh();

    void generateTetras();      // read data from tetgen and populate particles and tetras
    void outputFrame(int i);    // write data to frame i
private:
    Particles<T,dim> mParticles;
    std::vector<Tetrahedron> mTetras;
};

template<class T, int dim>
TetraMesh<T,dim>::TetraMesh(std::string s) : Mesh<T,dim>(s) {}

template<class T, int dim>
TetraMesh<T,dim>::~TetraMesh(){}

template<class T, int dim>
void TetraMesh<T,dim>::generateTetras(){}

template<class T, int dim>
void TetraMesh<T,dim>::outputFrame(int i){}