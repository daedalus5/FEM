#pragma once

#include <iostream>
#include <string>

template<class T, int dim>
class Mesh{
public:
    std::string filepath;   // the mesh filepath
    Mesh(std::string s);
    virtual ~Mesh();
};

template<class T, int dim>
Mesh<T,dim>::Mesh(std::string s) : filepath(s) {}

template<class T, int dim>
Mesh<T,dim>::~Mesh(){}