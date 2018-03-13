#pragma once

#include <iostream>
#include <string>

#include "Mesh.h"
#include "Particles.h"
#include "Tetrahedron.h"

template<class T, int dim>
class TetraMesh : public Mesh<T,dim>{
public:
    static double k;
    static double nu;

    TetraMesh(std::string s);
    virtual ~TetraMesh();

    void generateTetras();      // read data from tetgen and populate particles and tetras
    void outputFrame(int frame);    // write data to frame
    void generateSimpleTetrahedron();

    Particles<T,dim> mParticles;
    std::vector<Tetrahedron<T,dim>> mTetras;
};

template<class T, int dim>
TetraMesh<T,dim>::TetraMesh(std::string s) : Mesh<T,dim>(s) {}

template<class T, int dim>
TetraMesh<T,dim>::~TetraMesh(){}

template<class T, int dim>
void TetraMesh<T,dim>::generateTetras(){
    std::ifstream instream; //input file stream

    // output to a poly file for debugging
    std::ofstream outFile;
    outFile.open("out.poly");
    if (!outFile) {
        std::cerr << "Unable to open file out.poly";
        exit(1);
    }
    outFile << "POINTS\n\n";


    std::ofstream outObject;
    outObject.open("out.obj");
    if (!outObject) {
        std::cout << "Unable to open file out.obj";
        exit(1);
    }
    outObject << "default\n\n";

    // .NODE FILE
    // 	list of vertices
            instream.open(this->filepath+".node");
            if (instream.fail())
            {
                    std::cout << "ERROR" << std::endl;
                    exit(1);
            }
            double x1, x2, x3, x4; // variables for parsing

            int numVerts;
            std::string line =  "";
            getline(instream, line);
            const char *l = &line[0];
            numVerts = atoi(l);

            for(int i = 0; i < numVerts; i++)
            {
                    instream >> x1 >> x2 >> x3 >> x4;
                    this->mParticles.positions.push_back(Eigen::Matrix<T,dim,1>(x2,x3,x4));
                    this->mParticles.velocities.push_back(Eigen::Matrix<T,dim,1>(0.0,0.0,0.0));
                    this->mParticles.masses.push_back(0.0);
                    this->mParticles.forces.push_back(Eigen::Matrix<T,dim,1>(0.0,0.0,0.0));
                    this->mParticles.drags.push_back(Eigen::Matrix<T,dim,1>(0.0,0.0,0.0));
                    this->mParticles.tets.push_back(0);

                    outFile << x1 << ": " << x2 << " " << x3 << " " << x4 << "\n";  // .poly
                    outObject << "v " << x2 << " " << x3 << " " << x4 << std::endl; // .obj
            }
            instream.close();

    // .ELE FILE
    //	list of tetrahedra
            instream.open(this->filepath+".ele");
            if (instream.fail())
            {
                   std::cout << "ERROR" << std::endl;
                   exit(1);
            }

            outFile << "\nPOLYS\n\n";

            getline(instream, line);
            const char *t = &line[0];
            int numTets = atoi(t);
            int a,b,c,d,e;
            std::vector<int> indices;


            for(int i = 0; i < numTets; i++)
            {
                    // indices of tetrahedron
                    instream >> a >> b >> c >> d >> e;
                    indices.push_back(b-1);
                    indices.push_back(c-1);
                    indices.push_back(d-1);
                    indices.push_back(e-1);

                    // create tetrahedron instance
                    Tetrahedron<T, dim> tet(indices);
                    this->mTetras.push_back(tet);
                    indices.clear();

                    outFile << (i * 3) + 1 <<": " << b << " " << c << " " << d << " " << b << "\n";
                    outFile << (i * 3) + 2 <<": " << b << " " << e << " " << d << "\n";
                    outFile << (i * 3) + 3 <<": " << c << " " << e << "\n";
            }
            instream.close();


// .FACE file
     // list of faces
            instream.open(this->filepath+".face");
            if (instream.fail())
            {
                   std::cout << "ERROR face" << std::endl;
                   exit(1);
            }

            getline(instream, line);
            const char *f = &line[0];
            numVerts = atoi(f);

            for(int i = 0; i < numVerts; i++)
            {
                    instream >> a >> b >> c >> d >> e;
                    //instream >> a >> b >> c >> d;

                    // output faces
                    outObject << "f " << b << " " << " " << c <<  " " << " " << d << std::endl;
            }
            instream.close();

    outFile << "\nEND";
    outFile.close();
    outObject.close();
}

template<class T, int dim>
void TetraMesh<T,dim>::generateSimpleTetrahedron() {

    this->mParticles.addParticle(Eigen::Matrix<T,dim,1>(1.0, 0.0, 1.0)); // 0
    this->mParticles.addParticle(Eigen::Matrix<T,dim,1>(-1.0, 0.0, 1.0)); // 1
    this->mParticles.addParticle(Eigen::Matrix<T,dim,1>(-1.0, 0.0, -1.0)); // 2
    this->mParticles.addParticle(Eigen::Matrix<T,dim,1>(1.0, 0.0, -1.0)); // 3

    this->mParticles.addParticle(Eigen::Matrix<T,dim,1>(1.0, 1.0, 1.0)); // 4
    this->mParticles.addParticle(Eigen::Matrix<T,dim,1>(-1.0, 1.0, 1.0)); // 5
    this->mParticles.addParticle(Eigen::Matrix<T,dim,1>(-1.0, 1.0, -1.0)); // 6
    this->mParticles.addParticle(Eigen::Matrix<T,dim,1>(1.0, 1.0, -1.0)); // 7

    std::vector<int> indices;

    indices.push_back(4);
    indices.push_back(1);
    indices.push_back(6);
    indices.push_back(3);
    this->mTetras.push_back(Tetrahedron<T, dim>(indices));

    indices.clear();
    indices.push_back(0);
    indices.push_back(1);
    indices.push_back(4);
    indices.push_back(3);
    this->mTetras.push_back(Tetrahedron<T, dim>(indices));

    indices.clear();
    indices.push_back(7);
    indices.push_back(4);
    indices.push_back(6);
    indices.push_back(3);
    this->mTetras.push_back(Tetrahedron<T, dim>(indices));

    indices.clear();
    indices.push_back(2);
    indices.push_back(1);
    indices.push_back(3);
    indices.push_back(6);
    this->mTetras.push_back(Tetrahedron<T, dim>(indices));

    indices.clear();
    indices.push_back(5);
    indices.push_back(6);
    indices.push_back(4);
    indices.push_back(1);
    this->mTetras.push_back(Tetrahedron<T, dim>(indices));
}

template<class T, int dim>
void TetraMesh<T,dim>::outputFrame(int frame){
    Partio::ParticlesDataMutable* parts = Partio::create();
       Partio::ParticleAttribute posH, vH, mH, fH;
       mH = parts->addAttribute("m", Partio::VECTOR, 1);
       posH = parts->addAttribute("position", Partio::VECTOR, 3);
       vH = parts->addAttribute("v", Partio::VECTOR, 3);
       fH = parts->addAttribute("f", Partio::VECTOR, 3);

       // iterate through the particles
          for (unsigned int i=0; i < this->mParticles.positions.size(); i++){
             int idx = parts->addParticle();
             float* m = parts->dataWrite<float>(mH, idx);
             float* p = parts->dataWrite<float>(posH, idx);
             float* v = parts->dataWrite<float>(vH, idx);
             float* f = parts->dataWrite<float>(fH, idx);
             m[0] = this->mParticles.masses[i];
             for (int k = 0; k < 3; k++)
                 p[k] = this->mParticles.positions[i][k];
             for (int k = 0; k < 3; k++)
                 v[k] = this->mParticles.velocities[i][k];
             for (int k = 0; k < 3; k++)
                 f[k] = this->mParticles.forces[i][k];
          }

          // write frames to .bgeo file
          std::string f = std::to_string(frame);
          std::string particleFile = "";
          if(f.length() == 1)
             particleFile = "frame000" + f +".bgeo";
          else if(f.length() == 2)
             particleFile = "frame00" + f +".bgeo";
          else if(f.length() == 3)
             particleFile = "frame0" + f +".bgeo";
          else
             particleFile = "frame" + f +".bgeo";

          particleFile = "output/" + particleFile;
          Partio::write(particleFile.c_str(), *parts);
          parts->release();
}
