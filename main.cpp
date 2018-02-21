#include "FEMSolver.h"
#include "globalincludes.h"

//#define TEST_EIGEN
//#define TEST_PARTIO

void checkIfEigenWorks() {
    MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    std::cout << m << std::endl;

    std::cout << "<<<< EIGEN LOADED SUCCESSFULLY >>>> " << std::endl;


}

template <class T, int dim>
void checkIfPartioWorks() {

    std::string particleFile = "test.bgeo";

    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, vH, mH;
    mH = parts->addAttribute("m", Partio::VECTOR, 1);
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    vH = parts->addAttribute("v", Partio::VECTOR, 3);
    for (int i=0; i<3; i++){
        int idx = parts->addParticle();
        float* m = parts->dataWrite<float>(mH, idx);
        float* p = parts->dataWrite<float>(posH, idx);
        float* v = parts->dataWrite<float>(vH, idx);
        m[0] = (T)(i + 1);
        for (int k = 0; k < 3; k++)
            p[k] = (T)(i + 1);
        for (int k = 0; k < 3; k++)
            v[k] = (T)(i + 100);
    }

    Partio::write(particleFile.c_str(), *parts);
    parts->release();

    std::cout << "<<<< PARTIO LOADED SUCCESSFULLY >>>> " << std::endl;

}


int main()
{

#ifdef TEST_EIGEN
    checkIfEigenWorks();
#endif

#ifdef TEST_PARTIO
    checkIfPartioWorks<float, 3>();
#endif

    // Cook My Jello!

    FEMSolver<float,3> solver(120);
    solver.initializeMesh();
    solver.cookMyJello();

}