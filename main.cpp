#include "FEMSolver.h"
#include "globalincludes.h"


int main()
{
    // Cook My Jello!

    FEMSolver<double,3> solver(10);
    solver.initializeMesh();
    solver.cookMyJello();

}
