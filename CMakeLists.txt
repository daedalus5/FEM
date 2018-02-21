cmake_minimum_required(VERSION 3.5)
# Find eigen in system directory and treat it as an interface library
find_package(Eigen3 REQUIRED)
add_library(eigen INTERFACE IMPORTED)# Generate the local Makefile
if(CREATE_LOCAL_MAKEFILES)
  file(RELATIVE_PATH CISPBA_RELATIVE_PATH ${CISPBA_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  configure_file(${CISPBA_SOURCE_DIR}/Scripts/Makefile.in ${CMAKE_CURRENT_SOURCE_DIR}/Makefile @ONLY)
endif()

add_cispba_executable(FEM main.cpp
        FEMSolver.cpp
        FEMSolver.h
        integrator/ForwardEuler.h
        integrator/ForwardEuler.cpp
        integrator/BaseIntegrator.h
        integrator/BaseIntegrator.cpp
        components/Spring.h
        components/Spring.cpp
        mesh/Mesh.h
        mesh/Particles.h
        mesh/TetraMesh.h
        mesh/Tetrahedron.h
        utility/FileHelper.cpp
        utility/FileHelper.h
        scene/shape.h
        scene/squareplane.h
        globalincludes.h)
        
target_include_directories(FEM SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR})

target_link_libraries(FEM partio)
