cmake_minimum_required(VERSION 3.5)
# Find eigen in system directory and treat it as an interface library
find_package(Eigen3 REQUIRED)
add_library(eigen INTERFACE IMPORTED)# Generate the local Makefile
if(CREATE_LOCAL_MAKEFILES)
  file(RELATIVE_PATH CISPBA_RELATIVE_PATH ${CISPBA_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  configure_file(${CISPBA_SOURCE_DIR}/Scripts/Makefile.in ${CMAKE_CURRENT_SOURCE_DIR}/Makefile @ONLY)
endif()

add_cispba_executable(example src/main.cpp
        src/FEMSolver.cpp
        src/FEMSolver.h
        src/integrator/ForwardEuler.h
        src/integrator/ForwardEuler.cpp
        src/integrator/BaseIntegrator.h
        src/integrator/BaseIntegrator.cpp
        src/components/Spring.h
        src/components/Spring.cpp
        src/mesh/Mesh.h
        src/mesh/Particles.h
        src/mesh/TetraMesh.h
        src/mesh/Tetrahedron.h
        src/utility/FileHelper.cpp
        src/utility/FileHelper.h
        src/globalincludes.h)
        
target_include_directories(FEM SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR})

target_link_libraries(FEM partio)
