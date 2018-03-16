# FEM Fatale
#### A jello simulation using the Finite Element Method
#### Demo video link: https://youtu.be/rfcBbrljoug

Contributors:  
Connie Chang - Collisions and rendering  
Marissa Como - I/O, tetgen, test cases, debugging  
Zach Corse - Force calculations, K matrix calculations  
Anantha Srinivas - Explicit and implicit integrators  

Implementation details:  
- Fixed Corotated elastic model  
- Even distribution of mass between tetrahedron vertices  
- Collisions using signed distance functions  
- Forward Euler integrator  
- (Work in progress) Backward Euler integrator 
- Tetgen file input reader
- OBJ output for rendering
- Rendering in Houdini  

Limitations and Struggles:  
- Implicit solver does not work...yet
- Instability in explicit solver as simulation goes on (>3 bounces) 

Resources:  
2012 SIGGRAPH course notes on FEM
Professor Kavan's Youtube videos
MPM Snow Simulation paper for implicit  
TA's Josh and Ziyin  
Andre  
Professor Jiang  
