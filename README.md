# FEM Fatale
#### A jello simulation using the Finite Element Method
#### Demo video link: 

Contributors:  
Connie Chang - Collisions and rendering  
Marissa Como - I/O, tetgen, test cases, debugging  
Zach Corse - Force calculations, matrix calculations  
Anantha Srinivas - Explicit and implicit integrators  

Implementation details:  
- Fixed Corotated elastic model  
- Even distribution of mass between tetrahedron vertices  
- Collisions using signed distance functions  
- Forward Euler integrator  
- (Work in progress) Backward Euler integrator  
- Rendering in Houdini  

Limitations and Struggles:  
- Physics is not realistic (gravity is very low and Young's modulus is extremely high)  
- Implicit solver does not work  
- Instability in explicit solver as simulation goes on  


Resources:  
SIGGRAPH paper  
Youtube video  
snow paper for implicit  
TA's  
Andre  
Professor Jiang  
