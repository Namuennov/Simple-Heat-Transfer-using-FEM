# Simple-Heat-Transfer-using-FEM
Program simulates heat transfer and solves Fourier Equation

Program can simulate heat transport with data inserted in code - it works only for rectangular grids with evenly distributed nodes.

Program can also do the same task, but with data read data from .txt file - it works even for nonstandard shaped grids. However, it is .txt file author who must guarantee correctness of the values written in input file. Format of the input data must be as given (can be seen in sample file "TestGrid.txt"):

SimulationTime value  
SimulationStepTime value  
Conductivity value  
Alfa value  
Tot value  
InitialTemp value  
Density value  
SpecificHeat value  
NodesNumber value  
ElementsNumber value  
*Node  
      1,  value, value  
      2, value, value  
      3, value, value  
      (...)  
     noNodes-1, value, value  
     noNodes,   value, value  
*Element, type=DC2D4  
 1,  value,  value,  value,  value  
 2,  value,  value,  value,  value  
 3,  value,  value,  value,  value  
 (...)  
 noElements-1, value,  value,  value,  value  
 noElements, value,  value,  value,  value  
*BC  
value, value, (...), value, value  

"value" means numerical value.  
In "*Node" Section there must be 2D coordinates of gird's nodes.  
In "*Element, type=DC2D4" Section there must be numbers of 4 grid nodes for each of gird's elements.  
In "*BC" Section there must be numbers of gird's nodes that are at the edge of the grid (that have boundary condition set on).
