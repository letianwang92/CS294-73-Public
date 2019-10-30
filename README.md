Simulate Crystal Growth with Phase Field Method
CS294-74 Final Project  
Fall 2015 
University of California, Berkeley  
[Latex Report](https://www.overleaf.com/3886930xwwgst)


# Abstract: 
Crystal growth is a common phenomena in chemical engineering and manufacturing (semicondutor and 3D printing). Simulation of crystal growth will be of great value to the scientific research and industrial application. However, it is not easy to precisely describe the fine features of crystal growth, due to the complicated coupling of heat transfer, mass transport and thermodynamics. Dendritic growth is one of the challenging topics which is considered as the benchmark of crystal growth simulation. The challenge of simulation lies on the precise definition of the interface.  Here we provide a simulation toolbox using phase filed method to define the interfaces. Some test investigations are given.

# Content

A list of included components can be found at [here](https://github.com/letianwang92/CS294-73-Public/tree/master/code/include)

[Dendritic Growth] (https://github.com/letianwang92/CS294-73-Public/tree/master/code/src/Dendritic)
* Classes required for generating temperature and phase field grid, propagate the required values.

[RectMDArray](https://github.com/letianwang92/CS294-73-Public/tree/master/code/src/RectArray)
* The framework of numerical grid is borrowed from the course materials

[Tests](https://github.com/letianwang92/CS294-73-Public/tree/master/code/src/unitTests)
* Unit test and performance test of the numerical algorithm


