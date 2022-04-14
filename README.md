# Symplectic_Integrators_EvalTask
This repository contains evaluation task solution regarding Geant4 - Symplectic Integrators for Google Summer of Code 2022

-------
All the files are written taking **SI unit** into account. Also important lines of code in the files have respectively comments on the reason and understanding of the code.

## Contents

  1. **Task1-3.m** - This file implements the electron's trajectory in an Electric and Magnetic Field, using *Runge-Kutta 4th Order Integrator* applied on the Ordinary Differential Equation in terms of electron's momentum (velocity). A 3D plot of electron's path with respect to time is also animated for better visualization. The file also calculates the time period `Tp` and radius of curvature of electron's trajectory `R` as required. Apart from the calculation, it also analyzes the output of the used integrator compared to analytical solution and gives a variable `err_percent`, which gives percentage error for each component of position and momentum for every timestep.
  2. **Task4.m** - This file calculates the position and momentum of the electron after **10**, **100**, **1000** and **10000** revolution and stores the data in the variables `Pos` and `Momen`.
  
      **Note** - For this file to work without errors, first file `Task1-3.m` should be executed and then this file should run.
  
  3. **High Order Integrators** - This folder contains 2 files, `RK_6_Order.m` and `DoPri5.m`, which respectively exectues **Runge-Kutta 6th Order** and **Dormant-Prince745 integrators** on the electron's ODE, these files contain a visualization as a 3D plot for the electron's path.
  4. **leapfrog.m** - This file implements **Leapfrog integrator** method from scratch taking into account of all the requirements specified in the task. Comments are made for better understanding on the working of the algorithm.
  5. **C++ files** - The folder contains header and source files for implementation of **Runge Kutta 4th Order**, **Leapfrog** and **DoPri5** method and outputs component-wise position and momentum of the electron at each timestep. This folder specifically contains 2 header files, `integrators.h` and `cross_p.h` which contains methods for above mentioned 3 integrators and methods for 1D and 2D vector operations respectively. Also it contains a C++ source file, `integrators.cpp` which implements the mentioned integration methods and prints the corresponding outputs.


## Important Results

* Expected Result of electron's trajectory using various numerical integrators:

    ![image](https://user-images.githubusercontent.com/69386934/161440164-2b53a04c-457d-45ba-8bfa-12927515ad38.png)

* Comparing Runge-Kutta 4th Order solution and Analytical solution. Here first three columns represent the relative error percentage in x, y and z positions with respect to analytical solution and last three columns represent the same for x, y and z momentum, whereas, the rows represent the timestep. It is observed that error percentage for inital timesteps for all the components are very negligible as seen in the figure. This implies that 4th Order RK method in this perticular case of electron's ODE works well for the choosen timestep of around 10^-12. (NAN values just represent 0/0 which is meaningless)

    ![image](https://user-images.githubusercontent.com/69386934/161437498-80faddd9-191e-4f36-b21c-6c5fd9c3ea61.png)
    
* Analyzing position and momentum after 10, 100, 1000 and 10000 turns for timestep of 0.25 times time period, it is observed that the path calculated for the electron is very different as compared theoretically as can be seen in the figure. The points shown here represents the 3D position of the electron calculated at every timestep where clearly we can see that taking timestep to be one-forth of time period giving 4 points for every revolution. This wierd behavior of the plot is due to taking highly inaccurate timestep which plots only 4 points for one full turn which 4th Order RK method can't use to calculate eficiently. Moreover it is also observed that due to the scaling of this timestep, momentum and positoin values after 10, 100, 1000 and 10000 turns converges towards 0 which is not actually what should happen. The error in the values is due to the timestep taken. 

    ![image](https://user-images.githubusercontent.com/69386934/161438498-e3a84840-f1b9-48ea-91ca-51be8d279444.png)
    
* When comparing High Order integrators from the 4th Order RK method solution, it is seen that error percentage of the solution from RK 6th Order is very small and 6th Order integrator is more accurate compared to 4th Order, taking all the parameters and timestep to be same for both. Below figure shows the component wise position and momentum error percentage for each timestep between 4th and 6th Order RK method. (The values represent relative error percetage of 6th Order RK method with reprect to  4th Order RK method where first 3 columns represent x, y and z positions and last 3 columns as x, y and z momentum, rows represent the timestep)

    ![image](https://user-images.githubusercontent.com/69386934/161439096-a9ccff94-4d4e-4fd1-b2f9-4b9bcd0366b8.png)  
    
* It can also be seen by comparing with analytical solution, that RK 6th Order is very accurate compared to RK 4th Order as shown from below figure. Below figure shows the percentage error of RK 6th Order solution with analytical solution. (The values represent relative error percentage of the 6th Order RK method with respect to the analytical solution)

    ![image](https://user-images.githubusercontent.com/69386934/161439426-4075ea42-6829-46b0-a7b2-688f670e8dae.png)
    
* Also comparing DoPri5 integrator with RK 4th Order, we can see that again there is very less error percentage between the two, but DoPri5 method is more accurate while calculating the kinematic parameters. This is because, DoPri5 method updates its timestep based on the absolute error if the position and momentum were calculated by 4th Order RK and 6th Order RK. One drawback is that it is very slow as compared to other numerical integration methods.

## Second Order Symplectic Integration for Hamiltonian Problems

### Stormer-Verlet Method

This method is used widely for calculating particle or any object trajectories by observing the governing hamiltonian dynamics of time evolution, where verlet method works quite impressive as for normal eular method integration order of global error is 1 whereas for midpoint method used by velocity verlet, order of error is 2. Also due to Stormer-Verlet method being time reversible in nature, that is we can get certain parameters like position of previous timestep using a futuretimestep parameter values, such property is also observed in Leapfrog method.

Reference - https://en.wikipedia.org/wiki/Verlet_integration

## Second Order Symplectic Integration for Motion in an Electomagnetic field

### Heun's Method (Runge-Kutta method)

This method is based on modified Eular's method which repeatively uses line tangent approximate for the function for a given time sample and then it uses that to estimate the slope and calculate the future position. Since motion in an electromagnetic field has a linear and circular motion combined, the path throughout is curved in nature in 3 dimension, and thus high precision estimation of the slope for a given timestep is very important to find future position of theh particle and hence heun's method can be used to integrate the particle's corresponding differential equation in an electromagnetic field.

Reference - https://www.hindawi.com/journals/ahep/2018/8621573/
