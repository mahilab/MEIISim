# MEIISim
Code for generating kinematics and dynamics of the **[Mahi Exo II][1]**, and simulating them.

## Overview
---
This repository has two main sections to work in. 
The first is the **[matlab folder][2]**. 
This folder contains information regarding the actual mathematics for the kinematics and dynamics of the MEII. 
The main files in this folder are the files [`MEII_Kinematics.m`][3], and [`MEII_Dynamics.m`][4]. 
These contain the information for generating code with respect to the kinematics and dynamics.

## Kinematics 
---
The kinematics is based on the papers *3RPS Kinematics Lee Shah* and *3RPS Dynamics Lee Shah* in the [Matlab/Papers folder][5]. 
These papers detail the methods for finding the important platform position, velocity, orientation, and angular velocity in closed form for FK, and mostly closed for for IK.
The math behind doing this for the non-symmetric Mahi Exo II is in the file [`MEII_Kinematics.m`][3], and the resulting functions, [`MEII_FK.m`][6] and [`MEII_IK.m`][7]

## Dynamics
---
The dynamics are loosely based on the *3RPS Dynamics Lee Shah* paper in the [Matlab/Papers folder][5]. 
The paper details most of the requirements to arrive at the full dynamics solution for a *symmetric* 3RPS mechanism. 
In this paper, the three prismatic joints of the 3RPS mechanism are considered the independent coordinates, and the three rotational coordinates (not the spherical joints) are considered three dependent coordinates.
With equations of motion (EOM) defined in terms of the those 6 coordinates, we can convert the EOMs to be in terms of only the independent coordinates using the result from the *Parallel Kinematics Ghorbel* paper found in the [Matlab/Papers folder][5]. 
This results in the full explicit EOMs (5x5 M, 5x5 V, 5x1 G) that describe the full system, but these expressions are extremely long and tedious to work with (need to do multithreading for them to be solved fast enough).

The EOMs can instead be solved directly using the result of the *Parallel Kinematics Ghorbel* paper with using the three independent prismatic coordinates in addition to 9 dependent coordinates (as opposed to 3 in the previous method). 
This allows us to express the EOMs in terms of matrices of size (14x14 M, 14x14 V, 14x1 G) and then convert these to explicit size (5x5 M, 5x5 V, 5x1 G) numerically using matrices defined in *Parallel Kinematics Ghorbel*.
Because this conversion happens numerically rather than symbolically, the computation time is extremely reduced (~10x).
In either version, the EOMs are output to the [Matlab/DynamicEqs][8] folder as .txt files.
These .txt files are then written to a c++ readable format using the file `correctequations.py` in the [matlab folder][2], which adds the appropriate includes and format for accessing the results of these matrices.
These files, {M, G, psi, psi_dq, psi_dq_dt}.hpp, are included in the cpp project which mainly consists of the files `MeiiModel.cpp/hpp` and either `dll.cpp` or `dll_virtual.cpp` in the [src][9] and [include][10] folders.
The files `dll.cpp` and `dll_virtual.cpp` compile into the dlls `meii_model.dll` and `virtual_meii.dll` respectively. 

## Unity Project
---
The unity project consists of a couple key components. 
The first main component is the 3D model of the Mahi Exo II which is generated from solidworks using the [solidworks obj exporter v2.0][11]. 
The obj files are arranged in the unity project so that there is an empty gameobject which controls the orientation of the next joint based on the orientation of the preceding joint. 
When the unity program is started, the dll is started using from the [MeiiScript.cs][12] which starts the simulation.
The orientation is then controlled through the script [MeiiScript.cs][12] which pulls the orientations of each of the joints from the appropriate dll, and writes those values to the correct controlling gameobject.
Two different versions have been built to use as a standalone .exe in its respective folder.

### **MEII Virtual**
The **Virtual** folder in the [Build subfolder][13] of the unity project can be used in combination with a virtual MEII object in the [MEII codebase][1].
This uses the `virtual_meii.dll` codebase to simulate and communicate between the sim and MEII script.
This will run a standard MEII script and show the output on the simulation of the MEII. 

### **MEII Tuner**
The **Tuner** folder in the [Build subfolder][13] of the unity project can be used in combination with a virtual MEII object in the [MEII codebase][1].
This uses the `meii_model.dll` codebase to simulate and communicate with the `tuner.exe` file located in the same folder.
This tuner can be used to play with control parameters and to adjust reference positions to see how the exo responds.

## Todo
---
- add ability to output different versions of dynamics
  - fully decomposed (currently the setup and very efficient)
  - fully composed (raw equations for 5x5 M, V, and 5x1 G and very inefficient)
- add in comments to the kinematics and dynamics file
- add in comments to the cpp files
- test ability to use params in cpp files instead of substituting in values early
  - this will enable ability to possibly learn params online?


[1]: https://github.com/mahilab/MEII
[2]: https://github.com/mahilab/MEIISim/tree/master/matlab
[3]: https://github.com/mahilab/MEIISim/blob/master/matlab/MEII_Kinematics.m
[4]: https://github.com/mahilab/MEIISim/blob/master/matlab/MEII_Dynamics.m
[5]: https://github.com/mahilab/MEIISim/tree/master/matlab/Papers
[6]: https://github.com/mahilab/MEIISim/tree/master/matlab/MEII_FK.m
[7]: https://github.com/mahilab/MEIISim/tree/master/matlab/MEII_IK.m
[8]: https://github.com/mahilab/MEIISim/tree/master/matlab/DynamicEqs
[9]: https://github.com/mahilab/MEIISim/tree/master/src
[10]: https://github.com/mahilab/MEIISim/tree/master/include
[11]: https://forum.solidworks.com/thread/54270
[12]: https://github.com/mahilab/MEIISim/blob/master/unity/MEIISim/Assets/MeiiScript.cs
[13]: https://github.com/mahilab/MEIISim/blob/master/unity/MEIISim/Build