This is a framework solving excitonic drift-diffusion equations.

The aim of this code is to solve electron and hole transportation problem in organic solar cell.
With certain modification, it is able to solve electric transportation in other devices.
The weak form of the governing equations are explained in detail in the file [XDD_Eqs.pdf](XDD_Eqs.pdf) in this folder.

# Dependencies

## A C++11-compliant compiler.

The code has been tested with: gcc (4.8.5+), clang (3.8+), and the Intel compiler (ver 16+). Any compiler that supports the C++11 standard should work.

## CMake 3.5 or above.

```bash
sudo apt-get install cmake
```

## PETSc 3.6.4

Requires ParMETIS, mumps, etc.
Before you do anything, please first set the variables:  (NOTE: when you see `/path/to_some_dir`, you have to replace this with real path on your computer)

```bash
export PETSC_DIR=/path/to/petsc-3.6.4    
export PETSC_ARCH=arch-linux2-opt
```

Recommended PETSc configuration options:

```bash
./configure --with-clanguage=cxx --download-scalapack=1 --download-mpich= --download-parmetis=1 --with-scalar-type=real --download-fblaslapack=1 --download-blacs= --download-metis=1 --download-mumps=1 --download-mpich=1
```

Configure and then follow the instructions at the terminal.

```bash
# last step is to put all enviorenment vaiable to bashrc
echo "export PETSC_DIR=/path/to/petsc-3.6.4" >> ~/.bashrc
echo "export PETSC_ARCH=arch-linux-opt" >> ~/.bashrc
source ~/.bashrc
```

## Install libconfig 1.5

```bash
wget http://www.hyperrealm.com/libconfig/libconfig-1.5.tar.gz
tar xvf libconfig-1.5.tar.gz
rm libconfig-1.5.tar.gz

cd libconfig-1.5   
export LIBCONFIG_DIR=`pwd`/install
./configure --prefix=$LIBCONFIG_DIR
make
make install
# Save LIBCONFIG_DIR and add the shared object file to LD_LIBRARY_PATH
echo "LIBCONFIG_DIR=$LIBCONFIG_DIR" >> ~/.bashrc
echo "LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$LIBCONFIG_DIR/lib" >> ~/.bashrc
source ~/.bashrc
```

## Install HDF5 (optional, but if you want to save storage space on disk then you need this)

```bash
wget https://support.hdfgroup.org/ftp/HDF5/current18/src/hdf5-1.8.19.tar
tar xvf hdf5-1.8.19.tar
rm hdf5-1.8.19.tar 
cd hdf5-1.8.19
CC=mpicc ./configure --enable-parallel --prefix=`pwd`/install
```

## Install TalyFEM 7.2.0

```bash
# If you have access to our git repository, do:
git clone https://bitbucket.org/baskargroup/taly_fem
cd taly_fem
git checkout 7.2.0  # get the correct version of talyFEM

# if you receive a zip file, unzip it.
tar xvf taly_fem.tar.gz
cd taly_fem

# Build
mkdir build
cd build

# if without HDF5
cmake ..
make -j8
# if HDF5 installed
export HDF5_ROOT=/path/to/hdf5/install
cmake -DENABLE_HDF5=YES ..

# Test
ctest -j4
```

# Compiling this excitonic-drift-diffusion code

Run the following commands in order.

```bash
# if you have access to our git repository, do:
git clone https://bitbucket.org/baskargroup/drift_diffusion_uncoupled
cd drift_diffusion_uncoupled

# if you instead received a zip file, unzip it.
tar xvf drift_diffusion_uncoupled.tar.gz
cd drift_diffusion_uncoupled

mkdir build
cd build
cmake .. -Dtalyfem_DIR=/path/to/talyfem/build
make -j4
```

If everything is correct, you will see an executable file "DD" has been generated (under build directory). This is our program.

# Using the Program

## Input files

### Morphology file

The morphology file stores the physical morphology as numbers: 0 for donor and 1 for acceptor.
The essential requirements of the morphology file are:

a) Only values between 0-1 are allowed. The interface of the donor and acceptor material is denoted by the number 0.5.

b) The top surface of the domain is connected to the hole collecting electrode and the bottom surface is connected to electron collecting electrode.

c) The pixel (or each physical point) sizes along all directions should be the same, eg. in 2D, each pixel represents a physical square area.

Current code accepts two types of morphology files:

a) '*.bits' file which stores the morphology in one single line.
The morphology should contain indicators either 0 (for donor) or 1 (for acceptor).
The top surface of the domain is supposed to be anode (hole collector) and the bottom
surface is the cathode. The indicators (0 or 1) are stored in the row-wise format (left to right, top to bottom).

For example, in a 4 pixel by 2 pixel morphology file you have:

    1 0 1 1 0 0 1 0

The corresponding morphology (or the physical system) in 2D is:

    ------------anode
     1 0
     0 0
     1 1
     1 0
    -----------cathode

b) The `*.plt` file stores the morphology on a finite element mesh, which means it has both 
information about pixel coordinates as well as node connection. The order of the morphology indicators are arranged in a different way.
It starts from the the bottom-left corner and first march along horizontal direction and then up. 
For example the morphology in previous example reads:

    TITLE ="t=0.000000"                                           | First three lines are overall mesh information:
    VARIABLES = "x"  "y"      "MS"                                | : 3 elements, 8 nodes, the first and second columns are coordinates, 
    ZONE N=8, E=3, F=FEPOINT, ET=QUADRILATERAL                    | last column is morphology
    0.00000  0.000000 1                                           <------First part stores the x , y, morphology indicator
    0.33333  0.000000 0
    0.00000  0.333333 1
    0.33333  0.333333 1
    0.00000  0.666666 0
    0.33333  0.666666 0
    0.00000  0.999999 1
    0.33333  0.999999 0
    1 2 4 3                                                       <------ Second part stores the node connection information.
    3 4 6 5                                                       !!!!!!!!NOTE:this second part is not mandatory. But it's good 
    5 6 8 7                                                       !!!!!!!!    to have it so you can plot&check your morphology.

NOTE: the coordinates don't need to be accurate because our code will calculate all the mesh information based on the the config file (see next subsection) - but it is important to keep it accurate so you can check the morphology before you run the simulations. It is suggested that for the first few tests, you check the morphology output in the result file to gain more confidence with the morphology reading.

### Config file

The config file stores all material properties and mesh information for calculation.
The config file has a fixed name: `config.txt`. An example config file is explained line by line below.

---------------------example of config file-----------------
```python
# if ifBoxGrid == True, you are using a mesh generated by talyFEM. Value option: True or False.
ifBoxGrid = True
# If "ifBoxGrid" is False, then the code grab an mesh generate by GMSH, otherwise the following line is ignored.
inputFilenameGrid = "mesh1.msh"
# Dimension of the physical problem. Value option: 2 for 2D 3 for 3D.
nsd = 2
# Keep the next two untouched.
nOfTS = 1
ifDD = False 
# There are two types of output file. If the output extention is ".h5" then the output is in binary format (save at least 50%) disk space.
# if the extention is ".plt", then the output is in ascii format. You need paraview to plot the binary format. "plt" format result can
# be open using Tecplot and paraview.
outputExtension = ".plt"
# Choices for finite element basis function: linear, quadratic or cubic
basisFunction = "linear"
# basisRelativeOrder determin theNumber of Gauss points for integration: if value is 0, then for linear basis function 4 Gauss points 
# are used and for quadratic 9 are used for cubic 16 are used; if 1 is used, for linear, 9 Gauss points are used, 16 for quadratic, 
# 25 for cubic; if -1 is used, then only 1 Gauss point is used for all kinds of elements. Value of 0 is recommended.
basisRelativeOrder = 1
# Keep dt unchanged, anyway we are dealing with a steady state problem.
dt = 0.001
#typeOfIC = 0
typeOfIC = 1
# Mobilities of electron (n), hole (p) and exciton (x). Unit: (m^2/Vs)
mu_n = 2e-7
mu_p = 1.5e-7
mu_x = 3.9e-9
# for example, electron mobiity in donor naterial will be mu_n*muRatio
# Mobility coefficient of particles as minority charge carrieres
muRatio = 0.01
# Realtive dielectric constants of donor (D) and acceptor (A)
eps_A = 3.9
eps_D = 3.
#potential difference between LUMO of acceptor and HOMO of donor, unit: eV.
E_g = 1.1
# Thickness of interfaical area, in nanometer
interfaceThk = 2
# Thickness of the domain, in meter
Height = 100e-9
# temperature, in K
T = 300
# Estimated open circuit voltage, so the calculation stops here
Voc = 0.69
# e-h separation distance, in meter
a = 1.8e-9
# lifetime of exciton (unit: s)
tau_x = 1e-6
# It is difficult to include 100% interfacial potential, so for most of the time, only a small fraction is added.
# use potentialFactor = 0.05 to added 5% real interfacial potential or 0 to neglected interfacial potential.
# It is suggested to keep it <0.08. 
potentialFactor = 0
# vapp_step (in eV) is the applied voltage step for Jv curve. At the turn of the JV curve, we need small steps to capture fill factor
vapp_step = 0.02
# vapp_ratio: Ratio of large applied voltage to fine calculation. At lower Vapp, we can use large step to save calculation time
# So at lower Vapp, the vapp_step really is vapp_step = vapp_step*vapp_ratio
# For example, with "vapp_step = 0.02, vapp_ratio = 5", you will start with Vapp step of 0.1eV(=0.02*5) 
# and when the JV curve change rapidly your Vapp step is 0.02eV
vapp_ratio = 5
# List the names of the parameters you are interested in and want to output. If not specified, then all available variables will be output.
# The most useful variables are [electricPotential, electron density, hole density, Exciton density, y-compoent electron current density
# y-component of hole current density, Dissociaiton rate, recombination rate, exiton generation rate, morphology] and they are listed below:
varListOutput = "[phi, n, p, Exciton, Jn_y, Jp_y, Dis, R, G, MS]"
# Define the finite element mesh with mesh size along each direction. If you are using quadratic element, these numbers have to be multiples 
# of 2; if you are using cubic element, these number have to be multiples of 3.
# Without weak boundary, (see below for "weak_bc"), for linear basis functionthe recommented mesh setup is 2elements per nm along x-direction and 500 along y-direction
Nelemx =100
Nelemy =394
# Mesh clustering option: if you want clustered mesh, beta=0.02 is recommended, otherwise use 0 for even mesh.
# Smaller beta value means more clustered mesh.
beta = 0.02
# Morphology file name
MSfile = "sawtoothwavy_80x86-2-40_trial3_test.bits"
# Value of (pixel number -1) in each direction. These are the numbers related to you input morphology.
# e.g, the coresponding values for previous example morphology, elmnx=1 and elmny=3
elmnx = 79
elmny = 85
# Exciton generation rate (unit: m^-3 s^-1)
Gx = 1e28
# if isPeriodic == True, then for the 2D case, your structure should be periodic structure along width direction.
# currently only tested with 2D morphology. If 'false', then the left and right boundaries has zero flux of all species.
isPeriodic = False
# Adding hole/electron transporting layers. This number is the ratio to the thickness of active layer, if not zero, you have to 
# set the thickness of the whole domain to be Height=(2*transportLayerThickness+1)*Height. This is to minimize the error
# caused by wrong material contact: donor contact with cathode and acceptor contact with anode. If there is not wrong material
# contact, you can set it to be zero. 
transport_layer_frac_thickness = 0.08  
# if weak_bc = 1 then weak boundary is used for electrode contacts, otherwise Neumann boundary condition is used.
# For now leave it as '0'
weak_bc = 0
# for 3phase model, set this to 1, otherwise, 0
threephase = 0
# if the morphology is not stored in binary 0-1 format, but a continuously change range, you can set the threshold to separate into to distinct phases
threshold = 0.5

```

--------------------the end of the example config file------------      


## Running the Program

Once you have the executable file "DD", the config file and a morphology file, you are ready to run the program. Put all the three into the same folder and run with (from a terminal, within the same folder where the input files are):

```bash
mpirun -np 16 ./DD  -pc_type lu -pc_factor_mat_solver_package mumps -snes_monitor -snes_converged_reason -snes_type newtonls -snes_rtol 1e-6
```

Here, `16` is the number of cores you want to use to run the simulation. You can change it to the number you want. The rest should stay unchanged.

Whenever you rerun, it is better to clean all existing result files.   


## Output Files

### text file: 'JVcurve_'+mophology_file_name

This file stores the JV curve. The first column is the applied voltage, the second column is the corresponding hole current density flows out of bottom surface, and the third column is the electron current density collected from the top surface. The last column is the output current density you should use when you report JV curves.

### text file: 'drg.txt'
It stores the domain integral of dissociation, recombination and generation (3 values and saved in this order) at short circuit condition.
The integral is done on a physical area of [active layer thickness] times [1cm]. In this case you can easily check the simulation result by checking if Jsc is close to Dissocation-Recombination (D-R).


### Contour result files

In the config file, you should specify the output variables that you are interested in (see [varListOutput](#markdown-Config%20file)). This will reduce the storage space needed for the contour files.

The contour result files are named following the convention:

```
result_"+appliedVoltage+"_"+inputMorphologyName+fileExtension
```

The most convenient way of viewing the results is using Tecplot. But if you are output as binary file, you have to use Paraview or Visit.

NOTE: when you plot the result file, you should notice that the coordinates of the domain is non-dimensionalized. The thickness of the domain is always 1. The unit of electron, hole, exciton density is `(m^-3 s^-1)`. The unit of the current density of electron and hole is `(mA cm^-2)`.
