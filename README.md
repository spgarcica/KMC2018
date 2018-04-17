# Kinetic Montecarlo MASM 2018 #

## Getting Started 

### Prerequisites

The following packages are mandatory to compile and run the program

```
make gcc
```

To plot the graphics python 2 and the following packages are required:

```
matplotlib numpy
```

### Compilation

Execute the following command to compilate the source code:

```
make compilation
```

### Data harvest

Modify the file **param.dat** to define the parameters of the simulation, then execute:

```
make datum
```

Two files are generated:

* **traj.xyz** contains the trajectory in a xyz format to visualize it with wxMacmolplt or vmd
* **ener.dat** contains the different energies of the system in every 

### Data plot

To plot the Energy and polarization graphics execute:

```
make plot
```

### Clean

To clean the data generated on the compilation execute:

```
make clean
```

## Author

* Sergio P. G.
