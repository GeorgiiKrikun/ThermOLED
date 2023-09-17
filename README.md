# ThermOLED

The 3D time-resolved simulation of organic light emitting diodes. Main simulation tool that I used to finish my PhD [dissertation](https://online.tugraz.at/tug_online/pl/ui/$ctx;lang=EN/wbAbs.showThesis?pThesisNr=72911&pOrgNr=2380) at TUGraz. 
It takes a description of multiple organic layers and outputs the steady state (if exists) for such configuration. Steady state includes: charge carrier distributions, temperature distribution, electric current, electric field distribution.
While the simulation is actually 3-dimensional, we found the most interesting results eventually come from one-dimensional distribution, hence the version here is customized for one dimension. If you want to start it for 3 dimensions, don't hesitate to contact me.

Simulation is written in Fortran90 and uses [OpenMPI](https://www.open-mpi.org/) for parallelisation on computational cluster.

## Compile
`make`

## Input
To run the simulation you need to first customize the layers. First number in input.dat corresponds to the number of layers. Then each layer description follows separated by newline. Each layer in turn is described by the space separated quantities that go as follows:
* pmu0 - hole mobility
* Width of hole density of states distribution. In kT units, T is assumed here to be 300K. During simulation, this value will get adjusted accordingly to temperature.
* LUMO - lowest unoccupied molecular orbital in eV units.
* nmu0 - electron mobility
* Width of electron density of states distribution. In kT units, T is assumed here to be 300K. During simulation, this value will get adjusted accordingly to temperature.
* HOMO - highest occupied molecular orbital in eV units.
* lattice - lattice constant. The average distance between molecules (SI units)
* epsR - relative electric permittivity
* rho - density of the layer (SI)
* cp - heat capacity of the layer (SI)
* kappa - thermal conductivity of the layer (SI)
* matbegin - where the layer starts (SI)
* matend - where the layer ends (SI)

Note that simulation is inherently unstable (it's idea was to capture the unstable nature of OLEDs) and can get into thermal runaway process quite quickly. Values you provide must therefore be at least realistic.

## Run
mpirun -np N_processes thermoled N_x N_y N_z N_cut_x N_cut_y N_cut_z V
Where 
1. N_x, N_y, N_z describe discretization along corresponding axes. I suggest for tests in 1D: 30 1 1
2. N_cut_x N_cut_y N_cut_z number of parrallel processing units in each direction. N_cut_x*N_cut_y*N_cut_z=N_processes
3. V voltage applied to the first layer (voltage on the other side is assumed to be 0.
