(experimental and learning by doing)

Overview:

The BOINC Project "Laser Induced Electron Emission" is basically a Solver for
the one-dimensional Schroedinger equation with electrons sitting in a simplified
solid state potential and getting nudged into freedom by a Laser field.

The main-app "app_opti_wrapp.cpp" is a behemoth wrapper which assembles and
executes various modules as requested in the allmighty config-file 
"liee_parameter.xml". The central modules are the potential and the solver, while
data from the simulation can be gathered by observers. The observers shall write
to outfiles, which get packed into an archive together with the log to be send
home to the BOINC-server.

On the BOINC-server runs the liee_scheduler with an optimization routine.
The idea and goal is to let the scheduler vary certain parameters (for instance
Laser-frequency). A single result from an observer (for instance total
tunnel current) is marked as the objective function. The scheduler will then
automatically generate work-units with "liee_parameter.xml" altered such in such
a way, that the optimization algorithm can progress in finding an extremum.

Currently implemented optimization algorithms are the straight forward 
Downhill-Simplex and a stochastic particle swarm optimizer (both not well 
tested). The tricky part for the optimization is to account for the 
asynchronous and delayed execution on those remote BOINC clients.