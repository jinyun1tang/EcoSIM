# EcoSIM

A biogeochemical modeling library spin off the ecosys model. 

## Building

make config CC=icc CXX=icpc FC=ifort

make install CC=icc CXX=icpc FC=ifort

Find the executable under ./local/bin

## Examples
Examples can be found in examples/run_dir. 
To run the example, for instance, dryland, goes to dryland, then type

../../../local/bin/ecosys.x dryland.namelist 

and output can be found at dryland_maize_outputs.
