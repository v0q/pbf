# Position Based Fluids (OOPG Assignment)
Position Based Fluids CPU (& rough GPU) implementation<br />
Position Based Fluids by Miles Macklin & Mathias Müller<br />
http://mmacklin.com/pbf_sig_preprint.pdf
https://github.com/NCCA/oopg1516-v0q

## Installation instructions:

# "Special" requirements:
OpenMP
GCC 5+
Might have to modify .pro file if something's not found

make clean<br />
qmake (-spec xxx)<br />
make<br />
./pbf<br />
<br />
Particle count can be modified by editing the for loop values inside the init() function in<br />
src/FluidSystem.cpp, due to the focus on the GPU implementation towards the end, no separate<br />
configuration file was implemented so the modifications to the simulation have to be done<br />
withing code.
<br />
Areas to modify to change the simulation<br />
FluidSystem.cpp<br />
    Constructor : BoundingBox min and max values and solver iterations<br />
    init() :  For loop max values to increase/decrease the particle count (note the initial positions relative to the bounding box size)<br />
    execute() : Time step in line 76 and m_bb.m_maxx value to match if bounding box values are changed
FluidSolver.cpp<br />
    Constructor : Solver parameters (everything except gravity and kernel constants)<br />
<br />

# Instructions:

Keyboard:<br />
<br />
1 - Toggle the simulation on/off<br />
2 - Toggle the wave machine on/off<br />
Escape - Exit the program<br />
