About
=====

LaDa is a python framework to control physics simulations, from DFT to
empirical pseudo-potentials to point ion electrostatics. It is designed to be
extremelly modular. It's goal is to provide the basic building blocks with
which any calculation can be constructred. We have used within many schemes.
In one, we would fire thousands of DFT calculations simultaneously, enumerating
all possible A2BX4 compounds within forty possible crystal structures. In a
more complex scheme, LaDa helped us to automatically determine for any crystal
structure the possible vacancy and substitutional point defects, as well as
their charge states, and then launch a series of individual DFT calculation
finally yielding their formation enthalpy. In third scheme, we
