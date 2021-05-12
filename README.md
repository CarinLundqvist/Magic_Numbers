# Magic_Numbers
Code used to model a few ultracold fermions in a harmonic oscillator. The code was used in the Bachelor Thesis "Patterns in interacting quantum gases" at Chalmers Univerisity of Technology, written by Viktor Bekassy, Olivia Heuts, Alex Lech, Carin Lundqvist and Erik Magnusson. All code was written by either Carin Lundqvist or Viktor Bekassy.

The building block of all scripts is the Quantum Harmonic Oscillator (QHO). It is necessary for all calculations. In the case of the Pauli crystals, it is sometimes necessary to use the 2D version.

The Pauli_Crystals folder contains a couple of different ways of simulating Pauli crystals. The probability densities use the single particle wave functions, wheras the algorithms used in Pauli_crystals.m requiers the full many-body wave function constructed in N_particles_1D and N_particles_2D. In order to see reall cool patterns in 2D, Pauli_distance.m is also required to perform the necessary image processing.

When perturbations are added, the tricky function is Wavefunction_corrected.m. A moderate understanding of perturbation theor is necessary to understand the code.
