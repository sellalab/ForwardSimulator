# ForwardSimulator
These is the code for the forward simulator from Simons et al (2014).
To complie: g++ main.cpp population.cpp BRand.cpp
Program expects 3 parameters - number of runs, selection coefficient, dominance coefficient (h=0 recessive, h=0.5 additive).
The program writes three result files. One for summary statistics, one with joint frequency spectrums for both populations and one for the sampled joined frequency spectrums for both populations.

Branches provide variations on the main program:
Amorim2017 - version used in Amorim et al (2017).

