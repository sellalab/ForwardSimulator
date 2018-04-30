# ForwardSimulator for Simons et al. (2018)

This version has been customized for Simons et al (2018).

Notable changes in this version: Demographic model changed to Schiffles and Durbin’s MSMC based inferred model. Mode of selection changed to underdominance.


### Compiling: 
g++ main.cpp population.cpp BRand.cpp
### Running:
Program expects 2 parameters - number of runs, selection coefficient (selection is always underdominant).  
The program writes one result fil with list of segregating mutation's age and frequency.  


## Branches
Other branches:
### Amorim2017 - version used in Amorim et al (2017). 
### Main - original version used in Simons et al (2014). 


