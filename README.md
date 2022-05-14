# Sortase_Polymer_Gillespie
Repository to simulate Sortase-mediated transpeptidation using a Gillespie Algorithm

## Modelling the Reaction Dynamics of Sortase Transpeptidase
Herein the Gillespie algorithm is used for a few simple chemical reactions. We'll focus on Sortase Transpeptidase in particular: The simulation of a sortase reaction. 

## Motivation 
Making a simulation that mimics a sortase reaction to get a better understanding of the process. 
## Layout 
Start with a simple reaction to see how the sortase reaction is set up intailly. The Sortase A from S. aureus has been optimized to act as a transpeptidase to catalyze the formation of a peptide bond between the C-terminal LPXTG sequence of one protein to an N terminal G on another protein.  Sortase A is also able to catalyze the reverse reaction, whereby a ligated product is broken apart into its constitutive components.
Our interest in modeling this reaction stems from the fact that the efficient cross-linking reaction for sortase is dependent on many different factors, including the initial concentrations of reactants, the temperature of the sortase reaction, and the inherent properties of the domains of the reactants.  Experimentally, this means that the optimal conditions for sortase reactions have to be ascertained experimentally for each unique sortase reaction.  (Li, Jess. et al "Optimization of sortase A ligation for flexible engineering of complex protein systems" 2020).

How is the sortase reaction set up? \
R1: A + B -> C  reaction constant kf \
R2: C -> A + B  reaction constant kr \
Slowly begin to add more values and terms to create a more complex reaction.\
Finally we automate the process.\
We must repeat the same process in order to simulatethe system to a maximum time (`tmax`) or maximum number of reactions (`nmax`) of our choice. Creating a reapeat reaction. Plotting the simulation help get a visual representation of the sortase reaction when A,B,C are present. Running the now automated code also allows to comfirm the reaction is at a steady state as multple simulations are ran all giving similar results.

**Stoichiometrically Equivalent Ratio of Reactants Appears to Favor Forwards Reaction When Forwards and Reverse Reaction Rates Are Held the Same**


**Model Substrate Death**


**Model Enzyme Building Long Chains**


**Model Enzyme Death With Building Long Chains**


## Installation 
Anaconda3, this will inculde Jupyter Notebook and Spyder.
## Tests

## Instructions 


### References

The tutorial is heavily based on Gillespie's original 1977 [paper](http://wwwf.imperial.ac.uk/~nsjones/gillespie_1977.pdf) and we will often  quote text directly from the paper. 

The Wikipedia [page](https://en.wikipedia.org/wiki/Mathematical_model) on Mathematical Models.


### Acknowledgements

This tutorial is adapted from the one made by Dr. Karin Sasaki (Centre for Biological Modelling at the European Molecular Biology Laboratory).
code found [here] https://github.com/karinsasaki/gillespie-algorithm-python
