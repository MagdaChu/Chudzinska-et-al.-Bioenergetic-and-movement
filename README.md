# Chudzinska-et-al.-Bioenergetic-and-movement
Supplementary material for Chudzinska et al. Combining bioenergetics and movement models to improve understanding of the population consequences of disturbance

The inputs are structured in four folders.

'Bioenergetic model' is the folder containing the bioenergetic models for the three species. See next section on how to use it

'3D graph SI' contains code and input file to generate Figure E1 from the Supplementary Information

'Movement' contains all input data and code to generate Figure 2 from the main paper

'Reports' contain two reports which may not be easy for a reader to find online

Code/Software

The bioenergetic models are written in R software. The code has the same structure for each species. SpeciesDEB_Control.R and SpeciesDEB_Params.R are the only two code needed user input. The SpeciesDEB_params.R file contain the list of all parameters need to run the models. The values in the codes are the values used in the final simulations. To see the list and meaning of each parameter, refer to the Supplementary Information. In SpeciesDEB_Control.R code, the user can define the settings and parameters for the disturbance scenario. Follow the commented lines in the code for details. This is the main script to start the simulations.
