# An-R-Based-Framework-for-Implementing-Large-Scale-Spatial-Models-of-Infectious-Diseases
Modelling the spread of infectious diseases in multiple populations

Demonstrated at the womENcourage conference 2015, Uppsala, Sweden

For use with other populations and sub-populations, shapefiles for those areas are needed. The contact matrix also needs to be changed. 
This matrix needs a contact rate for each population, but also for between each population. In this example it's 40x40 - 1600 values.

graphPlotMatrices shows all compartments, but depending on the number of subpopulations, it can be difficult to seperate them all.
The Susceptible stocks are in brown/yellow colour, Exposed in green, Infected in blue/purple, and Recovered in Pink.

mapPlotMatrices can show any compartment, here it is set for Infected. By dividing Infected by the Initial Population we can measure
the proportion of people infected. The code saves the images for each map in the folder specified, and in this example, creates 5001
images. Using software to transform them into videos makes it easier to visualise the spread of the infection.
