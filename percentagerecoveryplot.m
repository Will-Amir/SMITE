%%%% Percentage recovery plots
parameters=readtable("runoutput.csv");
maximums = [max(parameters.temperature) max(parameters.RCO2) max(parameters.biomass_tot)];
disp(maximums)