clear; clc;

possible_weights = 1:1:2000;

jp10_rho = 940; %kg/m^3
htp_rho = 1440; %kg/m^3

total_weight = possible_weights./2.205; %kg

mass_ox = (8/9) .* total_weight;
mass_fuel = (1/9) .* total_weight;

vol_ox = mass_ox ./ htp_rho;
vol_fuel = mass_fuel ./ jp10_rho;

total_vol = vol_ox+vol_fuel;

total_vol_in = total_vol.*61020;

plot(possible_weights,total_vol_in)

