function [Initial, Final, explosives, Strucutre, Avionics, Engine, Fuel_Ox] = Weights()
    % These weights were dervied from spreadsheets
    % LBS
    explosives = 50;
    casing = 150;
    
    Fuselage = 161.6;
    Fins = 38.17;
    Nose = 90.9;
    Strucutre = Fuselage + Fins + Nose;
    
    Avionics = 12.5;
    
    Engine = 41.67;
    Fuel_Ox = 74.5+595.24;
    Tank = 73;
    Nozzle = 45;
    Propulsion = Fuel_Ox + Tank + Nozzle + Engine;

    Initial = Strucutre + Avionics + Propulsion + explosives + casing;
    Final = Initial - Fuel_Ox;
end