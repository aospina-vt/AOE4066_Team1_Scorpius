function [Initial, Final, Warhead, Strucutre, Avionics, Engine, Fuel_Ox] = Weights()
    %LBS
    Warhead = 200;
    Strucutre = 291;
    Avionics = 19;
    Engine = 87;
    Fuel_Ox = 1400;

    Initial = Strucutre + Avionics + Engine + Fuel_Ox + Warhead;
    Final = Initial - Fuel_Ox;

end