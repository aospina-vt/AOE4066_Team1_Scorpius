function T_LE = leadingEdgeTemp(Mach, Altitude, cp, T_init, r, z, time)
    T = atmosia(Altitude)

    k = 3.58e-6 * (717/(T0 + 225)) * (T0/492)^1.5
    h = k*Nnu;

    T_LE = T_r - (T_r - T_init).*(exp(-h*time/(cp*z)))

end

