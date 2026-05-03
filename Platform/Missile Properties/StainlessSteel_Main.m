function results = StainlessSteel_Main()

    %All values gotten from material properties or simulations in Abaqus
    %for the simulated loads and pressures in Ansys
    results.yieldStress_ksi = 31;
    results.ultimateStress_ksi = 73;
    results.tempLimit_K = 1073;


    results.maxStress_ksi = 4.5;
    results.factorSafety = results.yieldStress_ksi / results.maxStress_ksi;

    results.failsUltimateLoad = false;

    results.maxTemp_K = 893.151;
    results.maxCumulativeTemp_K = 900.56;

    results.fatigueLife_cycles = 135662;
    results.requiredLife_cycles = 100000;

    results.actualProofPressure_psi = 715.2;
    results.requiredProofPressure_psi = 500;

    results.maxOperatingPressure_psi = 335.32;
    results.minimumProofPressure_psi = 500;

    results.maxDeflection_in = 0.0034;
    results.allowableDeflection_in = 0.05;

    results.crackedStress_ksi = 18.922;
    results.allowableCrackedStress_ksi = 22;

    results.fractureToughness_ksi_sqrtin = 177.2;
    results.planeStrainToughness_ksi_sqrtin = 180;

    results.fabricationWithinTolerance = true;
    results.fabricationReliable = true;

    results.internalPressure_psi = 423;
    results.safeInternalPressure_psi = 500;

end