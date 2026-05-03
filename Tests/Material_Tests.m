classdef Material_Tests < matlab.unittest.TestCase
    % Ensures the values obtained throguh ABAQUS for the strucutre and ANSYS for the heating is
    % compliant with the matterial properties. 
    
    methods (Test)

        function FactorSafety(testCase)
            results = StainlessSteel_Main();

            testCase.verifyGreaterThanOrEqual(results.factorSafety, 1.5)
        end

        function UltimateMissionLoad(testCase)
            results = StainlessSteel_Main();

            testCase.verifyFalse(results.failsUltimateLoad)
        end

        function ExtremeTemp(testCase)
            results = StainlessSteel_Main();

            tempLimit_K = 1073;

            testCase.verifyLessThanOrEqual(results.maxTemp_K, tempLimit_K)
        end

        function HeatingEffects(testCase)
            results = StainlessSteel_Main();

            tempLimit_K = 1073; 

            testCase.verifyLessThanOrEqual(results.maxCumulativeTemp_K, tempLimit_K)
        end

        function FatigueLife(testCase)
            results = StainlessSteel_Main();

            testCase.verifyGreaterThanOrEqual(results.fatigueLife_cycles, ...
                                             results.requiredLife_cycles)
        end

        function ProofPressure(testCase)
            results = StainlessSteel_Main();

            testCase.verifyGreaterThanOrEqual(results.actualProofPressure_psi, ...
                                             results.requiredProofPressure_psi)
        end

        function MinimumProofPressure(testCase)
            results = StainlessSteel_Main();

            testCase.verifyGreaterThanOrEqual(results.minimumProofPressure_psi, ...
                                             1.1 * results.maxOperatingPressure_psi)
        end

        function Stiffness(testCase)
            results = StainlessSteel_Main();

            testCase.verifyLessThanOrEqual(results.maxDeflection_in, ...
                                          results.allowableDeflection_in)
        end

        function LoadingWithCracks(testCase)
            results = StainlessSteel_Main();

            testCase.verifyLessThanOrEqual(results.crackedStress_ksi, ...
                                          results.allowableCrackedStress_ksi)
        end

        function FractureToughness(testCase)
            results = StainlessSteel_Main();

            testCase.verifyGreaterThanOrEqual(results.fractureToughness_ksi_sqrtin, ...
                                             0.60 * results.planeStrainToughness_ksi_sqrtin)
        end

        function FabricationReliability(testCase)
            results = StainlessSteel_Main();

            testCase.verifyTrue(results.fabricationWithinTolerance)
            testCase.verifyTrue(results.fabricationReliable)
        end

        function InternalPressure(testCase)
            results = StainlessSteel_Main();

            testCase.verifyLessThanOrEqual(results.internalPressure_psi, ...
                                          results.safeInternalPressure_psi)
        end

    end
end