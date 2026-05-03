classdef Cruise_Metrics_Test < matlab.unittest.TestCase

    methods (Test)

        function Cruise_speed(testCase)
            [cruise_Mach, ~] = Scorpius_Trajectory();
            testCase.verifyGreaterThan(cruise_Mach, 3)
        end

        function MaxRange(testCase)
            [~, maxRange_nm] = Scorpius_Trajectory();
            testCase.verifyGreaterThan(maxRange_nm, 500)
        end

    end
end