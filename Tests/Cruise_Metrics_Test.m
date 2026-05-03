classdef Cruise_Metrics_Test < matlab.unittest.TestCase
    % Uses the missile trajectory file to compute the cruise speed, range,
    % and air launch condition to ensure it conforms to reqs. 
    
    methods (Test)

        function Cruise_speed(testCase)
            [~, ~, ~, cruise_Mach, ~] = Scorpius_Trajectory();
            testCase.verifyGreaterThan(cruise_Mach, 3)
        end

        function MaxRange(testCase)
            [~, ~, ~, ~, maxRange_nm] = Scorpius_Trajectory();
            testCase.verifyGreaterThan(maxRange_nm, 500)
        end

        function AirLaunchCondition(testCase)
            [launch_mode, launch_Mach, launch_alt_ft, ~, ~] = Scorpius_Trajectory();
            testCase.verifyEqual(launch_mode, "AIR")
            testCase.verifyGreaterThanOrEqual(launch_Mach, 0.85)
            testCase.verifyEqual(launch_alt_ft, 30000)
        end

    end
end