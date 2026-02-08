classdef Cruise_Speed_Test < matlab.unittest.TestCase
    methods (Test)
        function cruise_speed_test_1(testCase)
            actualspeed = Cruise_Speed_Function();
            testCase.verifyGreaterThan(actualspeed, 3)
        end
    end
end