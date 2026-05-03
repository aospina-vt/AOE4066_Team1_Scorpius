classdef hit_the_target < matlab.unittest.TestCase
    methods (Test)
        function cep_REQ(testCase)
            MissDistance;
            CEP = distance * 1.744; % Converts mean miss-distance to CEP radius;
            testCase.verifyLessThan(CEP,5);
        end
    end

end