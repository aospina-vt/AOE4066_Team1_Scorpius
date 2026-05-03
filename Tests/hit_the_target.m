classdef hit_the_target < matlab.unittest.TestCase
    % Takes the mean miss-distance computed by the Simulink and converts it
    % to CEP by the conversion factor of 1.744. Confirms it is lower than 5
    % which, per the report, is the minimum CEP we need to hit HDBT targets
    % and naval vessles. 
    
    methods (Test)
        function cep_REQ(testCase)
            MissDistance;
            CEP = distance * 1.744; % Converts mean miss-distance to CEP radius;
            testCase.verifyLessThan(CEP,5);
        end
    end

end