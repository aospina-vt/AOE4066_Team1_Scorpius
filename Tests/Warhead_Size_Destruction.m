classdef Warhead_Size_Destruction < matlab.unittest.TestCase
    % Ensures the amount of warhead matterial is greater than 30 lbs, which
    % is the minimum needed to acheive the blastover pressure as described
    % in the explosive calc. 
    
    methods (Test)
        function warheadDestructionSize(testCase)
            [~,~,explosive_weight] = Weights();
            
            testCase.verifyGreaterThan(explosive_weight, 30);
        end
    end

end