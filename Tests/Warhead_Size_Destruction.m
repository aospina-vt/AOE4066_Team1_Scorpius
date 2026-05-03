classdef Warhead_Size_Destruction < matlab.unittest.TestCase
    methods (Test)
        function warheadDestructionSize(testCase)
            [~,~,explosive_weight] = Weights();
            
            testCase.verifyGreaterThan(explosive_weight, 30);
        end
    end

end