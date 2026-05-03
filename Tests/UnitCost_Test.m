classdef UnitCost_Test < matlab.unittest.TestCase
    methods (Test)
        function FlyAwayCost(testCase)
            cost = cost_tests();
            
            testCase.verifyLessThan(cost, 500000);
        end
    end

end