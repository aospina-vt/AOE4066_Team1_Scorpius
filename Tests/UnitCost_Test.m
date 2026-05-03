classdef UnitCost_Test < matlab.unittest.TestCase
    % Tests flyaway cost form the cost-curve is less than the RFP req of
    % 500,000. 
    
    methods (Test)
        function FlyAwayCost(testCase)
            cost = cost_tests();
            
            testCase.verifyLessThan(cost, 500000);
        end
    end

end