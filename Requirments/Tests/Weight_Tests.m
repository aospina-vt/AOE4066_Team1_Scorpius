classdef Weight_Tests < matlab.unittest.TestCase
    methods (Test)
        function HIMARS_Fighter(testCase)
            [weight,~] = Weights();
            testCase.verifyLessThan(weight, 2000)
        end
        function submarine(testCase)
            [weight,~] = Weights();
            testCase.verifyLessThan(weight, 3500)
        end
        function mark41VLS(testCase)
            [weight,~] = Weights();
            testCase.verifyLessThan(weight, 5000)
        end
    end

end