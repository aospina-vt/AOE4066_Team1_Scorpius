classdef LaunchBox_Size_Test < matlab.unittest.TestCase
    methods (Test)
        function internalAir(testCase)
            [L,W,H] = Dimensions();
            testCase.verifyLessThan(L, 144)
            testCase.verifyLessThan(W, 15)
            testCase.verifyLessThan(H, 15)
        end
    end

end