classdef LaunchBox_Size_Test < matlab.unittest.TestCase
    methods (Test)
        function internalAir(testCase)
            [L,W,H] = Dimensions_Boxed();
            testCase.verifyLessThanOrEqual(L, 144)
            testCase.verifyLessThanOrEqual(W, 15)
            testCase.verifyLessThanOrEqual(H, 15)
        end
    end

end