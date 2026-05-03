classdef LaunchBox_Size_Test < matlab.unittest.TestCase
    % Ensures our dimensions hit the requriments for the F-35, and thus all
    % launch boxes. 
    
    methods (Test)
        function internalAir(testCase)
            [L,W,H] = Dimensions_Boxed();
            testCase.verifyLessThanOrEqual(L, 144)
            testCase.verifyLessThanOrEqual(W, 15)
            testCase.verifyLessThanOrEqual(H, 15)
        end
    end

end