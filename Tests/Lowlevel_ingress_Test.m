classdef Lowlevel_ingress_Test < matlab.unittest.TestCase
    % The missile does NOT use a low level ingress in its profile, so this
    % test is false. 
    
    methods (Test)
        function lowlevel_ingress(testCase)
            low_level_ingress = false;

            testCase.verifyTrue(low_level_ingress);
        end
       
    end 

end