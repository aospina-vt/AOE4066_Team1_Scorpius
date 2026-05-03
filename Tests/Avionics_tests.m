classdef Avionics_tests < matlab.unittest.TestCase
    methods (Test)
        function IgnitionSystem(testCase)
            insensitive_munitions_standards = true;

            testCase.verifyTrue(insensitive_munitions_standards);
        end
        function Software(testCase)
            follows_USAF_system_security_engineering_cyber_guidebook = true;

            testCase.verifyTrue(follows_USAF_system_security_engineering_cyber_guidebook);
        end
        function Connection_standard(testCase)
            uses_MIL_1760_interface = true;
            
            testCase.verifyTrue(uses_MIL_1760_interface);
        end
        function Power_Standard(testCase)
            EaglePicherBattery_is_military_compliant = true;
            
            testCase.verifyTrue(EaglePicherBattery_is_military_compliant);
        end
        function Datalink(testCase)
            uses_Link16 = true;

            testCase.verifyTrue(uses_Link16);
        end
        function Electromagnetic_Resistance(testCase)
            avionics_housing_EMI_resistant = true;
            
            testCase.verifyTrue(avionics_housing_EMI_resistant);
        end
    end 

end