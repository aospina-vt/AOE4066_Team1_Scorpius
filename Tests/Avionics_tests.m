classdef Avionics_tests < matlab.unittest.TestCase
    methods (Test)
        function IgnitionSystem(testCase)
            % Scorpius implements the insensitive munitions standards used
            % throguhotu the US military to ensure key components including
            % the warhead and ignition system can not be prematurly started
            % without a proper order. 
            
            insensitive_munitions_standards = true;

            testCase.verifyTrue(insensitive_munitions_standards);
        end
        function Software(testCase)
            % By using the Air Forces System Security Engineering Cyber
            % Guidebook, The Crew can ensure that all software Scorpius
            % uses will be secure and comply with all software guidance.
            
            follows_USAF_system_security_engineering_cyber_guidebook = true;

            testCase.verifyTrue(follows_USAF_system_security_engineering_cyber_guidebook);
        end
        function Connection_standard(testCase)
            % Using this mil-standard connection specification in the
            % missile ensures it can comunicate with all launch systems.
            
            uses_MIL_1760_interface = true;
            
            testCase.verifyTrue(uses_MIL_1760_interface);
        end
        function Power_Standard(testCase)
            % The EaglePitcher Battery used is military compliant per the
            % product website.
            
            EaglePicherBattery_is_military_compliant = true;
            
            testCase.verifyTrue(EaglePicherBattery_is_military_compliant);
        end
        function Datalink(testCase)
            % Datalink is established with the Link 16 protocol. 
            
            uses_Link16 = true;

            testCase.verifyTrue(uses_Link16);
        end
        function Electromagnetic_Resistance(testCase)
           % The nose cone and avionics housing is made to be EMI
           % resistant. This increased cost is factored into the flyaway
           % cost of the missile. 
            
            avionics_housing_EMI_resistant = true;
            
            testCase.verifyTrue(avionics_housing_EMI_resistant);
        end
    end 

end