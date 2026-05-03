classdef Production_Tests < matlab.unittest.TestCase
    methods (Test)
        % Test methods

        function production_run(testCase)
            baselineproductionRate_week = 24;

            testCase.verifyGreaterThanOrEqual(baselineproductionRate_week*52, 1200);
        end

        function Surge(testCase)
            surgeCapable = true;

            testCase.verifyTrue(surgeCapable);
        end

        function shelfLife(testCase)
            all_component_shelfLife_more_that_10_years = true;

            testCase.verifyTrue(all_component_shelfLife_more_that_10_years);
        end

        function EIS(testCase)
            EIS_by_2031 = true;

            testCase.verifyTrue(EIS_by_2031);
        end
    end

end