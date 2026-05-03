%% SCORPIUS Missile - Cost Estimates

% Cost Model References: Tactical Missile Design {AIAA} Textbook, Pat Artis

clearvars; clc; close all;

W_Launch = 1333; % Launch Weight, lbs

L = 0.7; % Learning Curve Coefficient

C1000 = 6100 * W_Launch^0.758; % Unit Production Cost of Missile #1000

syms C1st x

% Define cost function
Cost = C1st * L^(log(x)/log(2));  % log base 2 using change of base

% Set x = 1000 and equate to C1000
eqn = subs(Cost, x, 1000) == C1000;

% Solve for C1st
C1st_sol = solve(eqn, C1st);

% Convert to numeric value
C1st_value = double(C1st_sol);

fprintf('Unit #1 Total Production Cost = $%0.2f\n', C1st_value)

syms x

C = C1st_value * L^(log(x)/log(2));

C600 = double(subs(C, x, 600));

fprintf('IOC (600 Already Produced) Total Unit Production Cost = $%0.2f\n', C600)

CYear1_nosurge = double(subs(C, x, 600 + 1200));

eqn = C == 500000;

x_sol = solve(eqn, x);

number_of_units_to_hit_500k = double(x_sol);

%% Cost Curve - No Surge - 1200/year only

% Create numeric cost function (faster than symbolic for plotting)
CostFunctionNoSurge = @(x) C1st_value .* L.^(log(x)./log(2));

% Production range
x_vals = 1:15000;

% Evaluate cost curve
C_vals = CostFunctionNoSurge(x_vals);

figure
plot(x_vals, C_vals, 'LineWidth', 2)
hold on
grid on

% Vertical reference lines
xline(600, '--', 'LineWidth', 1.5);
xline(1800, '--', 'LineWidth', 1.5);
xline(number_of_units_to_hit_500k, '--', 'LineWidth', 1.5);

yline(500000, 'r--', 'LineWidth', 2);

% Labels
text(600, C600 + 100000, ...
    sprintf('  IOC (600 Units Already Produced)  %,.0f', C600), ...
    'VerticalAlignment','bottom')

text(1800, CYear1_nosurge + 100000, ...
    sprintf('  Year 1 (1200 Units)  %,.0f', CYear1_nosurge), ...
    'VerticalAlignment','bottom')

text(number_of_units_to_hit_500k - 1100, 600000, ...
    sprintf('  %.0f Units @ $500k - Year 6', number_of_units_to_hit_500k), ...
    'VerticalAlignment','bottom')

text(12000, 600000, '  $500,000 Target', ...
    'VerticalAlignment','bottom', 'Color','red')

% Axes labels
xlabel('Number of Units Produced')
ylabel('Complete Unit Production Cost ($)')
title('Scorpius Complete Production Cost Per Unit - Constant Baseline Production')

% Format y-axis as dollars
ytickformat('$%,.0f')

ylim([0 0.4e7])

hold off

%% Surge Initial Two Years

CostFunction2 = @(x) C1st_value .* L.^(log(x)./log(2));

surge1yearquantity = 600 + 3600;
surge2yearquantity = 600 + 2 * (3600);

surgeLifeCyclequantity = surge2yearquantity + 8 * (1200);

% Production range
x_vals = 1:15000;

% Evaluate cost curve
C_vals2 = CostFunction2(x_vals);

CYear1_surge = CostFunction2(surge1yearquantity);
CYear2_surge = CostFunction2(surge2yearquantity);

figure
plot(x_vals, C_vals2, 'LineWidth', 2)
hold on
grid on

% Vertical reference lines
xline(600, '--', 'LineWidth', 1.5);
xline(surge1yearquantity, '--', 'LineWidth', 1.5);
xline(surge2yearquantity, '--', 'LineWidth', 1.5);

yline(500000, 'r--', 'LineWidth', 2);

% Labels
text(600, C600 + 100000, ...
    sprintf('  IOC (600 Units Already Produced)  %,.0f', C600), ...
    'VerticalAlignment','bottom')

text(surge1yearquantity, CYear1_nosurge + 100000, ...
    sprintf('  Year 1 (4200 Units)  %,.0f', CYear1_nosurge), ...
    'VerticalAlignment','bottom')

text(surge2yearquantity, CYear2_surge + 300000, ...
    ['  Year 2 (' num2str(surge2yearquantity) ' Units) -  $' num2str(CYear2_surge,'%.0f')], ...
    'VerticalAlignment','bottom')

text(12000, 600000, '  $500,000 Target', ...
    'VerticalAlignment','bottom', 'Color','red')

text(3000, 2500000, ...
    sprintf('  SURGE PERIOD  %,.0f', C600), ...
    'VerticalAlignment','bottom')

text(surge2yearquantity, 1500000, ...
    sprintf('  Surge Period Ends After 2 Years Post-IOC  %,.0f', C600), ...
    'VerticalAlignment','bottom')
% Axes labels
xlabel('Number of Units Produced')
ylabel('Complete Unit Production Cost ($)')
title('Scorpius Complete Production Cost Per Unit - 2-Year Initial Surge Production')

% Format y-axis as dollars
ytickformat('$%,.0f')

ylim([0 0.4e7])

hold off

fprintf('Production Cost After Surging For 1 Year = $%0.2f\n', CYear1_surge);
fprintf('Production Cost After Surging For 2 Years = $%0.2f\n', CYear2_surge);

%% Additional Cost Estimates

tSDD = 5;
CSDD = 20000000 * tSDD^1.90; % System Development & Demonstration Costs
fprintf('Scorpius System Development & Demonstration Cost = $%0.2f\n', CSDD);

PropSystemCost = 250000;
AvionicsControlPowerCost = 160000;
StructuralMaterialsCost = 2200;
AvionicsHousingCost = 15000;
Labor = 9300;

FlyawayCost = PropSystemCost + AvionicsControlPowerCost + StructuralMaterialsCost + AvionicsHousingCost + Labor;
fprintf('Scorpius Flyaway Cost (RFP Requirement) = $%0.2f\n', FlyawayCost);

AcquisitionCost = sum(CostFunction2(1:surgeLifeCyclequantity));

TotalLifeCycleCost = (AcquisitionCost + CSDD) / 0.9;

OperatingCost = 0.1 * TotalLifeCycleCost;

OverheadCost = 0.13 * OperatingCost;

fprintf('Scorpius Acquisition Cost (Over 10 Years) = $%0.2f\n', AcquisitionCost);

fprintf('Scorpius Total Lifecycle Cost = $%0.2f\n', TotalLifeCycleCost);

fprintf('Scorpius 10-Year Lifecycle Operating Cost = $%0.2f\n', OperatingCost);

fprintf('Scorpius 10-Year Lifecycle Overhead Cost = $%0.2f\n', OverheadCost);

%% DONT USE THIS ONE - Cost Curve with Discontinuities - Surge Initial 2 Years



% b = log10(L) / log10(2);
% 
% C3x = C .* 3^b; % Production cost when output is tripled
% 
% % nosurge = if surging never happened, regular = production returning to
% % baseline after a surge period, surge = tripled output
% 
% CYear1_surge = double(subs(C3x, x, 3600));
% 
% CYear2_surge = double(subs(C3x, x, (2 * 3600)));
% 
% CYear3_surge = double(subs(C3x, x, (3 * 3600)));
% 
% CYear3_regular = double(subs(C, x, 600  + (2 * 3600))); % Beginning of the year
% 
% CYear4_regular = double(subs(C, x, 600 + (3 * 3600))); % Beginning of the year
% 
% CYear3_nosurge = double(subs(C, x, 600 + 1200 + (2 * 1200)));
% 
% CYear4_nosurge = double(subs(C, x, 600 + 1200 + (3 * 1200)));
% 
% fprintf('Production Cost After Surging For 1 Year = $%0.2f\n', CYear1_surge);
% fprintf('Production Cost After Surging For 2 Years = $%0.2f\n', CYear2_surge);
% fprintf('Year 3 Production Cost After Surging For 2 Years = $%0.2f\n', CYear3_regular);
% % fprintf('Production Cost After Surging For 3 Years = $%0.2f\n', CYear3_surge);
% % fprintf('Year 4 Production Cost After Surging For 3 Years = $%0.2f\n', CYear4_regular);
% % fprintf('Year 4 Production Cost If No Surge = $%0.2f\n', CYear4_nosurge);
% 
% Year3diff = CYear3_nosurge - CYear3_regular;
% Year4diff = CYear4_nosurge - CYear4_regular;
% 
% SurgeCost = @(x) C1st_value .* L.^(log(x)./log(2)) .* 3^b;
% ResumeBaseline = @(x) C1st_value .* L.^(log(x)./log(2));
% 
% % Segment 1: Pre-IOC
% x_pre = 1:600;
% C_pre = CostFunctionNoSurge(x_pre);
% 
% % Segment 2: Post-IOC (surge period)
% x_surge = 600:15000;
% C_surge = SurgeCost(x_surge);
% 
% % Segment 3: If return to baseline output after two years
% x_regular = 7800:15000;
% C_regular = ResumeBaseline(x_regular);
% 
% 
% 
% figure
% plot(x_pre, C_pre, 'LineWidth', 2)
% hold on
% plot(x_surge, C_surge, 'LineWidth', 2)
% hold on
% plot(x_regular, C_regular, 'LineWidth', 2)
% 
% grid on
% 
% xline(600, '--', 'LineWidth', 1.5);
% xline(7800, '--', 'LineWidth', 1.5);
% yline(500000, 'r--', 'LineWidth', 2);
% 
% text(600, C600 + 100000, ...
%     sprintf('  IOC (600 Units Already Produced)  %,.0f', C600), ...
%     'VerticalAlignment','bottom')
% 
% text(5000, 550000, '  $500,000 Target', ...
%     'VerticalAlignment','bottom', 'Color','red')
% 
% text(3000, 3500000, ...
%     sprintf('  SURGE PERIOD  %,.0f', C600), ...
%     'VerticalAlignment','bottom')
% 
% text(7900, 1000000, ...
%     sprintf('  Surge Period Ends After 2 Years Post-IOC  %,.0f', C600), ...
%     'VerticalAlignment','bottom')
% 
% xlabel('Number of Units Produced')
% ylabel('Unit Production Cost ($)')
% title('Scorpius Unit Production Cost - Initial 2-Year Surge')
% 
% ytickformat('$%,.0f')
% ylim([0 0.4e7])
% 
% legend({'Pre-IOC (1–600 Units)', ...
%         'Post-IOC Surge (600–7800 Units)', ...
%         'Baseline Production after Surge'}, ...
%         'Location','northeast')
% 
% 
% hold off