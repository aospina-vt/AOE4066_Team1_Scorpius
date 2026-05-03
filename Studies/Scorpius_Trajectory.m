%% ==================== INTEGRATED MISSILE TRAJECTORY: 75K FT ====================
close all; clear; clc;

%% ==================== 1. PARAMETERS & MISSION CONFIG ====================
launch_mode = "AIR"; 
target_alt  = 75000; 

params.width        = (14.8/12)/2;
params.height       = (12.27/12)/2;
params.rollAng      = 0;
params.L_fuse       = 11.726;
params.S_front      = 1.34646;
params.Diameter     = 2 * sqrt(params.width * params.height);
params.lnd          = 3.4833 / params.Diameter; 
params.ld           = params.L_fuse / params.Diameter;
params.S_wing       = 6;  
params.AR           = 4^2 / params.S_wing;  
params.A_nozzleExit = pi * params.Diameter^2 / 4;

Isp = 290; T_max = 2000; g = 32.2; W_0 = 1332; m_0 = W_0/g;  
n_max = 4; cruise_Mach = 3.5; max_fuel = 670; max_aoa = 40; 
m_empty = m_0 - (max_fuel/g);

min_throttle = 0.30; 
T_min = T_max * min_throttle; 

%% ==================== 2. MISSION SIMULATION ====================
total_x = []; total_h = []; total_V = []; total_gamma = []; total_m = [];

% -------------------- PHASE 1: CLIMB --------------------
h = (strcmp(launch_mode,"AIR")*30000); x = 0; gamma = (strcmp(launch_mode,"GROUND")*90); m = m_0;
V = (strcmp(launch_mode,"AIR")*0.8*1100) + (strcmp(launch_mode,"GROUND")*0); 
dt = 0.1;

m_start_climb = m;
fprintf('Phase 1: Climb...\n');
while h < target_alt && m > m_empty
    [L, D, ~] = computeForces(h, V, params, 3);
    T = T_max; 
    
    if strcmp(launch_mode, "AIR")
        n_aero = L/(m*g); 
        g_dot_max = rad2deg(max(0, n_max - n_aero)*g/max(V,1));
        g_dot = (gamma < 35)*min(g_dot_max, (35-gamma)/dt) + (h > target_alt*0.95)*(-gamma/2);
    else
        if h > 1000 && gamma > 80, gamma = 80; end
        g_dot = (V > 50)*rad2deg((L - m*g*cosd(gamma))/(m*max(V,10)));
        if h > target_alt*0.95, g_dot = 0; gamma = gamma*0.9; end
    end
    
    V_dot = (T - D)/m - g*sind(gamma);
    V = V + V_dot*dt; gamma = gamma + g_dot*dt; h = h + V*sind(gamma)*dt; x = x + V*cosd(gamma)*dt; 
    m = m - (T/(Isp*g))*dt;
    
    total_x(end+1)=x; total_h(end+1)=h; total_V(end+1)=V; total_gamma(end+1)=gamma; total_m(end+1)=m;
    if h > target_alt && abs(gamma) < 2, break; end
end
fuel_climb = (m_start_climb - m) * g;

% -------------------- PHASE 2: POWERED CRUISE (dh=0) --------------------
m_start_cruise = m;
fprintf('Phase 2: Powered Cruise...\n');
dt = 0.5; currAOA = 2; h = target_alt;
while m > m_empty
    [L, D, ~] = computeForces(h, V, params, currAOA);
    currAOA = currAOA + 0.2 * ((m*g - L) / (m*g)); 
    currAOA = max(0, min(currAOA, max_aoa));

    [~, a_at_h, ~, ~] = atmoscoesa(h * 0.3048);
    a_s = a_at_h * 3.28084;
    V_target = cruise_Mach * a_s;
    
    if V < V_target
        T = T_max; 
    else
        T = max(D, T_min); 
    end
    
    V_dot = (T - D)/m; 
    m = m - (T/(Isp*g))*dt; 
    V = V + V_dot*dt; x = x + V*dt;
    
    total_x(end+1)=x; total_h(end+1)=h; total_V(end+1)=V; total_gamma(end+1)=0; total_m(end+1)=m;
end
fuel_cruise = (m_start_cruise - m) * g;

% -------------------- PHASE 3: LEVEL UNPOWERED GLIDE (dh=0) --------------------
fprintf('Phase 3: Level Glide...\n');
while true
    [L, D, ~] = computeForces(h, V, params, currAOA);
    currAOA = currAOA + 0.2 * ((m*g - L) / (m*g)); 
    
    if currAOA >= max_aoa || L < m*g
        break; 
    end 
    
    V = V + (-D/m)*dt; 
    x = x + V*dt;
    
    total_x(end+1)=x; total_h(end+1)=h; total_V(end+1)=V; total_gamma(end+1)=0; total_m(end+1)=m;
end

% -------------------- PHASE 4: DIVE --------------------
fprintf('Phase 4: Glide Down...\n');
gamma_dive_target = -70; dt_glide = 0.02; gamma_g = 0; currAOA = 0; 

while h > 0 && V > 100
    [L, D, ~] = computeForces(h, V, params, currAOA);
    Weight = m * g;
    gamma_dot = ((L - Weight * cosd(gamma_g)) / (m * max(V,10))) * (180/pi);

    gamma_err = gamma_dive_target - gamma_g;
    currAOA = currAOA + 0.05 * gamma_err - 0.2 * gamma_dot; 
    currAOA = max(-10, min(currAOA, 15)); 

    gamma_g = gamma_g + gamma_dot * dt_glide;
    V_dot = (-D - Weight * sind(gamma_g)) / m;
    V = V + V_dot * dt_glide;
    h = h + (V * sind(gamma_g)) * dt_glide;
    x = x + (V * cosd(gamma_g)) * dt_glide;

    total_x(end+1)=x; total_h(end+1)=h; total_V(end+1)=V; total_gamma(end+1)=gamma_g; total_m(end+1)=m;
end

%% ==================== PLOTTING ====================
figure('Color','w','Position', [100, 100, 900, 700]);

% Find burnout point
% We look for where mass is essentially at the limit
burnout_idx = find(total_m <= m_empty + 0.001, 1);

subplot(2,1,1); 
plot(total_x/5280, total_h, 'b', 'LineWidth', 2); hold on;
if ~isempty(burnout_idx)
    x_burn = total_x(burnout_idx)/5280;
    h_burn = total_h(burnout_idx);
    plot(x_burn, h_burn, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    text(x_burn, h_burn + 3000, 'Burnout', 'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
grid on; ylabel('Altitude (ft)'); title('Missile Trajectory');

% Add Fuel Consumption Info Box
fuel_str = {sprintf('Total Fuel: %.1f lb', max_fuel), ...
            sprintf('Climb Fuel: %.1f lb', fuel_climb), ...
            sprintf('Cruise Fuel: %.1f lb', fuel_cruise), ...
            sprintf('Remaining: %.1f lb', (m - m_empty)*g)};
annotation('textbox', [0.15, 0.6, 0.2, 0.1], 'String', fuel_str, 'FitBoxToText', 'on', 'BackgroundColor', 'w');

subplot(2,1,2); 
plot(total_x/5280, total_V, 'r', 'LineWidth', 2); hold on;
if ~isempty(burnout_idx)
    plot(total_x(burnout_idx)/5280, total_V(burnout_idx), 'ko', 'MarkerFaceColor', 'k');
end
grid on; ylabel('Velocity (ft/s)'); xlabel('Range (miles)');
title('Velocity Profile');

%% ==================== FUNCTIONS ====================
function [F_lift, F_drag, LD_real] = computeForces(alt, V, params, AOA)
    alt_m_capped = max(0, min(84852, alt));
    [~, a, ~, rho] = atmoscoesa(alt_m_capped * 0.3048);
    rho = rho * 0.00194032; 
    a   = a * 3.28084;      
    M = V./max(a, 1); 
    q = 0.5 * rho .* V.^2;

    if M > 1.0
        CD0_wave = 3.6 ./ ((params.lnd .* (M - 1)) + 3);
        CD0_base = (1 - params.A_nozzleExit / params.S_front) * (0.25 / M);
    else
        CD0_wave = 0; 
        CD0_base = (1 - params.A_nozzleExit / params.S_front) * (0.12 + 0.13 * M^2);
    end
    CD0_body = 0.053 .* params.ld .* (M ./ max(q * params.L_fuse, 0.01)).^0.2;
    CD0 = (CD0_wave + CD0_base + CD0_body);
    CN_wing = computeCN_Wing(M, AOA, params.AR, params.S_wing, params.S_front);
    CN = computeCN(params.width, params.height, params.ld, params.rollAng, AOA) + CN_wing;
    CD_total = CD0 .* cosd(AOA) + CN .* sind(AOA);
    CL = CN .* cosd(AOA) - CD0 .* sind(AOA);
    F_drag = CD_total .* q .* params.S_front;
    F_lift = CL .* q .* params.S_front;
    LD_real = CL ./ max(CD_total, 0.001);
end

function CN = computeCN(width, height, ld, rollAng, AOA)
    CN_roll = abs((width/height)*cosd(rollAng)) + abs((height/width)*sind(rollAng));
    CN = CN_roll * (sind(2*AOA)*cosd(AOA/2) + 2*ld*sind(AOA)^2);
end

function CN_w = computeCN_Wing(M, AOA, AR, Sw, Sref)
    M_crit = sqrt(1 + (8 / (pi * AR))^2);
    alpha_rad = deg2rad(AOA);
    area_ratio = Sw / Sref;
    term_nonlinear = 2 * (sin(alpha_rad)^2);
    if M > M_crit
        term_linear = (4 * (sin(alpha_rad)*cos(alpha_rad))) / sqrt(max(M^2 - 1, 0.01));
    else
        term_linear = (pi * AR / 2) * (sin(alpha_rad)*cos(alpha_rad));
    end
    CN_w = (term_linear + term_nonlinear) * area_ratio;
end

maxRange_nm = max(total_x) / 6076.12; % Convert to nautical miles


