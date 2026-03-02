%% ==================== OPTIMIZED MISSILE TRAJECTORY ====================

close all; clear; clc;

%% ==================== PARAMETERS ====================
params.rollAng = 0;
params.AOA = 0;
params.L_fuse = 11.726;
params.S_front = 1.34646;
params.Diameter = 1.14102;
params.lnd = 3.4833/params.Diameter; 
params.ld = params.L_fuse/params.Diameter;
params.width  = (16/12)/2;
params.height = (12.27/12)/2;
d_nozzleExit = params.Diameter; 
params.A_nozzleExit = pi * d_nozzleExit^2 / 4;  

Isp = 290;       % s
T_max = 2500;    % lbf
g = 32.2;        % ft/s^2
W_0 = 1522;      % lb
m_0 = W_0/g;     % slugs

%% ==================== ATMOSPHERE TABLE ====================
h_grid = 0:100:120000;  % ft
a_grid = zeros(size(h_grid));
for ii = 1:length(h_grid)
    [~, a_tmp, ~, ~] = atmoscoesa(h_grid(ii)*0.3048);
    a_grid(ii) = a_tmp * 3.28084;  % m/s -> ft/s
end

dh = h_grid(2)-h_grid(1);  % step size for direct indexing

%% ==================== PRE-COMPUTED CONSTANTS ====================
nozzle_ratio = params.A_nozzleExit / params.S_front;
width_ratio = params.width / params.height;

%% ==================== TARGET ALTITUDES ====================
target_alts = 30000:5000:100000;
kick_angles = 4:0.05:8;
kick_alts   = 100:100:10000;
total_time = 1500;

final_results = nan(length(target_alts), 9);

%% ==================== CLIMB OPTIMIZATION ====================

climb_histories_x = cell(length(target_alts),1);
climb_histories_h = cell(length(target_alts),1);

for a = 1:length(target_alts)
    target_alt = target_alts(a);
    min_fuel_at_this_alt = inf;
    best_solution = [];
    found_valid = false;
    
    fprintf('Optimizing for Target Altitude: %d ft...\n', target_alt);

    for k_ang = kick_angles
        for k_alt = kick_alts
            
            hc_temp_hist = [];
            xc_temp_hist = [];

            if k_alt > target_alt  % early rejection
                continue
            end
            
            % --- Reset Initial Conditions ---
            V = 0.1; gamma = 90; h = 0; x = 0; m = m_0;
            dt = 0.5; not_kicked = true; AOA = 0;
            i = 0;

            % --- Run Climb ---
            while i <= total_time
                [L, D, ~] = computeForces(h, V, params, AOA);

                % Direct atmosphere lookup
                idx = real(max(1, min(length(h_grid), round(h/dh)+1)));
                a_sound = a_grid(idx);

                Mach = V / a_sound;

                % Thrust limit
                if Mach >= 3
                    T_req = (D + m*g*sind(gamma)) / cosd(AOA);
                    T = max(T_max*0.5, min(T_req, T_max));
                else
                    T = T_max;
                end

                % Physics derivatives (unchanged)
                V_dot = ((T*cosd(AOA) - D)/m) - g*sind(gamma);
                gamma_dot = 0;
                if V > 0.1
                    gamma_dot = rad2deg(((T*sind(AOA)+L) - (m*g*cosd(gamma))) / (m*V));
                end

                % Integrate
                V = V + V_dot*dt;
                if (h > k_alt && not_kicked)
                    gamma = gamma + gamma_dot*dt - k_ang;
                    not_kicked = false;
                    AOA = 3;
                else
                    gamma = gamma + gamma_dot*dt;
                end
                h = h + (V*sind(gamma))*dt;
                x = x + (V*cosd(gamma))*dt;
                m = m - (T/(Isp*g))*dt;

                hc_temp_hist(end + 1) = h;
                xc_temp_hist(end + 1) = x;

                if h > target_alt+5000 || h < -10
                    break
                end

                % Level-off success
                if h >= target_alt && gamma>-10 && gamma<10
                    fuel_used = (m_0 - m)*g;
                    if fuel_used < min_fuel_at_this_alt
                        min_fuel_at_this_alt = fuel_used;
                        best_solution = [k_ang, k_alt, V, gamma, h, x, m];
                        found_valid = true;
                        xc_best_hist = xc_temp_hist;
                        hc_best_hist = hc_temp_hist;
                    end
                    break
                end

                i = i + dt;
            end
        end
    end

    % Store results
    if found_valid
        final_results(a,:) = [target_alt, min_fuel_at_this_alt, best_solution];
        climb_histories_x{a} = xc_best_hist;
        climb_histories_h{a} = hc_best_hist;
        fprintf('  > Success! Min Fuel: %.2f lb (Kick: %.1f deg at %d ft)\n', ...
            min_fuel_at_this_alt, best_solution(1), best_solution(2));
    else
        final_results(a,:) = nan(1,9);
        final_results(a,1) = target_alt;
        fprintf('  > Failed to find level-off trajectory.\n');
    end
end

%% ==================== CLIMB FUEL PLOT ====================
figure('Color','w');
valid_data = final_results(~isnan(final_results(:,2)), :);
plot(valid_data(:,1)/1000, valid_data(:,2), '-o','LineWidth',2);
grid on; xlabel('Cruise Altitude (kft)'); ylabel('Climb Fuel (lb)');
title('Minimum Fuel to Reach Level Flight vs Altitude');

%% ==================== CRUISE ANALYSIS ====================
mission_dist = 500;               % nm
mission_dist_ft = mission_dist*6076.12;

mission_results = zeros(size(final_results,1),4);  % CruiseAlt, FuelClimb, FuelCruise, FuelTotal
glide_results = zeros(size(final_results,1),2); % Altitude , Glide nm


figure('Color','w'); hold on; grid on;
xlabel('Downrange Distance (nm)'); ylabel('Altitude (kft)');
title('Cruise Trajectories at Different Altitudes');

for a = 1:size(final_results,1)
    [~,a_sound,~,~] = atmoscoesa(target_alts(a)*0.3048);
    
    if (final_results(a,2) > 925), continue, end

    % --- Extract climb terminal state ---
    cruise_alt = final_results(a,1);
    fuel_climb = final_results(a,2);

    V     = final_results(a,5);  % Use V from best_solution
    gamma = final_results(a,6);
    h     = final_results(a,7);
    x     = final_results(a,8);
    m     = final_results(a,9);

    % Minimum structural + payload massuel
    max_fuel_mass = 925/g;
    m_empty = (m_0 - max_fuel_mass);

    ingress_dist_ft = mission_dist_ft*0.1;
    cruise_dist_ft = mission_dist_ft - (x);
    if cruise_dist_ft <= 0, continue, end

    dt = 0.05;  % larger step for speed
    x_start = x; 
    m_start = m;

    max_steps = ceil(cruise_dist_ft/(V*dt)) + 1000;
    x_hist = zeros(1,max_steps);
    h_hist = zeros(1,max_steps);
    vdot_hist = zeros(1,max_steps);
    T_hist = zeros(1,max_steps);
    D_hist = zeros(1,max_steps);
    V_hist = zeros(1,max_steps);

    step = 1;
    currAOA = 0;  % initial guess
    step_counter = 0;

    % --- Cruise loop ---
    while x <= (mission_dist_ft - ingress_dist_ft)
        step_counter = step_counter + 1;

        %Stop if fuel exhausted
        if m <= m_empty
            disp('Not enough fuel to complete segment');
            break; 
        end

        % --- Solve trim AOA only if lift deviates significantly ---
        if mod(step_counter,10) == 0
            trimFun = @(alpha) getLift(alpha,h,V,params)-m*g;
            try
                alpha_trim = fzero(trimFun, [-5,10]);
            catch
                L_min = getLift(h,V,params, -5);
                L_max = getLift(h,V,params, 10);
                if L_min > m*g
                    alpha_trim = -2;
                else
                    alpha_trim = 2;
                end
            end
            currAOA = alpha_trim;
        end

        % --- Compute forces ---
        [L,D,~] = computeForces(h,V,params,currAOA);

        D_hist(step) = D;

        % --- Thrust to balance drag ---

        if (V < 3 * a_sound)
            T = T_max;
        else
            T = max(T_max * 0.3,min(D/cosd(currAOA),T_max));
        end
        % --- Equations of motion ---
        V_dot = (T*cosd(currAOA)-D)/m;
        h_dot = 0;
        m_dot = -T/(Isp*g);

        % --- Integrate ---
        V = V + V_dot*dt;
        h = h + h_dot*dt;
        x = x + V*dt;
        m = m + m_dot*dt;

        % --- Store history ---
        x_hist(step) = x;
        h_hist(step) = h;
        vdot_hist(step) = V_dot;
        T_hist(step) = T;
        V_hist(step) = V;
        step = step + 1;

        % Break if cruise completed
        if (x - x_start) >= cruise_dist_ft
            disp('Mission segment complete');
            break; 
        end
    end

    % Trim unused preallocated space
    x_hist = x_hist(1:step-1);
    h_hist = h_hist(1:step-1);

    % --- Compute fuel ---
    fuel_cruise = (m_start - m)*g;
    fuel_total  = fuel_climb + fuel_cruise;

    mission_results(a,:) = [cruise_alt, fuel_climb, fuel_cruise, fuel_total];


    %% ==================== GLIDE SIMULATION ====================

    if isnan(final_results(a,2))
        continue
    end

    % --- Initialize glide state ---
    gamma = 0;       % start level
    T = 0;           % engine off
    dt_glide = 0.05;

    x_glide_start = x;
    h_glide = h;
    V_glide = V;
    m_glide = m;

    xg_hist = [];
    hg_hist = [];

    while h_glide > 0 && V_glide > 50

        % Atmosphere
        [~,a_sound,~,rho] = atmoscoesa(h_glide*0.3048);
        rho = rho*0.00194032;
        a_sound = a_sound*3.28084;

        % Trim AOA for lift = weight component
        trimFun = @(alpha) getLift(h_glide,V_glide,params,alpha) ...
                       - m_glide*g*cosd(gamma);
        try
            AOA_glide = fzero(trimFun, [-5 25]);
        catch
            AOA_glide = 2;
        end

        % Forces
        [L,D,~] = computeForces(h_glide,V_glide,params,AOA_glide);

        % Equations of motion
        V_dot = (-D/m_glide) - g*sind(gamma);
        gamma_dot = rad2deg((L - m_glide*g*cosd(gamma)) / (m_glide*V_glide));
        h_dot = V_glide*sind(gamma);
        x_dot = V_glide*cosd(gamma);

        % Integrate
        V_glide = V_glide + V_dot*dt_glide;
        gamma   = gamma + gamma_dot*dt_glide;
        h_glide = h_glide + h_dot*dt_glide;
        x       = x + x_dot*dt_glide;

        % Store
        xg_hist(end+1) = x;
        hg_hist(end+1) = h_glide;
    end

    x_climb = climb_histories_x{a};
    h_climb = climb_histories_h{a};

    glide_distance_ft = x - x_glide_start;
    glide_distance_nm = glide_distance_ft / 6076.12;

    glide_results(a,:) = [cruise_alt, glide_distance_nm];

    % --- Plot trajectory ---
    if a <= 60000
        plot([x_climb/6076.12, x_hist/6076.12 , xg_hist/6076.12], [h_climb/1000, h_hist/1000, hg_hist/1000], 'LineWidth',1.5);
    end
end

% --- Legend ---
legend(arrayfun(@(alt) sprintf('%dkft',alt/1000),mission_results(:,1),'UniformOutput',false));

%% ==================== CRUISE FUEL PLOTS ====================
valid = mission_results(:,4) > 0;

% Total mission fuel
figure('Color','w');
plot(mission_results(valid,1)/1000, mission_results(valid,4), '-o','LineWidth',2);
grid on; xlabel('Cruise Altitude (kft)'); ylabel('Total Fuel (lb)');
title('Total Mission Fuel vs Cruise Altitude');

% Cruise fuel only
figure('Color','w');
plot(mission_results(valid,1)/1000, mission_results(valid,3), '-s','LineWidth',2);
grid on; xlabel('Cruise Altitude (kft)'); ylabel('Cruise Fuel Only (lb)');
title('Cruise Fuel vs Cruise Altitude');

figure('Color','w');
plot(xg_hist/6072, ...
     hg_hist/1000, '-o','LineWidth',2);
grid on;
xlabel('Cruise Altitude (kft)');
ylabel('Glide Distance (nm)');
title('Glide Distance vs Cruise Altitude');

%plot(mission_results(:,1)/1000, D_values); % D_values = drag at cruise
%hold on; plot(mission_results(:,1)/1000, 0.5*T_max*ones(size(D_values)),'r--');
%% ==================== FUNCTIONS ====================
function L = getLift(alt,V,params,AOA)
    [L,~,~] = computeForces(alt,V,params,AOA);
end

function [F_lift,F_drag,LD_real] = computeForces(alt,V,params,AOA)
    [~,a,~,rho] = atmoscoesa(alt*0.3048);
    rho = rho*0.00194032;
    a   = a*3.28084;

    M = V./a;
    q = 0.5*rho.*V.^2;

    if M > 1
        CD0_wave = 3.6./((params.lnd.*(M-1))+3);
    else
        CD0_wave = 0;
    end
    CD0_base = (0.25./M).*(1 - params.A_nozzleExit/params.S_front);
    CD0_body = 0.053.*params.ld.*(M./(q*params.L_fuse)).^0.2;
    CD0 = (CD0_wave + CD0_base + CD0_body)*1.7;

    CN = computeCN(params.width,params.height,params.ld,params.rollAng,AOA);% *0.38;
    CD_total = CD0.*cosd(AOA)+CN.*sind(AOA);
    CL = (CN.*cosd(AOA)-CD0.*sind(AOA))*0.5;

    F_drag = CD_total.*q.*params.S_front;
    F_lift = CL.*q.*params.S_front;
    LD_real = CL./CD_total;
end

function CLa = getCLa(width,height,L_fuse,rollAng)
    dA = 0.5;
    Aplus  = 2 + dA;
    Aminus = 2 - dA;
    CN_plus  = computeCN(width,height,L_fuse,rollAng,Aplus);
    CN_minus = computeCN(width,height,L_fuse,rollAng,Aminus);
    CLa = (CN_plus-CN_minus)/deg2rad(2*dA);
end

function CN = computeCN(width,height,ld,rollAng,AOA)
    CN_roll = abs((width/height)*cosd(rollAng)) + abs((height/width)*sind(rollAng));
    CN = CN_roll*(abs(sind(2*AOA)*cosd(AOA/2))+2*ld*sind(AOA)^2);
end