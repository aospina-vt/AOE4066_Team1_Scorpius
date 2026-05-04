function [best_method, results] = ClimbTradeStudy(target_alt, show_plots)
    % evaluate_climb_methods compares different climb strategies for Scorpius.
    % 
    % Inputs:
    %   target_alt - Target cruise altitude (e.g., 75000)
    %   show_plots - boolean (true/false) to display comparison plots
    %
    % Outputs:
    %   best_method - The name of the winning strategy (string)
    %   results     - A table containing raw data and calculated scores

    if nargin < 2, show_plots = true; end
    if nargin < 1, target_alt = 75000; end

    %% ==================== 1. PARAMETERS & CONFIG ====================
    const_angles = 30:5:90; 
    case_names = [string(const_angles), "Gravity Turn"];
    results = table(); 

    % Airframe Parameters
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
    n_max = 4; max_fuel = 670; m_empty = m_0 - (max_fuel/g);

    if show_plots
        fig = figure('Color','w','Position', [50, 50, 1200, 800]);
        subplot(2,2,1); hold on; grid on; xlabel('Range (mi)'); ylabel('Altitude (ft)');
        title('Trajectory Comparison: Smooth Level-off');
    end

    %% ==================== 2. ITERATIVE SIMULATION ====================
    for i = 1:length(case_names)
        current_case = case_names(i);
        
        % Initial State
        h = 0; x = 0; gamma = 90; m = m_0; V = 1.0; dt = 0.1;
        total_x = []; total_h = [];
        h_capture_start = target_alt * 0.85; 
        
        while h < target_alt + 100 && m > m_empty
            % --- 1. GUIDANCE LOGIC ---
            if current_case == "Gravity Turn"
                if h < 2000
                    gamma_cmd = 90; aoa_steer = 0;
                elseif h < 3500
                    gamma_cmd = 75; % Early pitch-over kick
                    aoa_steer = (gamma_cmd - gamma) * 0.5;
                else
                    aoa_steer = 0; % Natural gravity-driven curve
                end
            else
                target_gamma = double(current_case);
                gamma_cmd = (h < 2000)*90 + (h >= 2000)*target_gamma;
                aoa_steer = (gamma_cmd - gamma) * 0.5;
            end
            
            % --- 2. LEVEL-OFF CAPTURE (Vz RATE CONTROL) ---
            if h > h_capture_start
                h_err = target_alt - h;
                Vz_target = max(0, h_err * 0.15); 
                Vz_actual = V * sind(gamma);
                aoa_steer = (Vz_target - Vz_actual) * 0.02;
                aoa_steer = aoa_steer + (0 - gamma)*0.1; 
            end

            % --- 3. PHYSICS & CONSTRAINTS ---
            aoa_steer = max(-8, min(aoa_steer, 15)); 
            [L, D, ~] = computeForces(h, V, params, aoa_steer);
            
            if (L/(m*g)) > n_max
                aoa_steer = aoa_steer * (n_max / (L/(m*g)));
                [L, D, ~] = computeForces(h, V, params, aoa_steer);
            end
            
            T = T_max;
            V_dot = (T - D)/m - g*sind(gamma);
            g_dot = rad2deg((L - m*g*cosd(gamma))/(m*max(V,10)));
            
            V = V + V_dot*dt; gamma = gamma + g_dot*dt; 
            h = h + V*sind(gamma)*dt; x = x + V*cosd(gamma)*dt; 
            m = m - (T/(Isp*g))*dt;
            
            total_x(end+1) = x; total_h(end+1) = h;
            if h >= target_alt - 20 && abs(gamma) < 1.0, break; end
        end
        
        fuel_used = W_0 - (m*g);
        results = [results; table(current_case, fuel_used, x/5280, V, ...
            'VariableNames',{'Method','Fuel_lb','Range_mi','Final_Vel'})];
            
        if show_plots
            subplot(2,2,1);
            plot(total_x/5280, total_h, 'LineWidth', 2, 'DisplayName', char(current_case));
        end
    end

    %% ==================== 3. SCORING LOGIC ====================
    % Weighting: 40% Fuel, 25% Range, 35% Velocity
    results.Fuel_Score = (min(results.Fuel_lb) ./ results.Fuel_lb) * 100;
    results.Range_Score = (results.Range_mi ./ max(results.Range_mi)) * 100;
    results.Vel_Score = (results.Final_Vel ./ max(results.Final_Vel)) * 100;
    results.Total_Score = (0.4 * results.Fuel_Score) + (0.25 * results.Range_Score) + (0.35 * results.Vel_Score);

    [~, best_idx] = max(results.Total_Score);
    best_method = results.Method(best_idx);

    %% ==================== 4. PLOTTING (OPTIONAL) ====================
    if show_plots
        subplot(2,2,1); yline(target_alt, '--k', 'Cruise Altitude'); legend('Location','Best');
        
        subplot(2,2,2);
        bar(categorical(results.Method), results.Fuel_lb, 'FaceColor', [0.8 0.2 0.2]);
        grid on; ylabel('Fuel Burn (lb)'); title('Cost: Fuel Consumed');

        subplot(2,2,3);
        bar(categorical(results.Method), results.Final_Vel, 'FaceColor', [0.2 0.7 0.2]);
        grid on; ylabel('Final Velocity (ft/s)'); title('Benefit: Level-off Speed');

        subplot(2,2,4);
        b = bar(categorical(results.Method), results.Total_Score);
        b.FaceColor = 'flat';
        b.CData(best_idx,:) = [0.1 0.1 0.8];
        grid on; ylabel('Score (0-100)'); title('Mission Performance Score');
    end
end

%% ==================== LOCAL FUNCTIONS ====================
function [F_lift, F_drag, LD_real] = computeForces(alt, V, params, AOA)
    alt_m_capped = max(0, min(84852, alt));
    [~, a, ~, rho] = atmoscoesa(alt_m_capped * 0.3048);
    rho = rho * 0.00194032; a = a * 3.28084;      
    M = V./max(a, 1); q = 0.5 * rho .* V.^2;
    if M > 1.0
        CD0_wave = 3.6 ./ ((params.lnd .* (M - 1)) + 3);
        CD0_base = (1 - params.A_nozzleExit / params.S_front) * (0.25 / M);
    else
        CD0_wave = 0; CD0_base = (1 - params.A_nozzleExit / params.S_front) * (0.12 + 0.13 * M^2);
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