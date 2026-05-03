function dp_psi = peak_overpressure_psi_cased_mode(W_charge_lb, TNTeq_factor, W_casing_lb, r_ft, p0_psi, mode)
%PEAK_OVERPRESSURE_PSI_CASED_MODE
% Returns only peak overpressure Δp [psi].
%
% mode:
%   "BARE_SPHERE" -> shapeFactor=1.0  (no casing -> no change)
%   "SOURCE_EQ"   -> shapeFactor=1.19 (uses your excerpt's factor as-written)

    if nargin < 6 || isempty(mode)
        mode = "BARE_SPHERE";
    end
    mode = upper(string(mode));

    if any([W_charge_lb, TNTeq_factor, r_ft, p0_psi] <= 0) || W_casing_lb < 0
        error("Inputs must be > 0 (W_casing_lb may be 0).");
    end

    % Shape/reference factor
    switch mode
        case "BARE_SPHERE"
            shapeFactor = 1.0;
        case "SOURCE_EQ"
            shapeFactor = 1.19;
        otherwise
            error('mode must be "BARE_SPHERE" or "SOURCE_EQ".');
    end

    % TNT equivalent (peak-pressure basis)
    W_TNTeq_lb = TNTeq_factor .* W_charge_lb;

    % Map TNT-eq to the pentolite-fit "charge weight"
    F_ref = 1.42;
    c_ref_lb = W_TNTeq_lb ./ F_ref;

    % Casing ratio and M' (per your snippet)
    mr = 0;
    if W_casing_lb > 0
        mr = W_casing_lb ./ W_charge_lb;
    end

    if mr >= 1
        Mprime = 1.0;
    else
        Mprime = mr/2;   % includes mr=0 case
    end

    wprime_over_w = shapeFactor .* (1 + mr .* (1 - Mprime)) ./ (1 + mr);

    % Effective charge weight in correlation
    c_eff_lb = c_ref_lb .* wprime_over_w;

    % Bare-charge overpressure correlation (psi)
    z = r_ft ./ (c_eff_lb.^(1/3));
    x = z .* (p0_psi.^(1/3));

    dp_over_p0 = 37.95 ./ x + 154.9 ./ (x.^2) + 203.4 ./ (x.^3) + 403.9 ./ (x.^4);
    dp_psi = dp_over_p0 .* p0_psi;
end

o1 = peak_overpressure_psi_cased_mode(38.8, 1.42, 0, 10, 6.8, "BARE_SPHERE"); %missile book example
s1 = peak_overpressure_psi_cased_mode(10, 1.42, 190, 32.8084, 14.7, "SOURCE_EQ");
s2 = peak_overpressure_psi_cased_mode(20, 1.42, 180, 32.8084, 14.7, "SOURCE_EQ");
s3 = peak_overpressure_psi_cased_mode(30, 1.42, 170, 32.8084, 14.7, "SOURCE_EQ");
s4 = peak_overpressure_psi_cased_mode(40, 1.42, 160, 32.8084, 14.7, "SOURCE_EQ");
s5 = peak_overpressure_psi_cased_mode(50, 1.42, 150, 32.8084, 14.7, "SOURCE_EQ");
s6 = peak_overpressure_psi_cased_mode(60, 1.42, 140, 32.8084, 14.7, "SOURCE_EQ");
s7 = peak_overpressure_psi_cased_mode(70, 1.42, 130, 32.8084, 14.7, "SOURCE_EQ");

function V_total_in3 = warhead_volume_from_weights(W_charge_lb, W_casing_lb, rho_charge_lb_in3, rho_steel_lb_in3)
%WARHEAD_VOLUME_FROM_WEIGHTS_IN3 Total volume (in^3) from weights and densities in inch-units.
%
% Inputs:
%   W_charge_lb         charge weight [lbm]
%   W_casing_lb         casing weight [lbm]
%   rho_charge_lb_in3   charge density [lbm/in^3]
%   rho_steel_lb_in3    steel density  [lbm/in^3] (default = 490/1728)
%
% Output:
%   V_total_in3         total volume [in^3]

    if nargin < 4 || isempty(rho_steel_lb_in3)
        rho_steel_lb_in3 = 490/1728; % lbm/in^3 (typical steel)
    end

    if any([W_charge_lb, W_casing_lb, rho_charge_lb_in3, rho_steel_lb_in3] <= 0)
        error("All inputs must be > 0.");
    end

    V_charge_in3 = W_charge_lb ./ rho_charge_lb_in3
    V_steel_in3  = W_casing_lb ./ rho_steel_lb_in3

    V_total_in3  = V_charge_in3 + V_steel_in3;
end

rho_pentolite = 0.0596;  % if you have ~SG=1.65 -> lbm/ft^3 (example)
vol = warhead_volume_from_weights(50, 150, rho_pentolite);

function penetration = calc_pen(volume, diameter, V,rho_target,sigma_target)
% Inputs:

%   volume           volume of cylinder [inch3]
%   diameter         diameter of cylinder [inch]
%   V                strike velocity [ft/s]
%   rho_target = 4.02 % slug/ft^3 (concrete)  15.22 for HY-100
%   sigma_target 720000 psf (concrete)
rho_pen=15.19; % slug/ft^3 (metal)

d=diameter;
l  = (4*volume) ./ (pi * diameter.^2)
ld =l/d %should be higher than 2

penetration=(((l/d)-1)*((rho_pen/rho_target).^(1/2))+3.67*((rho_pen/rho_target).^(2/3))*(((rho_target*(V)^2)/sigma_target).^(1/3)))*(d/12);
%ft penetration

end

concretepenetration = calc_pen(1.367905766333379e+03, 9.5, 1350,4.02,720000)
shippenetration = calc_pen(1.367905766333379e+03, 9.5, 1350,15.22,14400000)

figure; hold on; grid on; axis equal;
xlabel('Explosive weight out of 200 lbs total warhead [lbs]');
ylabel('Blast overpressure experienced from 10m away [psi]');
title('Explosive weight vs Blast overpressure from 10m away');

% Straight track and predicted no-evasion end point
plot(10:10:70, [s1,s2,s3,s4,s5,s6,s7]);