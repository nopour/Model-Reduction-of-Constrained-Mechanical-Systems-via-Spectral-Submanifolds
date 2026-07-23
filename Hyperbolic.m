
clear all; close all; clc;

%% 0. Journal Plotting Settings
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultLineLineWidth', 1.5);

%  Color Palette for Multiple Orders
c_o3   = [0.4660, 0.6740, 0.1880]; % Green
c_o5   = [0.9290, 0.6940, 0.1250]; % Yellow/Orange
c_o7   = [0.4940, 0.1840, 0.5560]; % Purple
c_o11  = [0.8500, 0.3250, 0.0980]; % Red
c_dae  = [0.0000, 0.4470, 0.7410]; % Blue
c_back = [0.5000, 0.5000, 0.5000]; % Gray for Backbone
c_surf = [0.85 0.85 0.85];         % Manifold Surface Color

%% 1. Setup Model Parameters
om1 = 2.0;   % First natural frequency
om2 = 4.0;   % Second natural frequency (1:2 resonance)
om3 = 5.0;   
zeta = 0.01; 
m = 1; 

alpha = 0.8; % Saddle curvature x
beta  = 0.5; % Saddle curvature y 
gamma = 0.1; % Cubic stiffness
f1 = 1;      % Forcing amplitude

%% 2. Exact DAE Matrices and Tensors
n = 3;
M = m * eye(n);
C = diag([2*zeta*om1, 2*zeta*om2, 2*zeta*om3]);
K = diag([om1^2, om2^2, om3^2]);

% State vector: z = [x1, x2, x3, v1, v2, v3, lambda]
B = [C, M, zeros(3,1); M, zeros(3,3), zeros(3,1); zeros(1,7)];
A = [-K, zeros(3,3), [0;0;-1]; zeros(3,3), M, zeros(3,1); [0,0,1], zeros(1,3), 0];

% Quadratic tensors (Constraint equation + Lagrange multiplier feedback)
subs2 = [1 1 7; 1 7 1; 2 2 7; 2 7 2; 7 1 1; 7 2 2];            
vals2 = [alpha; alpha; -beta; -beta; -alpha; beta];
F2 = sptensor(subs2, vals2, [7, 7, 7]);

% Cubic tensors (Material nonlinearity)
subs3 = [1 1 1 1; 2 2 2 2; 3 3 3 3];
vals3 = [-gamma; -gamma; -gamma];
F3 = sptensor(subs3, vals3, [7, 7, 7, 7]);

Fnl = {F2, F3};
Fext = [f1; 0; 0; 0; 0; 0; 0];

%% 3. Dynamical System Setup
DS = DynamicalSystem();
set(DS, 'B', B, 'A', A, 'fnl', Fnl);
set(DS.Options, 'Emax', 6, 'Nmax', 10, 'notation', 'multiindex');

% Harmonic Forcing
epsilon = 0.02;
kappas = [-1; 1];
coeffs = [Fext, Fext] / 2;
DS.add_forcing(coeffs, kappas, epsilon);

%% 4. Linear Spectral Analysis & 4D Master Subspace
[V, D, W] = DS.linear_spectral_analysis();

S = SSM(DS);
set(S.Options, 'reltol', 0.5, 'notation', 'multiindex');
% 4D Subspace explicitly chosen to capture 1:2 internal resonance
resonant_modes = [1, 2, 3, 4]; 
S.choose_E(resonant_modes);

%% 5. Autonomous SSM Computation & Transient Validation
order_sim = 9; 
disp(['--> Computing 4D SSM Whisker O(', num2str(order_sim), ')...']);
t_ssm_start = tic;
[W0, R0] = S.compute_whisker(order_sim);
t_ssm = toc(t_ssm_start);

% INITIAL CONDITION
tf = 60; nsteps = 3000;
q0 = [0.15*exp(1i*0.5); 0.15*exp(-1i*0.5); 0.01*exp(1i*0.1); 0.01*exp(-1i*0.1)]; 
z0_SSM = reduced_to_full_traj(0, q0, W0);

% SSM Simulation
t_sim_start = tic;
traj = transient_traj_on_auto_ssm(DS, resonant_modes, W0, R0, tf, nsteps, 1:7, [], q0);
t_sim_ssm = toc(t_sim_start);

% Direct ODE15s Integration using Index-1 Baumgarte Stabilization
t_dae_start = tic;
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
[tInt, zInt] = ode15s(@(t,x) saddle_ode_core(t, x, 0, 2*pi, alpha, beta, gamma, om1, om2, om3, zeta, m), ...
                      [0 tf], z0_SSM(1:6), options);
t_dae = toc(t_dae_start);

% Extract Lagrange Multiplier (lambda) from DAE
lambda_Int = zeros(length(tInt), 1);
for i = 1:length(tInt)
    [~, lam] = saddle_ode_with_lambda(tInt(i), zInt(i,:)', 0, 2*pi, alpha, beta, gamma, om1, om2, om3, zeta, m);
    lambda_Int(i) = lam;
end

%% 6. Generate Outputs & Tables
x1_ssm_interp  = interp1(traj.time, traj.phy(:,1), tInt, 'spline');
x2_ssm_interp  = interp1(traj.time, traj.phy(:,2), tInt, 'spline');
x3_ssm_interp  = interp1(traj.time, traj.phy(:,3), tInt, 'spline');
lam_ssm_interp = interp1(traj.time, traj.phy(:,7), tInt, 'spline');

err_x1_max = max(abs(zInt(:,1) - x1_ssm_interp)); err_x1_rms = rms(zInt(:,1) - x1_ssm_interp);
err_x2_max = max(abs(zInt(:,2) - x2_ssm_interp)); err_x2_rms = rms(zInt(:,2) - x2_ssm_interp);
err_x3_max = max(abs(zInt(:,3) - x3_ssm_interp)); err_x3_rms = rms(zInt(:,3) - x3_ssm_interp);
err_lm_max = max(abs(lambda_Int  - lam_ssm_interp)); err_lm_rms = rms(lambda_Int - lam_ssm_interp);

fprintf('\n=======================================================================\n');
fprintf('        Q1 JOURNAL TABLE: TRANSIENT TRAJECTORY ERROR ANALYSIS          \n');
fprintf('=======================================================================\n');
fprintf('| State Variable | Max Absolute Error | RMS Error      | Relative L2  |\n');
fprintf('|----------------|--------------------|----------------|--------------|\n');
fprintf('| x_1 (Mode 1)   | %8.4e         | %8.4e     | %8.4e   |\n', err_x1_max, err_x1_rms, norm(zInt(:,1) - x1_ssm_interp)/norm(zInt(:,1)));
fprintf('| x_2 (Mode 2)   | %8.4e         | %8.4e     | %8.4e   |\n', err_x2_max, err_x2_rms, norm(zInt(:,2) - x2_ssm_interp)/norm(zInt(:,2)));
fprintf('| x_3 (Const.)   | %8.4e         | %8.4e     | %8.4e   |\n', err_x3_max, err_x3_rms, norm(zInt(:,3) - x3_ssm_interp)/norm(zInt(:,3)));
fprintf('| Lambda (Force) | %8.4e         | %8.4e     | %8.4e   |\n', err_lm_max, err_lm_rms, norm(lambda_Int - lam_ssm_interp)/norm(lambda_Int));
fprintf('=======================================================================\n\n');

%% 7. Visualizations: Manifolds & Transients 

% --- FIG 1: 3D Configuration Manifold ---
max_x1 = max(abs(zInt(:,1))) * 1.2;
max_x2 = max(abs(zInt(:,2))) * 1.2;
x1_span = linspace(-max_x1, max_x1, 60);
x2_span = linspace(-max_x2, max_x2, 60);
[X1_grid, X2_grid] = meshgrid(x1_span, x2_span);
X3_grid = alpha * X1_grid.^2 - beta * X2_grid.^2; 

fig1 = figure('Name', 'Configuration Manifold', 'Position', [100, 100, 800, 600], 'Color', 'w');
surf(X1_grid, X2_grid, X3_grid, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', c_surf); hold on;
camlight left; lighting gouraud; 
plot3(traj.phy(:,1), traj.phy(:,2), traj.phy(:,3), '-', 'Color', c_o11, 'LineWidth', 2.5);
plot3(zInt(:,1), zInt(:,2), zInt(:,3), '--', 'Color', 'k', 'LineWidth', 1.5);
xlabel('Displacement $x_1$', 'FontSize', 15); ylabel('Displacement $x_2$', 'FontSize', 15); zlabel('Displacement $x_3$', 'FontSize', 15);
title(['Exact Invariance on the Configuration Manifold $\mathcal{M}$ ($\mathcal{O}(' num2str(order_sim) ')$)'], 'FontSize', 16);
leg1 = legend('Constraint: $x_3 = \alpha x_1^2 - \beta x_2^2$', '4D SSM Trajectory', 'Full DAE Integration', 'Location', 'best');
set(leg1, 'EdgeColor', 'none', 'Color', 'none');
grid on; view([-35, 45]); set(gca, 'Box', 'on'); pbaspect([1 1 0.8]);

% --- FIG 2: 2D Transient & Energy Transfer ---
fig2 = figure('Name', 'Time Histories', 'Position', [150, 150, 850, 450], 'Color', 'w');
plot(traj.time, traj.phy(:,1), '-', 'Color', c_o11, 'LineWidth', 2); hold on;
plot(traj.time, traj.phy(:,2), '-', 'Color', c_o3, 'LineWidth', 2);
plot(tInt, zInt(:,1), 'k--', 'LineWidth', 1.2);
plot(tInt, zInt(:,2), 'b--', 'LineWidth', 1.2);
xlabel('Time $t$ [s]', 'FontSize', 15); ylabel('Modal Displacements', 'FontSize', 15);
title('1:2 Internal Resonance: Energy Beating Phenomenon', 'FontSize', 16);
leg2 = legend('$x_1$ (SSM)', '$x_2$ (SSM)', 'DAE Exact $x_1$', 'DAE Exact $x_2$', 'Location', 'northeast', 'NumColumns', 2);
set(leg2, 'EdgeColor', 'none', 'Color', 'none');
grid on; set(gca, 'Box', 'on'); xlim([0 tf]);


%% 8. FORCED RESPONSE CURVES (FRC) & BACKBONE TRACING FOR O(3, 5, 7, 11)
disp('--> Extracting Forced Response Curves (SSM) for Orders 3, 5, 7, 11...');

outdof = [1, 2]; 
freqRange = [1.4, 2.9]; 
eps_val = 0.02;

DS.add_forcing(coeffs, kappas, eps_val);
S = SSM(DS); S.choose_E(resonant_modes);
set(S.FRCOptions, 'coordinates', 'cartesian', 'initialSolver', 'fsolve', 'outdof', outdof);
set(S.contOptions, 'PtMX', 200, 'h_max', 0.05); % Increased PtMX for continuity
z0_guess = 1e-4 * ones(4, 1); 
set(S.FRCOptions, 'method', 'continuation ep', 'z0', z0_guess);

orders_list = [3, 5, 7, 11];
FRC_all = cell(length(orders_list), 1);
BB_all  = cell(length(orders_list), 1);

S2 = SSM(DS); S2.choose_E([1, 2]); % 2D Subspace for Backbone
set(S2.FRCOptions, 'outDOF', 1);

for i = 1:length(orders_list)
    ord = orders_list(i);
    disp(['Extracting SSM FRC and Backbone O(', num2str(ord), ')...']);
    try
        FRC_all{i} = S.extract_FRC('freq', freqRange, ord);
        BB_all{i}  = S2.extract_backbone([1, 2], freqRange, ord, 0.4);
    catch ME
        warning(['Failed at Order ', num2str(ord)]);
        disp(ME.message);
    end
end
%%
%%
%%

S.extract_FRC('freq', freqRange, orders_list)

S2.extract_backbone([1, 2], freqRange, orders_list, 0.4)
%%
%%
%%
%% 9. Full Order Model (FOM) Continuation via COCO
disp('--> Starting COCO Continuation for FOM validation...');
try
    omega_start = 1.4; T = 2*pi/omega_start;
    psp = [eps_val; omega_start];

    [~, x0_trans] = ode15s(@(t,x) saddle_ode_core(t, x, eps_val, omega_start, alpha, beta, gamma, om1, om2, om3, zeta, m), [0 200*T], zeros(6,1), options);
    [ttor, xtor]  = ode15s(@(t,x) saddle_ode_core(t, x, eps_val, omega_start, alpha, beta, gamma, om1, om2, om3, zeta, m), [0 T], x0_trans(end,:)', options);

    prob = coco_prob();
    prob = coco_set(prob, 'ode', 'autonomous', false);
    prob = coco_set(prob, 'cont', 'NAdapt', 2, 'h_max', 0.1, 'PtMX', 400); % Large PtMX for full [1.4, 2.9] sweep
    funcs = {@(t,x,p) saddle_ode_coco(t, x, p, alpha, beta, gamma, om1, om2, om3, zeta, m)};
    coll_args = [funcs, {ttor, xtor, {'eps' 'Omega'}, psp}];
    prob = ode_isol2po(prob, '', coll_args{:});

    [data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
    maps = data.coll_seg.maps;
    prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, [], 'zero', 'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(2))]);

    prob = coco_add_func(prob, 'amp1', @amplitude, struct('dof',1,'zdim',6), 'regular', 'x1', 'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);
    prob = coco_add_func(prob, 'amp2', @amplitude, struct('dof',2,'zdim',6), 'regular', 'x2', 'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);

    cont_args = {1, {'Omega' 'po.period' 'x1' 'x2' 'eps'}, freqRange};
    bd_coco = coco(prob, 'fom_freq_resp', [], cont_args{:});

    om_coco = coco_bd_col(bd_coco, 'Omega');
    x1_coco = coco_bd_col(bd_coco, 'x1');
    x2_coco = coco_bd_col(bd_coco, 'x2');
    eig_coco = coco_bd_col(bd_coco, 'eigs');
    stab_coco = all(abs(eig_coco) < 1 + 1e-6, 1);
    coco_success = true;
catch ME
    warning('COCO Continuation failed.');
    coco_success = false;
end

%% 10. PLOT: 2D & 3D FRCs (SSM vs COCO) 
% Map orders to colors
c_ord = {c_o3, c_o5, c_o7, c_o11};

% --- FIG 3: 2D FRC (Amplitude vs Frequency) & Backbone ---
fig3 = figure('Name', '2D FRC & Backbone', 'Position', [250, 250, 900, 650], 'Color', 'w'); 
ax2D = axes('Parent', fig3); hold(ax2D, 'on'); grid(ax2D, 'on');

h_dummy_leg_2D = []; legend_strs_2D = {};

for i = 1:length(orders_list)
    % Backbone
    if ~isempty(BB_all{i})
        plot(ax2D, BB_all{i}.Omega, BB_all{i}.Aout, '-.', 'Color', c_back, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    
    % FRC
    if ~isempty(FRC_all{i})
        FRC_O = FRC_all{i};
        if isfield(FRC_O(1), 'Aout_frc'), af='Aout_frc'; else, af='Aout'; end
        if isfield(FRC_O(1), 'om'), fq='om'; else, fq='Omega'; end
        
        for b = 1:length(FRC_O)
            if isempty(FRC_O(b).(af)), continue; end % Safeguard against empty branches
            
            om_vec = FRC_O(b).(fq); 
            a1_vec = FRC_O(b).(af)(:,1); 
            if isfield(FRC_O(b), 'stability'), st_vec = logical(FRC_O(b).stability);
            elseif isfield(FRC_O(b), 'st'), st_vec = logical(FRC_O(b).st); else, st_vec = true(size(a1_vec)); end
            
            a1_s = a1_vec; a1_s(~st_vec) = NaN;
            a1_u = a1_vec; a1_u(st_vec) = NaN;
            
            plot(ax2D, om_vec, a1_s, '-', 'Color', c_ord{i}, 'LineWidth', 2.0, 'HandleVisibility', 'off');
            plot(ax2D, om_vec, a1_u, '--', 'Color', c_ord{i}, 'LineWidth', 2.0, 'HandleVisibility', 'off');
        end
    end
    
    idx = length(h_dummy_leg_2D) + 1;
    h_dummy_leg_2D(idx) = plot(ax2D, NaN, NaN, '-', 'Color', c_ord{i}, 'LineWidth', 2.0);
    legend_strs_2D{idx} = ['SSM $\mathcal{O}(', num2str(orders_list(i)), ')$'];
end

% Add COCO Validation
if coco_success
    om_c = om_coco(1:4:end); x1_c = x1_coco(1:4:end); st_c = stab_coco(1:4:end); % Less dense for clean plot
    plot(ax2D, om_c(st_c), x1_c(st_c), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax2D, om_c(~st_c), x1_c(~st_c), 'kx', 'MarkerSize', 6, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    idx = length(h_dummy_leg_2D) + 1;
    h_dummy_leg_2D(idx) = plot(ax2D, NaN, NaN, 'ko', 'MarkerFaceColor', 'w'); legend_strs_2D{idx} = 'FOM (COCO) Stable';
    idx = length(h_dummy_leg_2D) + 1;
    h_dummy_leg_2D(idx) = plot(ax2D, NaN, NaN, 'kx', 'LineWidth', 1.5); legend_strs_2D{idx} = 'FOM (COCO) Unstable';
end

xlabel(ax2D, 'Forcing Frequency $\Omega$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel(ax2D, 'Amplitude $\|x_1\|_{\infty}$', 'Interpreter', 'latex', 'FontSize', 16);
title(ax2D, 'Forced Response \& Backbone: SSM vs Full System', 'Interpreter', 'latex', 'FontSize', 16);
legend(ax2D, h_dummy_leg_2D, legend_strs_2D, 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 12);
set(ax2D, 'Box', 'on', 'GridAlpha', 0.2);

% --- FIG 4: 3D FRC (Internal Resonance) ---
fig4 = figure('Name', '3D FRC (Internal Resonance)', 'Position', [300, 300, 900, 700], 'Color', 'w'); 
ax3D = axes('Parent', fig4); hold(ax3D, 'on'); grid(ax3D, 'on');

h_dummy_leg_3D = []; legend_strs_3D = {};

for i = 1:length(orders_list)
    if ~isempty(FRC_all{i})
        FRC_O = FRC_all{i};
        if isfield(FRC_O(1), 'Aout_frc'), af='Aout_frc'; else, af='Aout'; end
        if isfield(FRC_O(1), 'om'), fq='om'; else, fq='Omega'; end
        
        for b = 1:length(FRC_O)
            if isempty(FRC_O(b).(af)), continue; end 
            
            om_vec = FRC_O(b).(fq); 
            a1_vec = FRC_O(b).(af)(:,1); 
            
            if size(FRC_O(b).(af), 2) >= 2
                a2_vec = FRC_O(b).(af)(:,2); 
            else
                a2_vec = zeros(size(a1_vec)); 
            end
            
            if isfield(FRC_O(b), 'stability'), st_vec = logical(FRC_O(b).stability);
            elseif isfield(FRC_O(b), 'st'), st_vec = logical(FRC_O(b).st); else, st_vec = true(size(a1_vec)); end
            
            a1_s = a1_vec; a1_s(~st_vec) = NaN; a2_s = a2_vec; a2_s(~st_vec) = NaN;
            a1_u = a1_vec; a1_u(st_vec) = NaN;  a2_u = a2_vec; a2_u(st_vec) = NaN;
            
            plot3(ax3D, om_vec, a1_s, a2_s, '-', 'Color', c_ord{i}, 'LineWidth', 2.0, 'HandleVisibility', 'off');
            plot3(ax3D, om_vec, a1_u, a2_u, '--', 'Color', c_ord{i}, 'LineWidth', 2.0, 'HandleVisibility', 'off');
        end
    end
    
    idx = length(h_dummy_leg_3D) + 1;
    h_dummy_leg_3D(idx) = plot3(ax3D, NaN, NaN, NaN, '-', 'Color', c_ord{i}, 'LineWidth', 2.0);
    legend_strs_3D{idx} = ['SSM $\mathcal{O}(', num2str(orders_list(i)), ')$'];
end

if coco_success
    x2_c = x2_coco(1:4:end); % Less dense
    plot3(ax3D, om_c(st_c), x1_c(st_c), x2_c(st_c), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot3(ax3D, om_c(~st_c), x1_c(~st_c), x2_c(~st_c), 'kx', 'MarkerSize', 6, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    idx = length(h_dummy_leg_3D) + 1;
    h_dummy_leg_3D(idx) = plot3(ax3D, NaN, NaN, NaN, 'ko', 'MarkerFaceColor', 'w'); legend_strs_3D{idx} = 'FOM (COCO) Stable';
    idx = length(h_dummy_leg_3D) + 1;
    h_dummy_leg_3D(idx) = plot3(ax3D, NaN, NaN, NaN, 'kx', 'LineWidth', 1.5); legend_strs_3D{idx} = 'FOM (COCO) Unstable';
end

legend(ax3D, h_dummy_leg_3D, legend_strs_3D, 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 12);
xlabel(ax3D, 'Forcing Frequency $\Omega$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel(ax3D, 'Amplitude $\|x_1\|_{\infty}$', 'Interpreter', 'latex', 'FontSize', 15);
zlabel(ax3D, 'Amplitude $\|x_2\|_{\infty}$', 'Interpreter', 'latex', 'FontSize', 15);
title(ax3D, '3D FRC: Modal Interaction via 1:2 Internal Resonance', 'Interpreter', 'latex', 'FontSize', 16);
view(ax3D, [-45, 35]); set(ax3D, 'Box', 'on', 'GridAlpha', 0.2);
%%
disp('--> Starting COCO Continuation for FOM validation...');
try
    omega_start = 1.4; T = 2*pi/omega_start;
    psp = [epsilon; omega_start];
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

    [~, x0_trans] = ode15s(@(t,x) saddle_ode_core(t, x, epsilon, omega_start, alpha, beta, gamma, om1, om2, om3, zeta, m), [0 150*T], zeros(6,1), options);
    [ttor, xtor]  = ode15s(@(t,x) saddle_ode_core(t, x, epsilon, omega_start, alpha, beta, gamma, om1, om2, om3, zeta, m), [0 T], x0_trans(end,:)', options);

    prob = coco_prob();
    prob = coco_set(prob, 'ode', 'autonomous', false);
    prob = coco_set(prob, 'cont', 'NAdapt', 2, 'h_max', 0.1, 'PtMX', 500); 
    funcs = {@(t,x,p) saddle_ode_coco(t, x, p, alpha, beta, gamma, om1, om2, om3, zeta, m)};
    coll_args = [funcs, {ttor, xtor, {'eps' 'Omega'}, psp}];
    prob = ode_isol2po(prob, '', coll_args{:});

    [data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
    maps = data.coll_seg.maps;
    prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, [], 'zero', 'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(2))]);
    prob = coco_add_func(prob, 'amp1', @amplitude, struct('dof',1,'zdim',6), 'regular', 'x1', 'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);
    prob = coco_add_func(prob, 'amp2', @amplitude, struct('dof',2,'zdim',6), 'regular', 'x2', 'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);

    cont_args = {1, {'Omega' 'po.period' 'x1' 'x2' 'eps'}, freqRange};
    bd_coco = coco(prob, 'fom_freq_resp', [], cont_args{:});

    om_coco = coco_bd_col(bd_coco, 'Omega');
    x1_coco = coco_bd_col(bd_coco, 'x1');
    x2_coco = coco_bd_col(bd_coco, 'x2');
    eig_coco = coco_bd_col(bd_coco, 'eigs');
    stab_coco = all(abs(eig_coco) < 1 + 1e-5, 1);
    coco_success = true;
catch ME
    warning('COCO Continuation failed.');
    coco_success = false;
end
%%
figure
% set up FRC options
set(S.FRCOptions, 'nCycle',500, 'initialSolver', 'fsolve');
set(S.contOptions, 'PtMX', 300, 'h_max', 0.1);
set(S.FRCOptions, 'omegaSampStyle', 'cocoBD');
set(S.FRCOptions, 'outdof',1,'method','continuation ep','p0',[]);
S.extract_FRC('freq',freqRange,[3,5,7, 11]);

if coco_success
    hold on; 
    skip = 2; 
    om_c = om_coco(1:skip:end); x1_c = x1_coco(1:skip:end); st_c = stab_coco(1:skip:end);
    
    plot(om_c(st_c), x1_c(st_c), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.0, 'DisplayName', 'FOM (COCO) Stable');
    plot(om_c(~st_c), x1_c(~st_c), 'kx', 'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'FOM (COCO) Unstable');
    
    legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12);
    
    xlabel('Forcing Frequency $\Omega$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('Amplitude $\|x_1\|_{\infty}$', 'Interpreter', 'latex', 'FontSize', 16);
    title('FRC with COCO', 'Interpreter', 'latex', 'FontSize', 16);
    grid on; set(gca, 'Box', 'on');
end


%% ========================================================================
%% HELPER FUNCTIONS
%% ========================================================================

function dy = saddle_ode_coco(t, x, p, alpha, beta, gamma, om1, om2, om3, zeta, m)
    epf = p(1,:); om  = p(2,:);
    [dy, ~] = saddle_ode_core(t, x, epf, om, alpha, beta, gamma, om1, om2, om3, zeta, m);
end

function [dy, lambda] = saddle_ode_core(t, x, epf, om, alpha, beta, gamma, om1, om2, om3, zeta, m)
    x1 = x(1,:); x2 = x(2,:); x3 = x(3,:); 
    v1 = x(4,:); v2 = x(5,:); v3 = x(6,:);
    
    F1 = epf.*cos(om.*t) - 2*zeta*om1.*v1 - (om1^2).*x1 - gamma.*x1.^3;
    F2 = -2*zeta*om2.*v2 - (om2^2).*x2 - gamma.*x2.^3;
    F3 = -2*zeta*om3.*v3 - (om3^2).*x3 - gamma.*x3.^3;
    
    g = x3 - alpha.*x1.^2 + beta.*x2.^2;
    gdot = v3 - 2*alpha.*x1.*v1 + 2*beta.*x2.*v2;
    c_stab = 100.*g + 20.*gdot; 
    
    num = F3 - 2*alpha*m.*v1.^2 - 2*alpha.*x1.*F1 + 2*beta*m.*v2.^2 + 2*beta.*x2.*F2 + m.*c_stab;
    den = 1 + 4*alpha^2.*x1.^2 + 4*beta^2.*x2.^2;
    lambda = num ./ den;
    
    dy = zeros(6, size(x, 2));
    dy(1,:) = v1; dy(2,:) = v2; dy(3,:) = v3;
    dy(4,:) = (F1 + 2*alpha.*x1.*lambda) ./ m;
    dy(5,:) = (F2 - 2*beta.*x2.*lambda) ./ m;
    dy(6,:) = (F3 - lambda) ./ m;
end

function [dy, lambda] = saddle_ode_with_lambda(t, x, epf, om, alpha, beta, gamma, om1, om2, om3, zeta, m)
    x1 = x(1); x2 = x(2); x3 = x(3); v1 = x(4); v2 = x(5); v3 = x(6);
    F1 = epf*cos(om*t) - 2*zeta*om1*v1 - (om1^2)*x1 - gamma*x1^3;
    F2 = -2*zeta*om2*v2 - (om2^2)*x2 - gamma*x2^3;
    F3 = -2*zeta*om3*v3 - (om3^2)*x3 - gamma*x3^3;
    g = x3 - alpha*x1^2 + beta*x2^2;
    gdot = v3 - 2*alpha*x1*v1 + 2*beta*x2*v2;
    c_stab = 100*g + 20*gdot; 
    num = F3 - 2*alpha*m*v1^2 - 2*alpha*x1*F1 + 2*beta*m*v2^2 + 2*beta*x2*F2 + m*c_stab;
    den = 1 + 4*alpha^2*x1^2 + 4*beta^2*x2^2;
    lambda = num / den;
    dy = zeros(6,1);
    dy(1:3) = [v1; v2; v3];
    dy(4) = (F1 + 2*alpha*x1*lambda) / m;
    dy(5) = (F2 - 2*beta*x2*lambda) / m;
    dy(6) = (F3 - lambda) / m;
end

function [data, y] = OmegaT(prob, data, u) %#ok<INUSL>
    y = u(1)*u(2) - 2*pi; 
end

function [data, J] = OmegaT_du(prob, data, u) %#ok<INUSL>
    J = [u(2), u(1)];
end

function [data, y] = amplitude(prob, data, u) %#ok<INUSL>
    xbp = reshape(u, data.zdim, []);
    y = max(abs(xbp(data.dof, :)), [], 2);
end

function [prob, status, xtr] = amplitude_remesh(prob, data, chart, old_u, old_V) %#ok<INUSD>
    [colldata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
    maps = colldata.coll_seg.maps;
    xtr    = [];
    prob   = coco_change_func(prob, data, 'uidx', uidx(maps.xbp_idx));
    status = 'success';
end