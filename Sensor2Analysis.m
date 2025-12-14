%% Glucose sensor Channel 2
%% Clear from previous run
clear;       % Clears all variables from the workspace (the kernel)
clc;         % Clears the Command Window
close all;   % Closes all open figure windows

%% Step Fitting 
% --- Parameters for Step Fitting ---
data = readtable('Ch2.txt');
t = data{:,1};   
ch2 = data{:,2}; 

% --- 1. Clean NaN values ---
valid_idx = ~isnan(ch2);
t = t(valid_idx);
ch2 = ch2(valid_idx);

% --- 3. Baseline correction ---
% Define Baseline Window that you want to average
time_baseline_start_str = '0:58:00';
parts = strsplit(time_baseline_start_str, ':');
H = str2double(parts{1}); M = str2double(parts{2}); S = str2double(parts{3});
time_baseline_start_s = H*3600 + M*60 + S;  
time_baseline_end_str = '1:00:49';
parts = strsplit(time_baseline_end_str, ':');
H = str2double(parts{1}); M = str2double(parts{2}); S = str2double(parts{3});
time_baseline_end_s = H*3600 + M*60 + S;

% Create index for the baseline period using the full time vector (t)
baseline_idx = (t >= time_baseline_start_s) & (t <= time_baseline_end_s);
% Calculate the Mean Baseline Current from the gain-corrected data
ch2_avg_baseline = mean(ch2(baseline_idx));

ch2_plot_baseline_corrected = ch2 - ch2_avg_baseline;%Subtract baseline
new_baseline_avg = mean(ch2_plot_baseline_corrected(baseline_idx));

% --- 6. Define Concentrations and Injection Timestamps ---
conc = [0.2 0.4 0.6 0.8 1.0 2.0 4.0 6.0 8.0 10.0 15.0 20.0];
time_stamps = {'1:00:49', '1:02:48', '1:04:54', '1:07:44', '1:09:43', '1:12:28', ...
               '1:15:13', '1:16:30', '1:17:17', '1:18:08', '1:20:44', '1:23:17'}; 
%Convert to seconds
injection_times_s = zeros(1, length(time_stamps)); 
for k = 1:length(time_stamps)
    time_str = time_stamps{k};
    parts = strsplit(time_str, ':');
    if length(parts) == 3
        H = str2double(parts{1}); M = str2double(parts{2}); S = str2double(parts{3});
        injection_times_s(k) = H*3600 + M*60 + S;
    end
end

% --- 7. STEP FITTING: Calculate Averages and Construct Flat-Line Signal ---
% Pad injection times with the end of the data to define the last plateau's end time
all_step_times = [injection_times_s, t(end)]; 
ch2_step_fit = zeros(size(t)); % Initialize the new step-fit vector
plateau_values = zeros(1, length(conc));
% Start the loop from the first injection time (k=1) up to the last injection (length(conc))
time_to_start_meas_s = 40; % Start averaging 30 seconds after injection (to allow settling)
time_buffer_s = 5;         % Stop averaging 5 seconds before the next injection
for k = 1:length(conc) 
    
    % Define the Measurement Window (for the k-th concentration plateau)
    step_start_time = all_step_times(k);
    step_end_time = all_step_times(k+1); 
    
    % Define the stable averaging region
    avg_start_time = step_start_time + time_to_start_meas_s;
    avg_end_time = step_end_time - time_buffer_s; 
    
    % Create the index for the averaging region within t
    avg_idx = (t >= avg_start_time) & (t < avg_end_time);
    
    % --- Robust Outlier Removal for Averaging ---
    data_for_avg = ch2_plot_baseline_corrected(avg_idx); % Use the original, unsmoothed data
    % Step A: Calculate the median and standard deviation of the window
    % The median is a robust estimate of the true average.
    M = median(data_for_avg);
    S = std(data_for_avg);
    % Step B: Define the clipping threshold- 3sd in this case
    N_SIGMA = 3.0; 
    % Step C: Create a new index that excludes outliers
    % We only average points that are within M +/- N_SIGMA*S
    outlier_free_idx = (data_for_avg >= (M - N_SIGMA*S)) & ...
                       (data_for_avg <= (M + N_SIGMA*S));                       
    % Step D: Calculate the final average using only the outlier-free data
    step_avg = mean(data_for_avg(outlier_free_idx));

    % Store the calculated average plateau value
    plateau_values(k) = step_avg; 
    
    % Find the index of the full step plateau (from injection k to injection k+1)
    plateau_idx = (t >= step_start_time) & (t < step_end_time);
    
    % Assign the calculated average to the full plateau region in the new vector
    ch2_step_fit(plateau_idx) = step_avg;

end

% --- 8. Plot the Original Data and the Step-Fit Signal ---
figure;
%Filter time for above 3000s
time_start = 3000;
plot_idx = t >= time_start; 
t_plot = t(plot_idx);
ch2_plot_baseline_corrected=ch2_plot_baseline_corrected(plot_idx);
ch2_step_fit=ch2_step_fit(plot_idx);

% Plot the ORIGINAL (baseline-corrected) data
plot(t_plot, ch2_plot_baseline_corrected, 'b-', 'DisplayName', 'Original Data (Noisy Steps)', 'LineWidth', 0.5); 
hold on;

% Plot the IDEALIZED FLAT-LINE Step Data
plot(t_plot, ch2_step_fit, 'r-', 'DisplayName', 'Idealized Flat-Line Steps', 'LineWidth', 2.5); 

% Add Vertical Lines (Injection Times)
for k = 1:length(injection_times_s)
    current_time = injection_times_s(k);
    current_conc = conc(k); 
    
    xline(current_time, 'k--', ... 
          sprintf('%.1f nM', current_conc), ...
          'LineWidth', 1.0, ...
          'LabelOrientation', 'aligned', 'FontWeight', 'bold'); 
end
hold off;

xlabel('Time (s)');
ylabel('Current (nA)');
title(['Sensor 2 Data and Step Fit']);
grid on

%% Callibration curve
% Produce a List of all values to be plotted
plateau_values=[new_baseline_avg, plateau_values];
conc = [0 0.2 0.4 0.6 0.8 1.0 2.0 4.0 6.0 8.0 10.0 15.0 20.0];

% Choose which values to be plotted (before plateau)
idx1 = conc <= 1;
x = conc(idx1);
y = plateau_values(idx1);
figure;

% Fit linear model (p(1) = slope, p(2) = y-intercept)
p = polyfit(x, y, 1);
slope_1 = p(1);
y_intercept_1 = p(2);

% --- LOD CALCULATION ---
% Find the Standard Deviation of the Noise Before injections
% Use the same time window as the one used for baseline average calculation 
noise_data = ch2(baseline_idx);
sigma_1 = std(noise_data);
% Calculate 3 times the standard deviation (3*SD)
three_sigma_1 = 3 * sigma_1;
% Calculate LOD concentration (x-intercept at 3*SD)
% x_LOD = (3*sigma - b) / m
x_lod_1 = (three_sigma_1 - y_intercept_1) / slope_1;

% --- PLOTTING ---
plot(x, y, 'bo', 'LineWidth', 2, 'DisplayName', 'Raw Data Points'); hold on;

% Line of Best Fit
eq1 = sprintf('y = %.4f x + %.4f', slope_1, y_intercept_1);
x_intercept = -y_intercept_1 / slope_1;
xfit_extended = linspace(x_intercept, max(x), 100);
yfit_extended = polyval(p, xfit_extended);
plot(xfit_extended, yfit_extended, 'b-', 'LineWidth', 2, 'DisplayName', ['LR: ',eq1, sprintf(' Sensitivity:(%.4f nA/mM)', slope_1)]);

% 3*SD Horizontal Line
plot(xlim, [three_sigma_1 three_sigma_1], 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('3 \\times SD (%.4f nA)', three_sigma_1));

% LOD Label at the intersection point
plot(x_lod_1, three_sigma_1, 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('LOD: %.4f mM', x_lod_1));

% Formatting
ylim([-0.1, 0.8]);
xlim([-0.5, 1.5]);
xlabel('Concentration of Glucose(mM)');
ylabel('Current (nA)');
title('Sensor 2 Callibration Curve');
grid on;
legend('Location', 'best');

%% Michaelis Menten curve
% --- Define the Michaelis-Menten Equation ---
% The equation is: V = (Vmax * [S]) / (Km + [S])
michaelis_menten = @(params, conc) (params(1) .* conc) ./ (params(2) + conc);
% params(1) is Vmax, params(2) is Km

% --- Fit the Data ---
% Initial parameter guesses are needed: Vmax, Km
Vmax_guess = max(plateau_values);
Km_guess = median(conc); 
p0 = [Vmax_guess, Km_guess];

% Perform the fit
[params, resnorm, residual, exitflag, output, lambda, jacobian] = ...
    lsqcurvefit(michaelis_menten, p0, conc, plateau_values);

Vmax = params(1);
Km = params(2);

% --- Generate and Plot the Fitted Curve---
figure;
% Generate a range of concentrations for plotting the fitted curve
conc_range = linspace(min(conc), max(conc), 500);
fitted_values = michaelis_menten(params, conc_range);
% Plot the fitted curve 
plot(conc_range, fitted_values, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Data');
hold on;
% Plot the raw data
plot(conc, plateau_values, 'bo', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Data');
hold on;

% 1. Plot the horizontal line for Vmax
Vmax_legend_label = sprintf('V_{max} = %.2f nA', Vmax);
plot(xlim, [Vmax Vmax], 'k--', 'LineWidth', 1, 'DisplayName', Vmax_legend_label);

% 2. Plot the horizontal line for 1/2 Vmax
light_blue = [0.6 0.8 1.0]; % RGB for Light Blue
half_Vmax = Vmax / 2; % The y-coordinate for Km
plot(xlim, [half_Vmax half_Vmax], '--','Color', light_blue,'LineWidth', 1, 'DisplayName', '1/2 V_{max}');

% 3. Plot the vertical line for Km
light_magenta = [1.0 0.6 1.0]; % RGB for Light Magenta
Km_legend_label = sprintf('K_{M} = %.2f mM', Km);
plot([Km Km], [0 half_Vmax], '--','Color', light_magenta, 'LineWidth', 1, 'DisplayName', Km_legend_label);

% 4. Mark the intersection point (Km, 1/2 Vmax)
plot(Km, half_Vmax, 'm*', 'MarkerSize', 5,'Color', [0.8 0.3 0.8], 'LineWidth', 1.5, 'HandleVisibility', 'off'); 


% Formatting 
xlabel('Concentration of Glucose (mM)');
ylabel('Current Generated (nA)');
title('Sensor 2 Michaelis-Menten Curve');
grid on;
legend('Location', 'best');
hold off;

% --- Calculate Sensitivity ---
% Define the parameters (using fitted parameters for calculation)
Vmax_sens = Vmax; % nA
Km_sens = Km;     % mM
S = 1.0;          % mM 

% Sensitivity = dV/d[S] at [S]. For M-M: (Vmax * Km) / (Km + [S])^2
michaelis_menten_sensitivity = @(Vmax, Km, S) (Vmax * Km) / (Km + S)^2;

% Calculate sensitivity
sensitivity = michaelis_menten_sensitivity(Vmax_sens, Km_sens, S);
