%% Glucose sensor Channel 2
%% Clear from previous run
clear;       % Clears all variables from the workspace (the kernel)
clc;         % Clears the Command Window
close all;   % Closes all open figure windows

%% Step Fitting 
data = readtable('Ch2.txt');
t = data{:,1};   
ch = data{:,2}; 

% --- 1. Clean NaN values ---
valid_idx = ~isnan(ch);
t = t(valid_idx);
ch = ch(valid_idx);

% --- 2. Baseline correction ---
% Define Baseline Window that you want to average
time_baseline_start_s = time2sec('0:58:00');
time_baseline_end_s   = time2sec('1:00:49');
% Create index for the baseline period
baseline_idx = (t >= time_baseline_start_s) & (t <= time_baseline_end_s);
% Calculate the Mean Baseline Current
ch_avg_baseline = mean(ch(baseline_idx));
ch_plot_baseline_corrected = ch - ch_avg_baseline;%Subtract baseline
new_baseline_avg = mean(ch_plot_baseline_corrected(baseline_idx));

% --- 3. Define Concentrations and Injection Timestamps ---
conc = [0.2 0.4 0.6 0.8 1.0 2.0 4.0 6.0 8.0 10.0 15.0 20.0];
time_stamps = {'1:00:49', '1:02:48', '1:04:54', '1:07:44', '1:09:43', '1:12:28', ...
               '1:15:13', '1:16:30', '1:17:17', '1:18:08', '1:20:44', '1:23:17'}; 
%Convert to seconds
injection_times_s = zeros(1, length(time_stamps)); 
for k = 1:length(time_stamps)
    injection_times_s(k)=time2sec(time_stamps{k});
end

% --- 4. Calculate Plateau Averages and Construct Step Signal ---
% Pad injection times with the end of the data to define the last plateau's end time
all_step_times = [injection_times_s, t(end)]; 
ch_step_fit = zeros(size(t)); 
plateau_values = zeros(1, length(conc));
time_to_start_meas_s = 40; % Start averaging 40 seconds after injection to allow settling
time_buffer_s = 5;         % Stop averaging 5 seconds before the next injection

% Loop through all injection times
for k = 1:length(conc) 
    
    % Define the Measurement Window 
    step_start_time = all_step_times(k);
    step_end_time = all_step_times(k+1); 
    
    % Define the stable averaging region
    avg_start_time = step_start_time + time_to_start_meas_s;
    avg_end_time = step_end_time - time_buffer_s; 
    
    % Create the index for the averaging region within t
    avg_idx = (t >= avg_start_time) & (t < avg_end_time);     
    data_for_avg = ch_plot_baseline_corrected(avg_idx); 

    % Remove outliers defined as 3sd away from mean
    % Step A: Calculate the median and standard deviation of the window
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
    
    % Find the index of the full step plateau
    plateau_idx = (t >= step_start_time) & (t <= step_end_time);
    
    % Assign the calculated average to the full plateau region in the new vector
    ch_step_fit(plateau_idx) = step_avg;
end

% --- 5. Plot the Original Data and the Step-Fit Signal ---
figure;
%Filter time for above 3000s
time_start = 3000;
plot_idx = (t >= time_start); 
t_plot = t(plot_idx);
ch_plot_baseline_corrected_cut=ch_plot_baseline_corrected(plot_idx);
ch_step_fit_cut=ch_step_fit(plot_idx);

% Plot the original data
plot(t_plot, ch_plot_baseline_corrected_cut, 'b-', 'DisplayName', 'Original Data (Noisy Steps)', 'LineWidth', 0.5); 
hold on;grid on;
hold off;
xlabel('Time (s)');
ylabel('Current (nA)');
title(['Sensor 2 Data']);

%Make another plot with the data together
figure;
%Filter time for above 3000s
time_start = 3000;
plot_idx = (t >= time_start); 
t_plot = t(plot_idx);
ch_plot_baseline_corrected_cut=ch_plot_baseline_corrected(plot_idx);
ch_step_fit_cut=ch_step_fit(plot_idx);

% Plot the original data
plot(t_plot, ch_plot_baseline_corrected_cut, 'b-', 'DisplayName', 'Original Data (Noisy Steps)', 'LineWidth', 0.5); 
hold on;

% Plot the step fit
plot(t_plot, ch_step_fit_cut, 'r-', 'DisplayName', 'Idealized Flat-Line Steps', 'LineWidth', 1.5); 

% Add Vertical Lines for Injection Times
for k = 1:length(injection_times_s)
    current_time = injection_times_s(k);
    current_conc = conc(k); 
    
    xline(current_time, 'k--', sprintf('%.1f nM', current_conc),'LineWidth', 1.0, ...
          'LabelOrientation', 'aligned', 'FontWeight', 'bold'); 
end

hold off;
xlabel('Time (s)');
ylabel('Current (nA)');
title(['Sensor 2 Data and Step Fit']);
grid on

%% Callibration curve
% Adding the baseline values to the data
plateau_values=[new_baseline_avg, plateau_values];
conc = [0 0.2 0.4 0.6 0.8 1.0 2.0 4.0 6.0 8.0 10.0 15.0 20.0];

% Choose which values to be plotted
idx1 = conc <= 1; %This is chosen as it gives a straight line
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
noise_data = ch(baseline_idx);
baseline_noise = std(noise_data);

% Calculate 3 times the standard deviation (3*SD)
three_sd = 3 * baseline_noise;

% Calculate LOD concentration (x-intercept at 3*SD)
% x_LOD = (3*sigma - b) / m
x_lod_1 = (three_sd - y_intercept_1) / slope_1;

% Plotting the Raw Data
plot(x, y, 'bo', 'LineWidth', 2, 'DisplayName', 'Raw Data Points'); hold on;

% Plotting the Line of Best Fit
eq1 = sprintf('y = %.4f x + %.4f', slope_1, y_intercept_1);
x_intercept = -y_intercept_1 / slope_1;
xfit_extended = linspace(x_intercept, max(x), 100);
yfit_extended = polyval(p, xfit_extended);
plot(xfit_extended, yfit_extended, 'b-', 'LineWidth', 2, 'DisplayName', ['LR: ',eq1]);

% 3*SD Horizontal Line
plot(xlim, [three_sd three_sd], 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('3 \\times SD (%.4f nA)', three_sd));

% LOD Label at the intersection point
plot(x_lod_1, three_sd, 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('LOD: %.4f mM', x_lod_1));

% Formatting
ylim([-0.1, 0.8]);
xlim([-0.1, 1.1]);
xlabel('Concentration of Glucose(mM)');
ylabel('Current (nA)');
title('Sensor 2 Callibration Curve');
grid on;
legend('Location', 'best');

%% Michaelis Menten curve
% The equationMichaelis-Menten Equation is: V = (Vmax * [S]) / (Km + [S])
michaelis_menten = @(params, conc) (params(1) .* conc) ./ (params(2) + conc);
% params(1) is Vmax, params(2) is Km

% Fit the Data to the equation
% Initial parameter guesses are needed: Vmax, Km
Vmax_guess = max(plateau_values);
Km_guess = median(conc); 
p0 = [Vmax_guess, Km_guess];

% Perform the fit
[params, resnorm, residual, exitflag, output, lambda, jacobian] = ...
    lsqcurvefit(michaelis_menten, p0, conc, plateau_values);
Vmax = params(1);
Km = params(2);

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

% Plot the horizontal line for Vmax
Vmax_legend_label = sprintf('V_{max} = %.2f nA', Vmax);
plot(xlim, [Vmax Vmax], 'k--', 'LineWidth', 1, 'DisplayName', Vmax_legend_label);

% Plot the horizontal line for 1/2 Vmax
light_blue = [0.6 0.8 1.0]; % RGB for Light Blue
half_Vmax = Vmax / 2; % The y-coordinate for Km
plot(xlim, [half_Vmax half_Vmax], '--','Color', light_blue,'LineWidth', 1, 'DisplayName', '1/2 V_{max}');

% Plot the vertical line for Km
light_magenta = [1.0 0.6 1.0]; % RGB for Light Magenta
Km_legend_label = sprintf('K_{M} = %.2f mM', Km);
plot([Km Km], [0 half_Vmax], '--','Color', light_magenta, 'LineWidth', 1, 'DisplayName', Km_legend_label);

% Mark the intersection point (Km, 1/2 Vmax)
plot(Km, half_Vmax, 'm*', 'MarkerSize', 5,'Color', [0.8 0.3 0.8], 'LineWidth', 1.5, 'HandleVisibility', 'off'); 

% Formatting 
xlabel('Concentration of Glucose (mM)');
ylabel('Current Generated (nA)');
title('Sensor 2 Michaelis-Menten Fit');
ylim([0 2.2]);
xlim([0 20.5]);
grid on;
legend('Location', 'best');
hold off;

% Calculate sensitivity based on curve
% Sensitivity = dV/d[S] at [S]. For M-M: (Vmax * Km) / (Km + [S])^2
michaelis_menten_sensitivity = @(Vmax, Km, S) (Vmax * Km) / (Km + S)^2;
% Calculate sensitivity at 1mM
S = 1.0;          % mM 
sensitivity = michaelis_menten_sensitivity(Vmax, Km, S);

%% t90 Response time calculation
% Filter plot range
plot_idx = (t > 3550) & (t < 3750); 
t_plot = t(plot_idx);
ch_plot_baseline_corrected_cut = ch_plot_baseline_corrected(plot_idx);
ch_step_fit_cut = ch_step_fit(plot_idx);

figure; 
% Plot the original data
plot(t_plot, ch_plot_baseline_corrected_cut, 'b-', 'DisplayName', 'Original Data', 'LineWidth', 0.5); 
hold on;

% Plot the step fit
plot(t_plot, ch_step_fit_cut, 'r-', 'DisplayName', 'Step Fit', 'LineWidth', 2); 

% Plot the bottom threashold at 3*SD+baseline mean
bottom_threshold_val = new_baseline_avg + (3 * baseline_noise);
yline(bottom_threshold_val, 'g-', 'DisplayName', '0mM mean + 3*SD', 'LineWidth', 1.5); 

% Plot the top threashold for a 90% change 
top_threshold_val = 0.9 * plateau_values(2);
yline(top_threshold_val, 'm-', 'DisplayName', '90% of 0.2mM mean', 'LineWidth', 1.5);

% Plot the adition time
xline(injection_times_s(1), 'y--', 'DisplayName', 'Injection time', 'LineWidth', 1.5); 

% Find the first index where the data crosses the top threshold
idx_crossing = find(ch_plot_baseline_corrected_cut >= top_threshold_val, 1, 'first');
if ~isempty(idx_crossing)
    % Extract the specific time of crossing
    t_90 = t_plot(idx_crossing);
    % Plot the vertical line (xline)
    xline(t_90, 'k--', 'LineWidth', 1.5,'DisplayName', 'Intersection point');
end

% Calculate the relative response time (t90 - injection time)
response_time = t_90 - injection_times_s(1);
% We plot a line to show response time at a height slightly above the threshold
span_height = top_threshold_val + 0.05; 
plot([injection_times_s(1), t_90], [span_height, span_height], 'k|-', 'LineWidth', 2.5, 'HandleVisibility', 'off');
% Add a label
t_mid = (injection_times_s(1) + t_90) / 2;
text(t_mid, span_height + 0.000001,sprintf('t_{90} = %.1f s', response_time), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontSize', 10, 'FontWeight', 'bold');

%Formating
xlabel('Time (s)');
ylabel('Current (nA)');
title('Sensor 2 Response Time');
grid on;
legend('Location', 'best');
hold off;

%% Function that converts string to seconds
function total_seconds = time2sec(time_str)
    % Converts 'H:M:S' string to total seconds
    parts = strsplit(time_str, ':');
    if length(parts) == 3
        H = str2double(parts{1}); 
        M = str2double(parts{2}); 
        S = str2double(parts{3});
        total_seconds = H*3600 + M*60 + S;
    else
        error('Time string must be in H:M:S format');
    end
end