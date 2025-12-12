%% Clear from previous run
clear;       % Clears all variables from the workspace (the kernel)
clc;         % Clears the Command Window
close all;   % Closes all open figure windows


%% Channel 1
data = readtable('Ch1.txt');
t = data{:,1};   
ch1 = data{:,2}; 

% --- 1. Clean NaN values ---
valid_idx = ~isnan(ch1);
t = t(valid_idx);
ch1 = ch1(valid_idx);

% --- 2. Gain correction ---
% 1. Convert the time '1:04:06' to seconds
gain_correction_time_str = '1:04:06';
parts = strsplit(gain_correction_time_str, ':');
H = str2double(parts{1});
M = str2double(parts{2});
S = str2double(parts{3});
gain_correction_time_s = H*3600 + M*60 + S; % Total seconds: 3846 s
% 2. Create the index for the gain correction period
gain_idx = t <= gain_correction_time_s;
% 3. Apply the gain correction factor (0.1) to the current values
ch1(gain_idx) = ch1(gain_idx) * 0.1;

% --- 3. Baseline correction ---
% Define Baseline Window (0:58:00 to 1:00:49 in seconds)
time_baseline_start_s = 3480; 
time_baseline_end_s = 3649; 
% Create index for the baseline period using the full time vector (t)
baseline_idx = (t >= time_baseline_start_s) & (t <= time_baseline_end_s);
% Calculate the Mean Baseline Current from the gain-corrected data
ch1_avg_baseline = mean(ch1(baseline_idx));


% --- 4. Filter data for time >= 3000 s---
time_start = 3000;
plot_idx = t >= time_start; % Creates a logical index for time >= 3000 s
% Apply the index to both time (t) and current (ch1)
t_plot = t(plot_idx);
ch1_plot = ch1(plot_idx)- ch1_avg_baseline; %Subtract baseline


% --- 5. Define Concentrations and Addition Timestamps ---
conc = [0.2 0.4 0.6 0.8 1.0 2.0 4.0 6.0 8.0 10.0 15.0 20.0]; % nM (12 values)
% Updated time stamps provided in HH:MM:SS format (12 values)
time_stamps = {'1:00:49', '1:02:48', '1:04:54', '1:07:44', '1:09:43', '1:12:28', ...
               '1:15:13', '1:16:30', '1:17:17', '1:18:08', '1:20:44', '1:23:17'};      
%Convert to seconds
injection_times_s = zeros(1, length(time_stamps)); 
for k = 1:length(time_stamps)
    time_str = time_stamps{k};
    parts = strsplit(time_str, ':');
    
    if length(parts) == 3
        H = str2double(parts{1}); % Hours
        M = str2double(parts{2}); % Minutes
        S = str2double(parts{3}); % Seconds
        
        % Convert to total seconds
        injection_times_s(k) = H*3600 + M*60 + S;
    end
end

% --- 6. Plot the filtered data ---
figure;
% Added DisplayName for the main plot
plot(t_plot, ch1_plot, 'b-', 'DisplayName', 'Current Response'); 
hold on;
% Add Vertical Lines. The loop index k directly corresponds to the concentration index
for k = 1:length(injection_times_s)
    current_time = injection_times_s(k);
    current_conc = conc(k); % Use k as the index
    
    % Use xline to draw a vertical line
    xline(current_time, 'r--', ...
          sprintf('%.1f nM', current_conc), ...
          'LineWidth', 1.5, ...
          'LabelOrientation', 'aligned'); 
end
hold off;

xlabel('Time (s)')
ylabel('Current (nA)')
title(['Sensor 1 calibration (Time > ' num2str(time_start) ' s)'])
grid on


%% Channel 2
data = readtable('Ch2.txt');
t = data{:,1};   
ch2 = data{:,2}; 

% --- 1. Clean NaN values ---
valid_idx = ~isnan(ch2);
t = t(valid_idx);
ch2 = ch2(valid_idx);

% --- 2. Baseline correction ---
% Define Baseline Window (0:58:00 to 1:00:49 in seconds)
time_baseline_start_s = 3480; 
time_baseline_end_s = 3649; 
% Create index for the baseline period using the full time vector (t)
baseline_idx = (t >= time_baseline_start_s) & (t <= time_baseline_end_s);
% Calculate the Mean Baseline Current from the gain-corrected data
ch2_avg_baseline = mean(ch2(baseline_idx));

% --- 3. Filter data for time >= 3000 s ---
time_start = 3000;
plot_idx = t >= time_start; % Creates a logical index for time >= 3000 s
% Apply the index to both time (t) and current (ch2)
t_plot = t(plot_idx);
ch2_plot = ch2(plot_idx)-ch2_avg_baseline;


% --- 4. Define Concentrations and Addition Timestamps ---
conc = [0.2 0.4 0.6 0.8 1.0 2.0 4.0 6.0 8.0 10.0 15.0 20.0]; % nM (12 values)
% Updated time stamps provided in HH:MM:SS format (12 values)
time_stamps = {'1:00:49', '1:02:48', '1:04:54', '1:07:44', '1:09:43', '1:12:28', ...
               '1:15:13', '1:16:30', '1:17:17', '1:18:08', '1:20:44', '1:23:17'};     
%Convert to seconds
injection_times_s = zeros(1, length(time_stamps)); 
for k = 1:length(time_stamps)
    time_str = time_stamps{k};
    parts = strsplit(time_str, ':');
    
    if length(parts) == 3
        H = str2double(parts{1}); % Hours
        M = str2double(parts{2}); % Minutes
        S = str2double(parts{3}); % Seconds
        
        % Convert to total seconds
        injection_times_s(k) = H*3600 + M*60 + S;
    end
end

% --- 5. Plot the filtered data ---
figure;
% Added DisplayName for the main plot
plot(t_plot, ch2_plot, 'b-', 'DisplayName', 'Current Response'); 
hold on;
% Add Vertical Lines. The loop index k directly corresponds to the concentration index
for k = 1:length(injection_times_s)
    current_time = injection_times_s(k);
    current_conc = conc(k); % Use k as the index
    
    % Use xline to draw a vertical line
    xline(current_time, 'r--', ...
          sprintf('%.1f nM', current_conc), ...
          'LineWidth', 1.5, ...
          'LabelOrientation', 'aligned'); 
end
hold off;

xlabel('Time (s)')
ylabel('Current (nA)')
title(['Sensor 2 calibration (Time > ' num2str(time_start) ' s)'])
grid on



%% Channel 3
data = readtable('Ch3.txt');
t = data{:,1};   
ch3 = data{:,2}; 

% --- 1. Clean NaN values ---
valid_idx = ~isnan(ch3);
t = t(valid_idx);
ch3 = ch3(valid_idx);

% --- 2. Gain correction ---
% 1. Convert the time '1:04:06' to seconds
gain_correction_time_str = '1:13:43';
parts = strsplit(gain_correction_time_str, ':');
H = str2double(parts{1});
M = str2double(parts{2});
S = str2double(parts{3});
gain_correction_time_s = H*3600 + M*60 + S; % Total seconds: 3846 s
% 2. Create the index for the gain correction period
gain_idx = t <= gain_correction_time_s;
% 3. Apply the gain correction factor (0.1) to the current values
ch3(gain_idx) = ch3(gain_idx) * 0.1;

% --- 3. Baseline correction ---
% Define Baseline Window (0:58:00 to 1:00:49 in seconds)
time_baseline_start_s = 3480; 
time_baseline_end_s = 3649; 
% Create index for the baseline period using the full time vector (t)
baseline_idx = (t >= time_baseline_start_s) & (t <= time_baseline_end_s);
% Calculate the Mean Baseline Current from the gain-corrected data
ch3_avg_baseline = mean(ch3(baseline_idx));

% --- 4. Filter data for time >= 3000 s---
time_start = 3000;
plot_idx = t >= time_start; % Creates a logical index for time >= 3000 
% Apply the index to both time (t) and current (ch1)
t_plot = t(plot_idx);
ch3_plot = ch3(plot_idx)-ch3_avg_baseline;


% --- 5. Define Concentrations and Addition Timestamps ---
conc = [0.2 0.4 0.6 0.8 1.0 2.0 4.0 6.0 8.0 10.0 15.0 20.0]; % nM (12 values)
% Updated time stamps provided in HH:MM:SS format (12 values)
time_stamps = {'1:00:49', '1:02:48', '1:04:54', '1:07:44', '1:09:43', '1:12:28', ...
               '1:15:13', '1:16:30', '1:17:17', '1:18:08', '1:20:44', '1:23:17'};    
%Convert to seconds
injection_times_s = zeros(1, length(time_stamps)); 
for k = 1:length(time_stamps)
    time_str = time_stamps{k};
    parts = strsplit(time_str, ':');
    
    if length(parts) == 3
        H = str2double(parts{1}); % Hours
        M = str2double(parts{2}); % Minutes
        S = str2double(parts{3}); % Seconds
        
        % Convert to total seconds
        injection_times_s(k) = H*3600 + M*60 + S;
    end
end

% --- 6. Plot the filtered data ---
figure;
% Added DisplayName for the main plot
plot(t_plot, ch3_plot, 'b-', 'DisplayName', 'Current Response'); 
hold on;

% Add Vertical Lines. The loop index k directly corresponds to the concentration index
for k = 1:length(injection_times_s)
    current_time = injection_times_s(k);
    current_conc = conc(k); % Use k as the index
    
    % Use xline to draw a vertical line
    xline(current_time, 'r--', ...
          sprintf('%.1f nM', current_conc), ...
          'LineWidth', 1.5, ...
          'LabelOrientation', 'aligned'); 
end
hold off;

xlabel('Time (s)')
ylabel('Current (nA)')
title(['Sensor 3 calibration (Time > ' num2str(time_start) ' s)'])
grid on

