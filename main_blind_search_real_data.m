% main_blind_search_real_data.m
% Conceptual implementation of the blind LEO beacon search algorithm
% designed to work with a user-provided IQ data recording.
% No toolboxes are used (only core MATLAB functions).

clear; clc; close all;

%% 1. USER-DEFINED PARAMETERS & DATA LOADING
% =========================================================================
% --- Provide your recording information here ---

% Path to your recorded data file (must be a .mat file)
filepath = 'path/to/your/signal.mat'; % <--- CHANGE THIS

% The name of the variable inside the .mat file that holds the IQ data
iq_variable_name = 'iq_data'; % <--- CHANGE THIS (e.g., 'rx_signal', 'data')

% The sampling frequency of your recording in Hz
fs = 2.4e6; % <--- CHANGE THIS (e.g., 2.4e6 for Iridium, 50e6 for OneWeb)

% The estimated period of the repetitive beacon signal in seconds
% This is a CRITICAL parameter.
T0 = 0.001; % <--- CHANGE THIS (e.g., 1ms for a 1kHz repetition rate)

% =========================================================================

% --- Load and Prepare Data ---
fprintf('Loading data from: %s\n', filepath);
if ~exist(filepath, 'file')
    error('Data file not found. Please check the `filepath`.');
end

try
    data_struct = load(filepath);
    if ~isfield(data_struct, iq_variable_name)
        error('The variable `%s` was not found in the .mat file.', iq_variable_name);
    end
    rx_signal = data_struct.(iq_variable_name);
catch ME
    error('Failed to load data. MATLAB error: %s', ME.message);
end

% --- Data Validation ---
if ~isvector(rx_signal) || ~isnumeric(rx_signal)
    error('Loaded data must be a numeric vector.');
end
if isreal(rx_signal)
    warning('The loaded signal is real. This algorithm is designed for complex IQ data. Results may be incorrect.');
end

% Ensure signal is a column vector
rx_signal = rx_signal(:);

fprintf('Data loaded successfully. Total samples: %d\n', length(rx_signal));

% --- Calculate Processing Parameters ---
L = round(fs * T0); % Length of the sequence in samples
num_blocks = floor(length(rx_signal) / L);

if num_blocks < 10 % Need a reasonable number of blocks to converge
    error('Not enough data for processing. Need at least 10 blocks of length L=%d samples.', L);
end

fprintf('Calculated block length L = %d samples.\n', L);
fprintf('Total number of blocks to process: %d\n\n', num_blocks);


%% 2. RECEIVER IMPLEMENTATION & PROCESSING
% --- Initialization
% Blind initialization: the first block is the first, noisy estimate
s_hat = rx_signal(1:L);
s_hat = s_hat / norm(s_hat); % Normalize initial estimate

% Kalman Filter State for Doppler Tracking: [freq_offset; freq_rate]
% Initialize with zero Doppler assumption
kf_state = [0; 0]; 
% Initialize with high uncertainty
P = diag([1e4^2, 1000^2]); % Covariance: [Hz^2, (Hz/s)^2]

% NCO State
nco_phase = 0;
nco_freq = 0;

% Storage for results
estimated_doppler = zeros(num_blocks, 1);
estimated_code_phase = zeros(num_blocks, 1);

fprintf('Starting blind search loop...\n');
% --- Processing Loop (process one block at a time)
for k = 1:num_blocks
    % Get the current block of the received signal
    current_block_idx = (k-1)*L + 1 : k*L;
    current_block = rx_signal(current_block_idx);

    % --- STEP I: Blind Doppler Tracking ---
    [freq_error, ~] = blind_doppler_discriminator(current_block, s_hat, fs);

    % --- STEP II: Kalman Filter Update ---
    [kf_state, P] = kalman_filter_update(kf_state, P, freq_error, T0);
    estimated_doppler(k) = kf_state(1); % The filtered Doppler estimate

    % Update the NCO frequency based on the KF output
    nco_freq = estimated_doppler(k);
    
    % --- STEP III: Doppler and Carrier Phase Wipe-off ---
    nco_time_vec = (0:L-1)' / fs;
    nco_signal = exp(-1i * (2 * pi * nco_freq * nco_time_vec + nco_phase));
    doppler_corrected_block = current_block .* nco_signal;
    
    % Update NCO phase for the next block (essential for phase continuity)
    nco_phase = mod(nco_phase + 2 * pi * nco_freq * T0, 2*pi);

    % --- STEP IV: Code Phase Tracking ---
    code_phase_offset = code_phase_tracker(doppler_corrected_block, s_hat);
    estimated_code_phase(k) = code_phase_offset;
    
    % Correct the code phase (circularly shift the block to align it)
    aligned_block = circshift(doppler_corrected_block, -code_phase_offset);
    
    % --- STEP V: Blind Beacon Estimation Update ---
    alpha = 0.05; % Update factor (like a learning rate). Smaller values mean slower, more stable updates.
    s_hat_new = aligned_block / norm(aligned_block); % Normalize current block
    s_hat = (1 - alpha) * s_hat + alpha * s_hat_new;
    s_hat = s_hat / norm(s_hat); % Re-normalize the beacon estimate

    if mod(k, 20) == 0 || k == 1
        fprintf('Block %d/%d: Est. Doppler = %.2f kHz, Est. Code Phase = %d samples\n', ...
            k, num_blocks, estimated_doppler(k)/1e3, code_phase_offset);
    end
end

fprintf('Processing complete.\n');

%% 3. PLOT RESULTS
time_axis = (1:num_blocks) * T0;

figure('Name', 'Blind Search Results', 'NumberTitle', 'off');

% Doppler Tracking Plot
subplot(3,1,1);
plot(time_axis, estimated_doppler/1e3, 'b', 'LineWidth', 1.5);
title('Estimated Doppler Frequency vs. Time');
xlabel('Time (s)');
ylabel('Frequency (kHz)');
grid on;
legend('Estimated Doppler');

% Code Phase Tracking Plot
subplot(3,1,2);
plot(time_axis, estimated_code_phase, 'r.');
title('Estimated Code Phase Offset vs. Time');
xlabel('Time (s)');
ylabel('Offset (samples)');
grid on;
ylim([-L/2, L/2]);
legend('Estimated Phase Offset');

% Final Estimated Beacon Plot
subplot(3,1,3);
plot(real(s_hat), 'b'); 
hold on; 
plot(imag(s_hat), 'r');
title('Final Estimated Beacon Signal (ŝ)');
xlabel('Sample Index');
ylabel('Amplitude');
legend('In-Phase (I)', 'Quadrature (Q)');
grid on;
xlim([0 L]);
