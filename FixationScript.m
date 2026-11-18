%% Load Dataset
clc
clearvars
load(fullfile('Data', 'mp4_fixation.mat')) 
%% Defining Variables
fs = 1000; % Sampling frequency (Hz)
num_participants = 19; % Number of particpants in dataset

% Frequency bands (Hz)
theta = [4 8];
alpha = [8 13];
beta = [13 30];
gamma = [30 80];
%% Preprocessing
% Designing a custom bandpass filter to increase efficiency and designing a
% bandstop filter to remove 60Hz spiking (notch filtering)

% Creating empty cells for overall dataset as well as for each frequency
% band
processed_data = cell(1,num_participants);
processed_theta = cell(1,num_participants);
processed_alpha = cell(1,num_participants);
processed_beta = cell(1,num_participants);
processed_gamma = cell(1, num_participants);

% Filter order
n = 100; 

% Designing bandpass filters for total dataset, and designing banstop
% filter to remove unwanted noise at 60Hz
bp_filter = designfilt('bandpassfir','FilterOrder', n, 'CutoffFrequency1', 0.5,'CutoffFrequency2',120,'SampleRate', fs);
bs_filter = designfilt('bandstopfir','FilterOrder', n, 'CutoffFrequency1', 58,'CutoffFrequency2', 62,'SampleRate', fs);

% Designing bandpass filters to isolate each frequency phase
theta_filter = designfilt('bandpassfir','FilterOrder', n,'CutoffFrequency1', theta(1),'CutoffFrequency2',theta(2),'SampleRate', fs);
alpha_filter = designfilt('bandpassfir','FilterOrder', n,'CutoffFrequency1', alpha(1),'CutoffFrequency2',alpha(2),'SampleRate', fs);
beta_filter = designfilt('bandpassfir','FilterOrder', n,'CutoffFrequency1', beta(1),'CutoffFrequency2',beta(2),'SampleRate', fs);
gamma_filter = designfilt('bandpassfir','FilterOrder', n,'CutoffFrequency1', gamma(1),'CutoffFrequency2',gamma(2),'SampleRate', fs);

% Implementing filters
for p = 1:num_participants
    % Filtering total dataset
    bandpassed_data = filtfilt(bp_filter, data{p});
    filtered_data = filtfilt(bs_filter, bandpassed_data);
    % Assigning filtered data to cell
    processed_data{p} = filtered_data;
    % Filtering to isolate each frequency phase
    processed_theta{p} = filtfilt(theta_filter,filtered_data);
    processed_alpha{p} = filtfilt(alpha_filter,filtered_data);
    processed_beta{p} = filtfilt(beta_filter,filtered_data);
    processed_gamma{p} = filtfilt(gamma_filter,filtered_data);
    % Sanity check to show that data has been processed
    disp(['Processed Participant ', num2str(p)])
end
%% Matching number of samples for participant data
% Each participant will have a sample size equal to the participant with
% the smallest sample size

% Variable is infinity so the first participants' sample size is always smaller
min_samples = inf;

for participants = 1:19 % Total number of participants
    [rows, ~] = size(processed_data{participants}); 
    % Focuses on total number of samples
    if rows < min_samples
        min_samples = rows;
        % Replacing the current value of min_samples with current
        % particpants samples, if they are smaller than the previous entry
    end
end

% Only taking into account the minimum amount of samples per participant
for p = 1:num_participants
    processed_data{p} = processed_data{p}(1:min_samples,:);
    processed_theta{p} = processed_theta{p}(1:min_samples,:);
    processed_alpha{p} = processed_alpha{p}(1:min_samples,:);
    processed_beta{p} = processed_beta{p}(1:min_samples,:);
    processed_gamma{p} = processed_gamma{p}(1:min_samples,:);
end
% This leaves total electrodes and only alters sample size
%% Fourier Transform of each electrode before computing average

% Creating empty cells to help with storage
power_all_participants = cell(1,num_participants);
power_all_participants_theta = cell(1,num_participants);
power_all_participants_alpha = cell(1,num_participants);
power_all_participants_beta = cell(1,num_participants);
power_all_participants_gamma = cell(1,num_participants);


% Computing Fast Fourier Transform (FFT) on dataset to get power spectrum
for p = 1:num_participants

    % Assigning new dataset to variable 'y' to avoid mixup of variables,
    % and then removing DC offset 
    Y = processed_data{p};
    Y = Y - mean(Y,1);

    % Computing FFT
    FFT_result = fft(Y);
    [total_samples, ~] = size(Y);
    % Rounds down so k is not a fraction
    % Keeping only the first half of information gained of the FFT because
    % FFT is symmetric - there is no new information gained
    k = floor(total_samples/2);

    % Finding both the magnitude and power spectrum from the FFT
    magnitude_all = (2/total_samples) * abs(FFT_result(1:k,:));
    power_all_participants{p} = magnitude_all.^2;
end
% Comments are sparse for the rest of the section: it is a repetetive
% process and the previous comments in lines 91 - 109 can also be used,
% just focusing on different frequency bands

% It was more effecient to run this section in this manner as opposed to
% creating a custom function

for p = 1:num_participants

    T = processed_theta{p};
    T = T - mean(T,1);

    FFT_theta = fft(T);
    [total_samples, ~] = size(T);

    magnitude_all_theta = (2/total_samples) * abs(FFT_theta(1:k,:));
    power_all_participants_theta{p} = magnitude_all_theta.^2;

end

for p = 1:num_participants

    A = processed_alpha{p};
    A = A - mean(A,1);

    FFT_alpha = fft(A);
    [total_samples, ~] = size(A);

    magnitude_all_alpha = (2/total_samples) * abs(FFT_alpha(1:k,:));
    power_all_participants_alpha{p} = magnitude_all_alpha.^2;

end

for p = 1:num_participants

    B = processed_beta{p};
    B = B - mean(B,1);

    FFT_beta = fft(B);
    [total_samples, ~] = size(B);

    magnitude_all_beta = (2/total_samples) * abs(FFT_beta(1:k,:));
    power_all_participants_beta{p} = magnitude_all_beta.^2;

end

for p = 1:num_participants

    G = processed_gamma{p};
    G = G - mean(G,1);

    FFT_gamma = fft(G);
    [total_samples, ~] = size(G);

    magnitude_all_gamma = (2/total_samples) * abs(FFT_gamma(1:k,:));
    power_all_participants_gamma{p} = magnitude_all_gamma.^2;

    
end
% Creating a frequency vector
f = (0:k-1)*(fs/total_samples);
num_freqs = length(f);
%% Calculate Power Spectrum Average across Participants
% Allows for cleaner visualisation of overall neural activity across
% participants

% Creating emptyy cells
average_power = cell(1,num_participants);
average_power_theta = cell(1,num_participants);
average_power_alpha = cell(1,num_participants);
average_power_beta = cell(1,num_participants);
average_power_gamma = cell(1,num_participants);

% Calculating the average total power spectrum and average frequency phase power spectrum for each participant
for p = 1:num_participants
    average_power{p} = mean(power_all_participants{p},2);
    average_power_theta{p} = mean(power_all_participants_theta{p},2);
    average_power_alpha{p} = mean(power_all_participants_alpha{p},2);
    average_power_beta{p} = mean(power_all_participants_beta{p},2);
    average_power_gamma{p} = mean(power_all_participants_gamma{p},2);
end

% Using the average power spectrums to find averages across participants
% Total dataset
PSD_matrix = cell2mat(average_power);
all_participants_mean = mean(PSD_matrix,2);

% Theta phase
theta_matrix = cell2mat(average_power_theta);
all_participants_mean_theta = mean(theta_matrix,2);

% Alpha phase
alpha_matrix = cell2mat(average_power_alpha);
all_participants_mean_alpha = mean(alpha_matrix,2);

% Beta phase
beta_matrix = cell2mat(average_power_beta);
all_participants_mean_beta = mean(beta_matrix,2);

% Gamma phase
gamma_matrix = cell2mat(average_power_gamma);
all_participants_mean_gamma = mean(gamma_matrix,2);
%% Visualisation
% Figure 1 will plot the average power spectrum across all participants,
% for both the total dataset and all frequency phases

figure(1)
% Tiled layout to improve figure appearance
t = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
nexttile([1 2]);
% Power spectrum across participants
plot(f,10*log10(all_participants_mean), 'w--', 'Linewidth', 1), title('Average Power Spectrum Across Participants')
xlabel('Frequency (Hz)'), ylabel('Power (dB)')
xlim([1 120]), ylim([0 75]) % Limits represent the passband cutoffs
grid on

% Theta oscillations
nexttile
plot(f, 10*log10(all_participants_mean_theta), 'g--', 'Linewidth', 1), title('Theta Phase')
xlabel('Frequency (Hz)'), ylabel('Power (dB)')
xlim([1 10]), ylim([0 100]) % Limits represent the passband cutoffs
grid on
% Alpha oscillations
nexttile
plot(f, 10*log10(all_participants_mean_alpha), 'r--', 'Linewidth', 1), title('Alpha Phase')
xlabel('Frequency (Hz)'), ylabel('Power (dB)')
xlim([5 20]), ylim([0 100]) % Limits represent the passband cutoffs
grid on
% Beta oscillations
nexttile
plot(f, 10*log10(all_participants_mean_beta), 'b--', 'Linewidth', 1), title('Beta Phase')
xlabel('Frequency (Hz)'), ylabel('Power (dB)')
xlim([10 35]), ylim([0 100]) % Limits represent the passband cutoffs
grid on
% Gamma oscillations
nexttile
plot(f, 10*log10(all_participants_mean_gamma),'y--', 'Linewidth', 1), title('Gamma Phase')
xlabel('Frequency (Hz)'), ylabel('Power (dB)')
xlim([30 85]), ylim([0 100]) % Limits represent the passband cutoffs
grid on
%% Power Spectrum of each Brain Region
% Grouping electrodes based off of their coordinates

% Defining matrices for each brain region
occipital   = zeros(num_freqs,num_participants);
frontal = zeros(num_freqs,num_participants);
temporal  = zeros(num_freqs,num_participants);
parietal   = zeros(num_freqs,num_participants);

for p = 1:19

    % Temporary variables for easier formatting
    loc_participant = locs{p};
    pap = power_all_participants{p};
    total_electrodes = size(pap,2);

    % Defining x, y and z as loc uses xyz coordinates
    x = loc_participant(:,1);
    y = loc_participant(:,2);
    z = loc_participant(:,3);

    % Normalize coordinates for more consistent mapping of electrodes
    x = x / max(abs(x));
    z = z / max(abs(z));

    % Assigning values to brain regions when they correspond with the
    % coordinates of said regions - based on Talairach system

    % Occipital region
    % Adds all data points to index when the z coordinate is less than -0.4
    occ_idx   = find(z < -0.4);
    % Removes invalid indicies
    occ_idx   = occ_idx(occ_idx <= total_electrodes);
    % Only calculate average if electrodes exist in this region
    if ~isempty(occ_idx)
        occipital(:,p) = mean(pap(:,occ_idx),2);
    end

    % Frontal Lobe
    % Adds all data points to index when the z coordinate is above than 0.4
    front_idx = find(z >  0.4);
    % Removes invalid indicies
    front_idx = front_idx(front_idx <= total_electrodes);
    % Only calculate average if electrodes exist in this region
    if ~isempty(front_idx)
        frontal(:,p) = mean(pap(:,front_idx),2);
    end

    % Temporal lobe
    % Adds all data points to index when the x coordinate is above 0.5 and
    % the z coordinate is below 0.3
    temp_idx  = find(abs(x) > 0.5 & abs(z) < 0.3);
    % Removes invalid indicies
    temp_idx  = temp_idx(temp_idx <= total_electrodes);
    % Only calculate average if electrodes exist in this region
    if ~isempty(temp_idx)
        temporal(:,p) = mean(pap(:,temp_idx),2);
    end

    % Parietal lobe
    % Adds all data points to index when the x coordinate is less or equal to 0.5 and
    % the z coordinate is less or equal to 0.4
    par_idx   = find(abs(z) <= 0.4 & abs(x) <= 0.5);
    % Removes invalid indicies
    par_idx   = par_idx(par_idx <= total_electrodes);
    % Only calculate average if electrodes exist in this region
    if ~isempty(par_idx)
        parietal(:,p) = mean(pap(:,par_idx),2);
    end
end

% Calculating average power spectrum in each brain region
average_occipital_electrodes   = mean(occipital,2);
average_frontal_electrodes = mean(frontal,2);
average_temporal_electrodes  = mean(temporal,2);
average_parietal_electrodes   = mean(parietal,2);

%% Visualisation 
% Figure 2 is a 2x2 subplot of the average power spectrum in all four regions
figure(2)

% Occipital Lobe
subplot(2,2,1)
plot(f,10*log10(average_occipital_electrodes), 'm--', 'LineWidth',2)
title('Occipital Lobe'), xlabel('Frequency (Hz)'), ylabel('Power (dB)')
xlim([1 120]), ylim([-25 100]) % Limits represent the upper passband cutoffs
grid on

% Temporal Lobe
subplot(2,2,2)
plot(f,10*log10(average_temporal_electrodes), 'r--', 'LineWidth',2)
title('Temporal Lobe'), xlabel('Frequency (Hz)'), ylabel('Power (dB)')
xlim([1 120]), ylim([-25 100]) % Limits represent the upper passband cutoffs
grid on

% Frontal Lobe
subplot(2,2,3)
plot(f,10*log10(average_frontal_electrodes),'g--', 'LineWidth',2)
title('Frontal Lobe'), xlabel('Frequency (Hz)'), ylabel('Power (dB)')
xlim([1 120]), ylim([-25 100]) % Limits represent the upper passband cutoffs
grid on

% Parietal Lobe
subplot(2,2,4)
plot(f,10*log10(average_parietal_electrodes),'y--','LineWidth',2)
title('Parietal Lobe'), xlabel('Frequency (Hz)'), ylabel('Power (dB)')
xlim([1 120]), ylim([-25 100]) % Limits represent the passband cutoffs
grid on 
% Title for entire subplot
sgtitle('Average Electrode Activity Across Brain Regions')