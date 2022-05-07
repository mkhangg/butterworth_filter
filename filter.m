clc;
clear;
close all;
warning('off','all')

%% 1. Audio Frequency Analysis
% 1a. Read audio file and get sampling frequency (Fs)
[data, Fs] = audioread('noisyaudio.wav');

% 1b. DFT of the retrived data
audio_fft = fft(data);

% 1c. Determine frequency axis values for plotting the DFT.
sz = size(data);
id = [1 : sz];
w = linspace(-Fs/2, Fs/2, sz(1));

figure;
width=1250;
height=400;
set(gcf,'position',[100,200,width,height])

% Plot noisy audio waveform
subplot(2, 3, 1);
plot(id, data)
grid on;
title('Noisy Audio Waveform');
xlabel('Time');
ylabel('Amplitude');

% 1d. Plot DFT magnitude vs frequency with 𝜔 = 0 at the center of the plot
subplot(2, 3, 2);
plot(w, fftshift(abs(audio_fft)));
grid on;
title('DFT of Noisy Audio');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% 1e. Plot DFT normalized log, where the max value of the DFT is 0 dB.
subplot(2, 3, 3);
data_db = mag2db(fftshift(abs(audio_fft)));
plot(w, data_db - max(data_db))
grid on;
title('Normalized Log of Noisy Audio');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

%% 2. Filter Design
% Obtaining the order and cutoff frequency of the Butterworth filter (See in report)
N = 17;
omega_c = 1987.5334;

%% 3. Filter Implementation
% 3a. Use butter() and filter() command to filter the noisy audio
[b,a] = butter(N, omega_c/(Fs/2));
filtered = filter(b, a, data);

% Plot the waveform of the filtered audio 
subplot(2, 3, 4);
plot(id, filtered);
grid on;
title('Filtered Audio Waveform');
xlabel('Time');
ylabel('Amplitude');

% Plot DFT magnitude of filtered audio vs frequency with 𝜔 = 0
filtered_fft = fft(filtered);
subplot(2, 3, 5);
plot(w, fftshift(abs(filtered_fft)));
grid on;
title('DFT of Filtered Audio');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% Plot DFT normalized log of filtered audio
subplot(2, 3, 6);
filtered_data_db = mag2db(fftshift(abs(filtered_fft)));
plot(w, filtered_data_db - max(abs(filtered_data_db)));
grid on;
title('Normalized Log of Filtered Audio');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

%% 2e. Plot the logarithmic gain of the frequency response
H = [];
for i = 1 : sz(1)
    H(i) = gain(i, omega_c, N);
end
figure;
plot(id, H);
grid on;
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');

%% Listen the filtered audio file and write to new .wav file
sound(filtered, Fs);
audiowrite('filteredaudio.wav', filtered, Fs);

%% Supplementary function
function H = gain(omega, omega_c, N)
    H = 20*log10(sqrt(1/(1+(omega/omega_c)^(2*N))));
end