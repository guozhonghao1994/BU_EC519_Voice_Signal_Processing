% Part2 - Pitch Detection
close all
clear 
clc
% Import the audio file and choose mode
filename = {'train1.wav' 'train2.wav' 'train3.wav' 'train4.wav' 'train5.wav' 'train6.wav' 'train7.wav'};
MODE = 1; % 0 for TRAIN, 1 for TEST

% WHITE GAUSSIAN NOISE
NOISE = 1; % 0 for original sound, 1 for noised sound
SNR = 0; % set the SNR here

if MODE == 0
    % TRAINING MODE
    for i=1:length(filename)
        [x,Fs] = wavread(cell2mat(filename(i)));
        % Filtering
        [B,A] = ellip(5,2,50,500*2/Fs);
        x = filter(B,A,x);
        pitchwatch(x,Fs)
        key = pitchwatch(x,Fs);
    end
    %figure
    %freqz(B,A)
else
    % TESTING MODE
    [x,Fs] = wavread('test.wav');
    if NOISE == 1
        x = awgn(x,SNR,'measured','db');
    end
    [B,A] = ellip(5,2,50,500*2/Fs);
    x = filter(B,A,x);
    pitchwatch(x,Fs)
    key = pitchwatch(x,Fs);
    %figure
    %freqz(B,A)
end

