% Plot the STFT spectrum and the LP spectrum of order p=4,10,18 of an frame
% on the same graph. Define the p(opt) by that.
% The window size is 30ms(240samples), no overlapping.
% We choose one voiced frame and one unvoiced frame to run analysis.
close all
clear
clc
FRAME = 1; % '0' for unvoiced
if FRAME == 0
    num = 67; % unvoiced frame #
else 
    num = 9; % voiced frame #
end
[x,fs] = audioread('s5.wav');
s = buffer(x,100);
data = s(num,:);
subplot(2,1,1)
plot(data);
xlabel('Sample')
ylabel('Amplitude')
title(sprintf('%d th frame',num));
w = hamming(240);
data_w = data.*w';
X = fft(data_w);
n = 0:(length(X)-1)/2;
f = fs*n/length(X);
X = X(1:length(X)/2);
X = abs(X);
STFT = 10*log10(X/max(X));
subplot(2,1,2)
plot(f,STFT);
xlabel('Frequency(Hz)')
ylabel('Log Magnitude(dB)')
hold on

% LP Spectrum
p = [4;10;18];
for j=1:3
    m = p(j);
    L = length(data_w);
    r = [];
    for i=0:m
        r = [r; sum(data_w(1:L-i).* data_w(1+i:L))];
    end
    RR = toeplitz(r(1:m));
    a = inv(RR)*r(2:m+1);
    al = [1;-a];
    G1 = sqrt(sum(al.*r));
    data_pre = filter([0 -al(2:end)'],1,data_w);
    err = data_w - data_pre;
    data_syn_err = filter(1,[1 al(2:end)'],err);
    [h,w] = freqz(1,[1 al(2:m+1)']);
    HH(:,j) = 10*log10(abs(h)/abs(max(h)));
end
ww = w*fs/pi/2;
plot(ww,HH(:,1),ww,HH(:,2),ww,HH(:,3));
legend('STFT','p=4','p=10','p=18');


    


