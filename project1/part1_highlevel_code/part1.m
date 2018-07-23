close all
clear
clc
% Read the file and calculate ShortTimeEnergy,SpectralCentroid
filename = 'test.wav';
[x,fs] = wavread(filename);
NOISE = 0; % 0 for oroginal signal; 1 for noised signal
SNR = 20; % The SNR can be 0dB,10dB or 20dB
if NOISE == 1
    x = awgn(x,SNR,'measured','db');
end
win = 0.020;
step = 0.020;
E = ShortTimeEnergy(x,win*fs,step*fs);
E = 10*log(E)-max(10*log(E));
C = SpectralCentroid(x,win*fs,step*fs,fs);
w_length = floor(length(x)/length(C));
subplot(3,1,1)
plot(x);xlabel('sample');ylabel('amplitude')
xlim([1 length(x)])
title(filename)
subplot(3,1,2)
plot(E,'r');xlabel('frame');ylabel('Log Energy')
xlim([1 length(E)])
subplot(3,1,3)
plot(C,'g');xlabel('frame');ylabel('Zero Crossing Rate')
xlim([1 length(C)])
eavg = mean(E(1:2));
esig = std(E(1:2));
zavg = mean(C(1:2));
zsig = mean(C(1:2));
E_TH = max(-28,eavg+3*esig);
C_TH = max(0.4,zavg+3*zsig);
SILENCE = -80;
% STARTPOINT
for i=1:length(E)-1
    if E(i)<=SILENCE && E(i+1)>=SILENCE
        STARTPOINT = i*w_length;
        break
    end
end
% ENDPOINT
for i=length(E):-1:2
    if E(i)<=SILENCE && E(i-1)>=SILENCE
        ENDPOINT = i*w_length;
        break
    end
end
% B1
B1 = [];
for i=1:length(E)-1
    if E(i)<=E_TH && E(i+1)>=E_TH
        B1 = [B1 i];
    end
end
% E1
E1 = [];
for i=1:length(E)-1
    if E(i)>=E_TH && E(i+1)<=E_TH
        E1 = [E1 i];
    end
end
% plot voiced parts
figure
subplot(2,1,1)
plot(x,'b');xlabel('sample');ylabel('amplitude')
xlim([1 length(x)]);hold on
for j=1:length(B1)
    temp = x(B1(j)*w_length:E1(j)*w_length);
    plot(B1(j)*w_length:E1(j)*w_length,temp,'r');hold on
end
ylim = get(gca,'Ylim');
plot([STARTPOINT,STARTPOINT],ylim,'m--');hold on
plot([ENDPOINT,ENDPOINT],ylim,'m--');hold on
legend('original','voiced')
if NOISE == 1
    filename = [filename '    SNR=' num2str(SNR)];
    title(filename)
else
title(filename)
end

% B2
B2 = [];
for i=1:length(C)-1
    if C(i)<=C_TH && C(i+1)>=C_TH
        B2 = [B2 i];
    end
end
% E2
E2 = [];
for i=1:length(C)-1
    if C(i)>=C_TH && C(i+1)<=C_TH
        E2 = [E2 i];
    end
end
% plot unvoiced parts
subplot(2,1,2)
plot(x);xlabel('sample');ylabel('amplitude')
xlim([1 length(x)]);hold on
for j=1:length(B2)
    temp = x(B2(j)*w_length:E2(j)*w_length);
    plot(B2(j)*w_length:E2(j)*w_length,temp,'g');hold on
end
ylim = get(gca,'Ylim');
plot([STARTPOINT,STARTPOINT],ylim,'m--');hold on
plot([ENDPOINT,ENDPOINT],ylim,'m--');hold on
legend('original','un-voiced')





