% part1 - Quantization
clear all
close all
clc
%% Set Parameters
filename = 'H.2.wav'; %s5/H.1/H.2
info_rate = 2400; %2400/3600/4800
start_frame = 30; 
end_frame = 30; 
frame_length = 330; % For H.1 & H.2, change this value to 330 because fs'=11kHz
[x,fs] = audioread(filename);
x = x(1:(length(x)-mod(length(x),frame_length))); % In case the file can not be divided evenly
load 'h1.mat' %pp5.mat/h1.mat/h2.mat(got by Praat)
pitch = h1; % load pitch value

%% Encoder
switch info_rate
    case 2400
        g_bit = [5 5 5 5 4 4 4 4 3 2];
        ampl = 5;
    case 3600
        g_bit = [10 10 10 10 8 8 8 8 7 6];
        ampl = 7;
    case 4800
        g_bit = [14 14 13 13 11 11 11 11 10 10];
        ampl = 9;
end

%% Quantize pitch file
pitch_q(1:length(pitch)) = 0;
for i = 1:length(pitch)
    if pitch(i) == 0
        pitch_q(i) = 0;
    else
        pitch_q(i) = quantiz(pitch(i)-55,1:1:2^7);
    end
end
%% Quantize V-T parameters
L = frame_length;
R = frame_length;
n = round((length(x)/R));
if strcmp(filename,'s5.wav')
    n = n-10;
end
p = 10;
g_max = [5.5 4 2.5 1.5 1.5 2.5 1.5 1.5 1.5 1.5]';
al = [];
G = [];
g_uq = [];
sig_q = [];
for k=1:n
x_frame = x(1+(k-1)*R:((k-1)*R+L));
al_frameu = lpc(x_frame,p);
r = zeros(p+1,1);
for i=0:p
    r(i+1) = x_frame(1:L-i)'*x_frame(1+i:L);
end
G_frame=sqrt(al_frameu*r);
al_frameu = -al_frameu(2:p+1);
al_frame(1:p,1)=0;
for i=1:p
    al_frame(i) = al_frameu(i);
end
% alpha to PARCOR
alpha(1:p,1:p)=0; 
alpha(1:p,p)=al_frame; 
kk(1:p,1)=0; kk(p)=alpha(p,p);
for i=p:-1:2
    for j=1:i-1
            alpha(j,i-1)=(alpha(j,i)+kk(i)*alpha(i-j,i))/(1-kk(i).^2);
    end
    kk(i-1)=alpha(i-1,i-1);
end
g(1:p,1)=0;
% PARCOR to log area
for i=1:p
    g(i,1)=log((1-kk(i,1))/(1+kk(i,1)));
end
sig_frame(1:p,1)=0;
for i=1:p
    if g(i)<0
        sig_frame(i)=0;
    else
        sig_frame(i)=1;
    end
end
G = [G G_frame];
al = [al al_frameu'];
sig_q = [sig_q sig_frame];
% amplitude parameter might be quantized with a ¦Ì-law quantizer
%using 4-5 bits per frame
g_uq = [g_uq g];
end
g_q(p,1:n) = 0;
k_g(1:p,1:n) = 0;
for i=1:n
    for k=1:p
        k_g(k,i)=compand(abs(g_uq(k,i)),100,g_max(k,1),'mu/compressor');
        g_q(k,i)=quantiz(k_g(k,i),0:g_max(k,1)/2^g_bit(1,k):g_max(k,1));
    end
end
G_mu=compand(G,100,1.5,'mu/compressor');
G_q=quantiz(G_mu,0:1/2^ampl:1.1);

%% Receiver(get synthesized speech file)
%[al_de,G_de,pitch_de] = vt_unquantize(pitch_q,g_q,G_q,sig_q,g_bit,ampl);
[col,row]=size(g_q);
    g_bint(1:col,1)=0;
    for i=1:col
        g_bint(i,1)=g_max(i,1)/2^g_bit(1,i);
    end
    G_bint=1/2^ampl;
    k_exp(1:col,1:row)=0; G_de(1:col,1:row)=0; G_de(1,1:row)=0; al_de(1:col,1:row)=0; k_exp(1:col,1:row)=0;
    for i=1:row
        G_de(1,i)=compand(G_q(1,i)*G_bint,40,1.5,'mu/expander');
        for k=1:col
            G_de(k,i)=compand(g_q(k,i)*g_bint(k,1),40,g_max(k,1),'mu/expander');
            if sig_q(k,i)==0
                G_de(k,i)=-1*G_de(k,i);
            end
            k_exp(k,i)=(1-exp(G_de(k,i)))/(1+exp(G_de(k,i)));
        end
        
        alpha(1:col,1:col)=0;
        for ii=1:col
            alpha(ii,ii)=k_exp(ii,i);
            if (ii > 1)
                for j=1:ii-1
                    alpha(j,ii)=alpha(j,ii-1)-k_exp(ii,i)*alpha(ii-j,ii-1);
                end
            end
        end
        al_de(:,i)=alpha(1:col,col);
    end
    
    % unquantize pitch
    pitch_de(1,1:length(pitch_q))=0;
    for i=1:length(pitch_q)
        if pitch_q(1,i)==0
            pitch_de(1,i)=0;
        else
            pitch_de(1,i)=pitch_q(1,i)+56;
        end
    end
%[x_syn_err_uq,x_err_uq] = syn_err_uq(filename,al,frame_length);
n = length(al);
x_syn_err_uq=[]; 
x_err_uq=[];
for k=1:n
    al_frame=al(:,k);
    x_frame=x(1+(k-1)*R:k*R);
    x_error=filter([1 -al_frame'],1,x_frame);
    x_syn_e_frame=filter(1,[1 -al_frame'],x_error);
    x_err_uq=[x_err_uq x_error'];
    x_syn_err_uq=[x_syn_err_uq x_syn_e_frame'];
end   
%[x_syn_err_q,x_err_q] = syn_err_q(filename,al_de,frame_length);
n = length(al_de);
x_syn_err_q=[];
x_err_q=[];
for k=1:n
    al_de_frame=al_de(:,k);
    x_frame=x(1+(k-1)*R:k*R);
    x_error=x_frame-filter([0 al_de_frame'],1,x_frame);
    x_syn_e_frame=filter(1,[1 -al_de_frame'],x_error);
    x_err_q=[x_err_q x_error'];
    x_syn_err_q=[x_syn_err_q x_syn_e_frame'];
end   
%[x_syn_bi,bi_s] = syn_bi_q(pitch_de,al_de,G_de,frame_length);
n=length(al_de);
    R=frame_length;
    pprev=0;
    bi_s=[];    
    for i=1:n;
        ppd=pitch_de(i);
        if(ppd == 0)
            bi_s=[bi_s randn(1,R)*0.05];
        else
            exc=zeros(1,R);
            if(pprev==0)
                exc(1)=1;
                loc=1;
                while (ppd+loc < R)
                    exc(loc+ppd)=1;
                    loc=loc+ppd;
                end
            else
               loc=loc-R;
               if (loc+ppd <1)
                   loc=1-ppd;
               end
               while(loc+ppd < R)
                   exc(loc+ppd)=1;
                   loc=loc+ppd;
               end
           end
           bi_s=[bi_s exc];
       end
       pprev=ppd;
    end
    
    x_syn_bi=[];
    for k=1:n
        bi_s_frame=bi_s(1+(k-1)*R:k*R);
        al_frame=al_de(:,k);
        X_syn_frame=filter(G_de(:,k),[1 -al_frame'],bi_s_frame);
        x_syn_bi=[x_syn_bi X_syn_frame];    
    end

%% plot
subplot(6,1,1)
temp = x(1+(start_frame-1)*frame_length:(end_frame*frame_length));
temp_err_q = x_syn_err_q(1+(start_frame-1)*frame_length:(end_frame*frame_length));
temp_err_uq = x_syn_err_uq(1+(start_frame-1)*frame_length:(end_frame*frame_length));
temp_bi = x_syn_bi(1+(start_frame-1)*frame_length:(end_frame*frame_length));
plot(1:(end_frame-start_frame+1)*frame_length,temp);
title('s(n)');xlabel('sample');ylabel('amplitude');

if start_frame == end_frame
    subplot(6,1,2)
    [h1,w1] = freqz(1,[1 -al(:,start_frame)']);
    ww = w1*fs/2/pi;
    H1 = 10*log10(abs(h1));
    h2 = freqz(1,[1 -al_de(:,start_frame)']);
    H2 = 10*log10(abs(h2));
    plot(ww,H1,'r',ww,H2,'b');legend('UQ','Q');title('Hq & Huq');
    xlabel('frequency/Hz');ylabel('log magnitude');
end

subplot(6,1,3)
plot(1:(end_frame-start_frame+1)*frame_length,x_err_uq(1+(start_frame-1)*frame_length:(end_frame*frame_length)),'r',1:(end_frame-start_frame+1)*frame_length,x_err_q(1+(start_frame-1)*frame_length:(end_frame*frame_length)),'b');
legend('UQ','Q');xlabel('sample');ylabel('amplitude');title('eHuq(n) & eHq(n)');

subplot(6,1,4)
plot(1:(end_frame-start_frame+1)*frame_length,temp,'--.');hold on
plot(1:(end_frame-start_frame+1)*frame_length,x_syn_err_q(1+(start_frame-1)*frame_length:(end_frame*frame_length)),'r');
legend('s(n)','s~(n)');xlabel('sample');ylabel('amplitude');title('s(n) & s~(n)');

subplot(6,1,5)
plot(bi_s(1+(start_frame-1)*frame_length:(end_frame*frame_length)));
xlabel('sample');ylabel('amplitude');title('binary source');

subplot(6,1,6)
plot(1:(end_frame-start_frame+1)*frame_length,temp,'b');hold on
plot(1:(end_frame-start_frame+1)*frame_length,x_syn_bi(1+(start_frame-1)*frame_length:(end_frame*frame_length))/max(x_syn_bi),'r');
legend('s(n)','s~(n)');xlabel('sample');ylabel('amplitude');title('s(n) & s~(n)');

%% SNR
sSNR_err_q = 10*log10(sum(temp.^2)/sum((temp-temp_err_q').^2))
sSNR_bi = 10*log10(sum(temp.^2)/sum((temp-temp_bi').^2))
wSNR_err_q = 10*log10(sum(x.^2)/sum((x-x_syn_err_q').^2))
wSNR_bi = 10*log10(sum(x.^2)/sum((x-x_syn_bi').^2))


