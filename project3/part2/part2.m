% Part 2 - CELP
close all
clear
clc
filename = 'H.2.wav';
[X,fs] = audioread(filename); % H.1.wav/H.2.wav
start_frame = 67;
end_frame = 67; % For s5.wav, no more than 90
frame_length = 330; % frame length in samples
info_rate = 4800; % 4800/9600
bits_lpc = [5 5 5 5 4 4 4 4 3 2];
switch info_rate
    case 4800
        bits_adpcb=7; % Pitch delay bits, fixed due to the range of humna pitch is 2^7=128
        bits_stocb=8; % Stochastic codebook bits
        bits_gain=5; % Stochastic and adaptive codebook gain bits
    case 9600
        bits_adpcb=7;
        bits_stocb=17;
        bits_gain=18;
end

%% Generate Codebook
rng('default');
stocb = randn(60,2^bits_stocb);
%% CELP Coding
    R=frame_length;
    n=floor((length(X))/R);
    if strcmp(filename,'s5')
        n = n-10;
    end
    
    % Prepare to record data for transmition
    adp_result(1:4,1:n)=0; adp_gain_result(1:4,1:n)=0;
    sto_result(1:4,1:n)=0; sto_gain_result(1:4,1:n)=0;
    
    % Preparef ro adaptve coding (LTP)
    adpcb_enc(1:frame_length,1)=0;
    
    % Prepare for speech synthesizing
    al(1:10,1:n)=0;
    
    % Record data for performance plot
    X_syn_uq(1:frame_length*n,1)=0; X_syn_uq_frame(1:frame_length,1)=0;
    d1_uq(1:frame_length*n,1)=0; e1_uq(1:frame_length*n,1)=0;
    d2_uq(1:frame_length*n,1)=0; e2_uq(1:frame_length*n,1)=0;
   
    for i=1:n
        X_frame=X(1+(i-1)*R:i*R);
    
        % Find VT area coefficient of VT
        al_frameu=lpc(X_frame,10);
        al_frame=-al_frameu(2:11)';
        al(:,i)=al_frame';
        
        % d[n], residual error from original speech
        d1_uq(1+(i-1)*frame_length:frame_length*i,1)=filter([1 -al_frame'],1,X_frame);

        % Zeroes for perceptual weighting filter
        al_perc=al_frame*0.85; % Zeroes for perceptual weighting filter
        
        % Prepare recording optimal value for subframes
        adp_result_frame(1:4,1)=0; adp_gain_result_frame(1:4,1)=0;
        sto_result_frame(1:4,1)=0; sto_gain_result_frame(1:4,1)=0; 
        e1_uq_frame(1:frame_length,1)=0; d2_uq_frame(1:frame_length,1)=0; e2_uq_frame(1:frame_length,1)=0;
        
        % LTP and STP
        for j=1:4
% STP
            X_ref=X_frame(1+(j-1)*60:j*60);
            
            % Generate adaptive codebook for current subframe
            % Reference: http://www.mathworks.com/matlabcentral/fileexchange/39038-celp-codec
            adpcb_subf(:,1:2^bits_adpcb)=toeplitz(adpcb_enc(2^bits_adpcb:2^bits_adpcb+59), ...
                flipud(adpcb_enc(1:2^bits_adpcb)));
            adpcb_ref=filter(1,[1 -al_frame'],adpcb_subf);
            
            % Find best match entry of adaptive codebook index
            % I am not using minimum MSE, cross corelation is used here
            adpcb_eng=sum(adpcb_ref.^2);
            adpcb_corr=X_ref'*adpcb_ref;
            kk = find(adpcb_corr == max(adpcb_corr));
            adpcb_index_subf = kk(1)+frame_length-2^bits_adpcb; 
            
            % Find best match adaptice codebook gain
            % Reference: http://www.mathworks.com/matlabcentral/fileexchange/39038-celp-codec
            adpgain_subf=abs(adpcb_corr(kk(1))/(adpcb_eng(kk(1))+10*eps));
            if adpgain_subf>1.4
              adpgain_subf=1.4;
            end
            
            % Final excitation in current subframe's STP
            excit_subf = adpgain_subf*adpcb_enc(frame_length-adpcb_index_subf+1:frame_length-adpcb_index_subf+60);
            
            % e[n], coding error past to LTP
            e1_uq_frame(1+(j-1)*60:60*j)=excit_subf;
            
% LTP
            % Prepare reference for LTP
            X_subf=X_ref-filter(1,[1 -al_frame'],excit_subf);
            
            % Apply perceptual weighting filter
            %X_kk=filter([1 -al_frame'],1,X_subf);
            %X_kk=filter([1 -al_frame'],[1 -al_perc'],X_kk);
            %X_ref=filter(1,[1 -al_frame'],X_kk);
            X_ref=X_subf;
            
            % Prepare stochastic codebook
            stocb_ref = filter(1,[1 -al_frame'],stocb);
            
            % Find best match entry of stochastic codebook index
            % Reference: http://markus-hauenstein.de/sigpro/codec/codec.shtml
            stocb_eng=sum(stocb_ref.^2);
            stocb_corr=X_ref'*stocb_ref;
            stocb_index_subf=find(stocb_corr==max(stocb_corr));
            
            % Find best match gain
            stogain_subf=stocb_corr(stocb_index_subf(1))/stocb_eng(stocb_index_subf(1));
             
            % e[n]', best match in stochastic codebook
            e2_uq_frame(1+(j-1)*60:60*j)=stogain_subf*stocb(:,stocb_index_subf(1));
            
            % Generate best match excitation for current subframe
            excit_subf = stogain_subf*stocb(:,stocb_index_subf(1))+adpgain_subf*adpcb_enc(frame_length ...
                  -adpcb_index_subf+1:frame_length-adpcb_index_subf+60);
                        
            % d[n]', final output from CELP
            d2_uq_frame(1+(j-1)*60:60*j)=excit_subf;
            
            % Update adaptive codebook for next frame
            adpcb_enc = [adpcb_enc(61:frame_length); excit_subf];    
            
            % Record for current subframe
            adp_result_frame(j)=adpcb_index_subf;
            adp_gain_result_frame(j)=adpgain_subf;
            sto_result_frame(j)=stocb_index_subf(1);
            sto_gain_result_frame(j)=stogain_subf;
            
            % Synthesize unquantized speech
            X_syn_uq_frame(1+(j-1)*60:60*j)=filter(1,[1 -al_frame'],excit_subf);
            
        end
        
        X_syn_uq(1+(i-1)*frame_length:frame_length*i,1)=X_syn_uq_frame;
        e1_uq(1+(i-1)*frame_length:frame_length*i,1)=e1_uq_frame;
        d2_uq(1+(i-1)*frame_length:frame_length*i,1)=d2_uq_frame; 
        e2_uq(1+(i-1)*frame_length:frame_length*i,1)=e2_uq_frame;
        
        % Record for current frame
        adp_result(1:4,i)=adp_result_frame; adp_gain_result(1:4,i)=adp_gain_result_frame;
        sto_result(1:4,i)=sto_result_frame; sto_gain_result(1:4,i)=sto_gain_result_frame;
      
    end

%% Quantization
g_max=[5 3.5 2.5 1.5 1.5 2.5 1 1.5 1.5 1]';
    p=10; 
    g_uqua=[]; sig_qua=[]; 
    adp_gain_result_qua(1:4,1:n)=0; sto_gain_result_qua(1:4,1:n)=0;
    alpha(1:p,1:p)=0; kk(1:p,1)=0;
    
    % Frame by frame processong
    for k=1:n
    al_frame=al(:,k);
    
    alpha(1:p,1:p)=0; alpha(1:p,p)=al_frame; 
    kk(1:p,1)=0; kk(p)=alpha(p,p);
    for i=p:-1:2
        for j=1:i-1
            alpha(j,i-1)=(alpha(j,i)+kk(i)*alpha(i-j,i))/(1-kk(i).^2);
        end
        kk(i-1)=alpha(i-1,i-1);
    end
    
    % Convert PARCOR to log area ratio
    g(1:p,1)=0;
    for i=1:p
        g(i,1)=log((1-kk(i,1))/(1+kk(i,1)));
    end
    
    % Produce sign array
    sig_frame(1:p,1)=0;
    for i=1:p
       if g(i)<0
           sig_frame(i)=0;
       else
           sig_frame(i)=1;
       end
    end

    % Produce sign and vocal track coefficient matrix, unquantized
    sig_qua=[sig_qua sig_frame];
    g_uqua=[g_uqua g];
    end
    
    % Quantize vocal track coefficient
    g_qua(p,1:n)=0; k_g(1:p,1:n)=0;
    for i=1:n
        for k=1:p
            k_g(k,i)=compand(abs(g_uqua(k,i)),40,g_max(k,1),'mu/compressor'); % Log companding
            g_qua(k,i)=quantiz(k_g(k,i),0:g_max(k,1)/2^(bits_lpc(1,k)-1):g_max(k,1));
        end
    end
    
% Quantize adaptive and stochastic codebook gains adp_gain_result,sto_gain_result
    for k=1:n
        for i=1:4
            %k_g(i,k)=compand(adp_gain_result(i,k),40,1.4,'mu/compressor'); % Log companding
            %adp_gain_result_qua(i,k)=quantiz(k_g(i,k),0:1.4/2^bits_gain:1.4);
            adp_gain_result_qua(i,k)=quantiz(adp_gain_result(i,k),0:1.4/2^bits_gain:1.4);
            sto_gain_result_qua(i,k)=quantiz(sto_gain_result(i,k),-0.06:0.12/2^bits_gain:0.06);
        end
    end
    
%% Unquantization
g_max=[5 3.5 2.5 1.5 1.5 2.5 1 1.5 1.5 1]';
    [pp,nn]=size(g_qua);
    g_bint(1:pp,1)=0;
    for i=1:pp
        g_bint(i,1)=g_max(i,1)/2^(bits_lpc(1,i)-1);
    end
    
    % Unquantization
    k_exp(1:pp,1:nn)=0; g_exp(1:pp,1:nn)=0; al_exp(1:pp,1:nn)=0; k_exp(1:pp,1:nn)=0;
    adp_gain_result_exp(1:4,1:nn)=0; sto_gain_result_exp(1:4,1:nn)=0;
    for i=1:nn
        % Unquzntize, expand g, convert g to PARCOR
        for k=1:pp
            g_exp(k,i)=compand(g_qua(k,i)*g_bint(k,1),40,g_max(k,1),'mu/expander');
            if sig_qua(k,i)==0
                g_exp(k,i)=-1*g_exp(k,i);
            end
            k_exp(k,i)=(1-exp(g_exp(k,i)))/(1+exp(g_exp(k,i)));
        end
        
        % PARCOR to polynomial corfficients
        alpha(1:pp,1:pp)=0;
        for ii=1:pp
            alpha(ii,ii)=k_exp(ii,i);
            if (ii > 1)
                for j=1:ii-1
                    alpha(j,ii)=alpha(j,ii-1)-k_exp(ii,i)*alpha(ii-j,ii-1);
                end
            end
        end
        al_exp(:,i)=alpha(1:pp,pp);
    end
    
% Unquantize adaptive and stochastic codebook gains adp_gain_result,sto_gain_result
    for k=1:nn
        for i=1:4
            %adp_gain_result_exp(i,k)=compand(adp_gain_result_qua(i,k)*1.4/2^bits_gain,40,1.4,'mu/expander');
            adp_gain_result_exp(i,k)=adp_gain_result_qua(i,k)*1.4/2^bits_gain;
            sto_gain_result_exp(i,k)=sto_gain_result_qua(i,k)*0.12/2^bits_gain-0.06;
        end
    end
    
%% Synthesizing
 X_syn_q(1:frame_length*n,1)=0; adpcb_dec(1:frame_length,1)=0;
    d1_q_frame(1:frame_length,1)=0; e1_q_frame(1:frame_length,1)=0;
    d2_q_frame(1:frame_length,1)=0; e2_q_frame(1:frame_length,1)=0;
    
    % Synthesizing
    for i=1:n
        X_syn_q_frame(1:frame_length,1)=0;
        excit_subf(1:60,1)=0;
        for j=1:4
            % e[n]
            e1_q_frame(1+(j-1)*60:60*j)=adp_gain_result_exp(j,i)*adpcb_dec(frame_length-adp_result(j,i)+1: ...
                frame_length-adp_result(j,i)+60);
            
            % e[n]'
            e2_q_frame(1+(j-1)*60:60*j)=sto_gain_result_exp(j,i)*stocb(:,sto_result(j,i));
            
            % Produce excitations using codebooks and gains
            excit_subf=sto_gain_result_exp(j,i)*stocb(:,sto_result(j,i))+ ...
                adp_gain_result_exp(j,i)*adpcb_dec(frame_length-adp_result(j,i)+1:frame_length-adp_result(j,i)+60);
            
            % d[n]' d[n]
            d2_q_frame(1+(j-1)*60:60*j)=excit_subf;
            d1_q_frame(1+(j-1)*60:60*j)=excit_subf;
            
            % Update adaptive codebook for next subframe
            adpcb_dec=[adpcb_dec(61:frame_length); excit_subf];
            
            % Generate synthesized speech
            X_syn_q_frame(1+(j-1)*60:60*j)=filter(1,[1 -al_exp(:,i)'],excit_subf);
        end
        
        d1_q(1+(i-1)*frame_length:frame_length*i,1)=d1_q_frame; e1_q(1+(i-1)*frame_length:frame_length*i,1)=e1_q_frame;
        d2_q(1+(i-1)*frame_length:frame_length*i,1)=d2_q_frame; e2_q(1+(i-1)*frame_length:frame_length*i,1)=e2_q_frame;
        X_syn_q(1+(i-1)*frame_length:frame_length*i,1)=X_syn_q_frame;
        
    end
    
%% Plot
subplot(7,1,1)
XX = (end_frame-start_frame+1)*frame_length;
temp = X((start_frame-1)*frame_length+1:end_frame*frame_length);
plot(1:XX,temp,'LineWidth',1.5);
xlabel('sample')
ylabel('amplitude')
title('s(n)')

if end_frame == start_frame
    subplot(7,1,2)
[h1,w1] = freqz(1,[1,-al(:,end_frame-start_frame+1)']);
ww = w1*fs/pi/2;
H1 = 10*log10(abs(h1));
h2 = freqz(1,[1 -al_exp(:,end_frame-start_frame+1)']);
H2 = 10*log10(abs(h2));
plot(ww,H1,'b',ww,H2,'r');
title('Hq')
legend('UQ','Q')
xlabel('frequency')
ylabel('amplitude')
end

subplot(7,1,3)
plot(1:XX,d1_uq((start_frame-1)*frame_length+1:end_frame*frame_length),'b',...
    1:XX,d1_q((start_frame-1)*frame_length+1:end_frame*frame_length),'r');
legend('UQ','Q')
title('d(n)')
xlabel('sample')
ylabel('amplitude')

subplot(7,1,4)
plot(1:XX,e1_uq((start_frame-1)*frame_length+1:end_frame*frame_length),'b',...
    1:XX,e1_q((start_frame-1)*frame_length+1:end_frame*frame_length),'r');
legend('UQ','Q')
title('e(n)')
xlabel('sample')
ylabel('amplitude')

subplot(7,1,5)
plot(1:XX,d2_uq((start_frame-1)*frame_length+1:end_frame*frame_length),'b',...
    1:XX,d2_q((start_frame-1)*frame_length+1:end_frame*frame_length),'r');
legend('UQ','Q')
title('n_(n)')
xlabel('sample')
ylabel('amplitude')

subplot(7,1,6)
plot(1:XX,e2_uq((start_frame-1)*frame_length+1:end_frame*frame_length),'b',...
    1:XX,e2_q((start_frame-1)*frame_length+1:end_frame*frame_length),'r');
legend('UQ','Q')
title('d_(n)')
xlabel('sample')
ylabel('amplitude')

subplot(7,1,7)
plot(1:XX,X_syn_uq((start_frame-1)*frame_length+1:end_frame*frame_length),'b',...
    1:XX,X_syn_q((start_frame-1)*frame_length+1:end_frame*frame_length),'r');
xlabel('sample')
ylabel('amplitude')
legend('UQ','Q')
title('s~(n)uq & s~(n)q')

%% SNR
wSNR_uq=10*log10(sum(X(1:length(X_syn_uq),1).^2)/sum((X(1:length(X_syn_q),1)-X_syn_uq).^2))
wSNR_q=10*log10(sum(X(1:length(X_syn_q),1).^2)/sum((X(1:length(X_syn_q),1)-X_syn_q).^2))
sSNR_uq=10*log10(sum(X(1+(start_frame-1)*frame_length:end_frame*frame_length,1).^2)/sum((X(1+(start_frame-1)*frame_length:end_frame*frame_length,1)-X_syn_uq(1+(start_frame-1)*frame_length:end_frame*frame_length,1)).^2))
sSNR_q=10*log10(sum(X(1+(start_frame-1)*frame_length:end_frame*frame_length,1).^2)/sum((X(1+(start_frame-1)*frame_length:end_frame*frame_length,1)-X_syn_q(1+(start_frame-1)*frame_length:end_frame*frame_length,1)).^2))