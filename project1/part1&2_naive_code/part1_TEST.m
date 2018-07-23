% Part1 - TESTING
% Run part1_TRAIN.m before testing!!!
% Add White Gaussian Noise
NOISE = 0; % 0 for original signal,1 for noised signal
SNR = 10; % Change the SNR here, 0dB,10dB,20dB available
if NOISE == 1
    [data_test,fs] = wavread('test.wav');
    data_temp = awgn(data_test,SNR,0);
    wavwrite(data_temp,fs,'test_noised.wav');
    [data_whole,data_o,LOCATION] = cutfile('test_noised.wav');
end
if NOISE == 0
    [data_whole,data_o,LOCATION] = cutfile('test.wav');
end
THRESHOLD = cell(1,10,2);
for point = 1:length(data_o)
    data_test = data_o{1,point};
    data_test = filter([1,-0.96],1,data_test); % Highpass filter
En_test = [];Zc_test = [];
i = w_length;
while i<=length(data_test)
    En_test = [En_test sum((data_test(i-w_length+1:i).*w_ham(1:w_length)).^2)];
    Zc_test = [Zc_test sum(abs(diff(sign(data_test(i-w_length+1:i)))))];
    i = i+w_shift;
end
fs_test = 10000;
n10msec = fs_test/w_shift;
Zc_test = Zc_test/(2*w_length)*n10msec;
En_test = 10*log(En_test)-max(10*log(En_test));

% Search the start point and endpoint(5 stages)
% B1
flag = 1;
c = 1;
B1_test = 1;
while(flag)
    while(En_test(c)<=ITR_TEST)
        c=c+1;
    end
    B1_test = c;
    flag = 0;
    for c=B1_test+1:B1_test+3
        if c>length(En_test)
            break;
        end
        if En_test(c)<ITU_TEST
            flag = 1;
            break;
        end
        if flag
            c = B1_test+1;
        else 
            break;
        end
    end
end
% E1
flag = 1;
c = length(En_test);
E1_test = c;
while(flag)
    while(En_test(c)<=ITR_TEST)
        c=c-1;
    end
    E1_test = c;
    flag = 0;
    for c=E1_test+1:-1:E1_test-3
        if c>length(length(En_test))
            break;
        end
        if En_test(c)<ITU_TEST
            flag = 1;
            break;
        end
    end
    if flag
        c=E1_test-1;
    else
        break;
    end
end
% B2
for i=B1_test:-1:B1_test-10
    if i<1
        break;
    end
    sumZ = 0;
    ind = [];
    if Zc_test(i)>IZCT_TEST
        sumZ = sumZ+1;
        ind = [i ind];
    end
end
if sumZ>=4
    B2_test = ind(i);
else
    B2_test = B1_test;
end
% E2
for i=E1_test:E1_test+10
    if i>length(En_test)
        break;
    end
    sumZ = 0;
    ind = [];
    if Zc_test(i)>IZCT_TEST
        sumZ = sumZ+1;
        ind = [ind i];
    end
end
if sumZ>=4
    E2_test = ind(end);
else
    E2_test = E1_test;
end
% Final Check
i = B2_test-1;
while(1)
    if i<1
        i = 1;
        break
    end
    if En_test(i)>ITR_TEST
        B2_test = i;
    else
        break;
    end
    i = i-1;
end
i = E2_test+1;
while(1)
    if i>length(En_test)
        i = length(En_test);
        break;
    end
    if En_test(i)>ITR_TEST
        E2_test = i;
    else
        break;
    end
    i = i+1;
end

% Plot the En/Zc of each phoneme 
% figure
% plot(En_test,'b')
% hold on
% plot(1:length(En_test),ITU_TEST,'r');hold on
% plot(1:length(En_test),ITR_TEST,'g');hold on
% ylim=get(gca,'Ylim');
% plot([B1_test,B1_test],ylim,'r--','LineWidth',2);hold on
% plot([E1_test,E1_test],ylim,'r--','LineWidth',2)
% xlabel('frame')
% ylabel('Er_test')
% if NOISE == 1
%     title([testname(1:4) '  WGN:' num2str(SNR) 'dB'])
% else
%     title(testname(1:4))
% end
% 
% figure
% plot(Zc_test,'r');hold on
% plot(1:length(Zc_test),IZCT_TEST,'g');hold on
% ylim=get(gca,'Ylim');
% plot([B2_test,B2_test],ylim,'b--','LineWidth',2);hold on
% plot([E2_test,E2_test],ylim,'b--','LineWidth',2);
% xlabel('frame')
% ylabel('Zr_test')
% if NOISE == 1
%     title([testname(1:4) '  WGN:' num2str(SNR) 'dB'])
% else
%     title(testname(1:4))
% end

% Fix the start/end point and plot V/UV parts
startpoint = min(B1_test,B2_test)*w_shift;
endpoint = max(E1_test,E2_test)*w_shift;
% figure
% plot(data_test);hold on
% data_v_test = data_test(startpoint:endpoint);
% plot(startpoint:endpoint,data_v_test,'r')
% xlabel('sample')
% ylabel('amplitude')
% legend('Original','V Region')
% if NOISE == 1
%     title([testname(1:4) '  Fs=11kHz  Resample=10kHz' '  WGN:' num2str(SNR) 'dB'])
% else
%     title([testname(1:4) '  Fs=11kHz  Resample=10kHz'])
% end

% Store startpoint and endpoint in THRESHOLD
THRESHOLD{1,point,1} = startpoint;
THRESHOLD{1,point,2} = endpoint;
end

% Plot V/UNV of the test case
figure
plot(data_whole,'b');hold on
tempdata = data_whole(THRESHOLD{1,1,1}:THRESHOLD{1,1,2});
plot(THRESHOLD{1,1,1}:THRESHOLD{1,1,2},tempdata,'r');hold on
for k=1:length(LOCATION)
        tempdata = data_whole(LOCATION(k)+THRESHOLD{1,k+1,1}:LOCATION(k)+THRESHOLD{1,k+1,2});
        plot(LOCATION(k)+THRESHOLD{1,k+1,1}:LOCATION(k)+THRESHOLD{1,k+1,2},tempdata,'r');hold on
end
xlabel('sample')
ylabel('amplitude')
legend('Original','V Region')
if NOISE == 1
     title(['test' '  Fs=11kHz  Resample=10kHz' '  WGN:' num2str(SNR) 'dB'])
else
     title(['test' '  Fs=11kHz  Resample=10kHz'])
end
