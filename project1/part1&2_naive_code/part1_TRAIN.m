% Part1 - TRAINING
% Start, Endpoint, V-UNV Detection
close all
clear
clc
% Input the audio file and split them apart by phoneme
filename = ['train1.wav';'train2.wav';'train3.wav';'train4.wav';'train5.wav';'train6.wav';'train7.wav'];
ITU_TRAIN = cell(7,10);
ITR_TRAIN = cell(7,10);
IZCT_TRAIN = cell(7,10);
THRESHOLD = cell(7,10,2);
% For each file, read and resample
for num=1:size(filename,1)
    [data_whole,data_o,LOCATION] = cutfile(filename(num,:)); % cut the num_th wav file to cells
for point = 1:length(data_o)
    data = data_o{1,point}; % extract each phoneme in the cell
% Design a highpass filter and act on data
data = filter([1,-0.9],1,data);

% Hamming window
w_length = 40;
w_shift = 10;
w_ham = hamming(w_length);

% Calculate the En and ZC
En = [];Zc = [];
i = w_length;
while i<=length(data)
    En = [En sum((data(i-w_length+1:i).*w_ham(1:w_length)).^2)];
    Zc = [Zc sum(abs(diff(sign(data(i-w_length+1:i)))))];
    i = i+w_shift;
end
fs = 10000;
n10msec = fs/w_shift;
Zc = Zc/(2*w_length)*n10msec;
En = 10*log(En)-max(10*log(En));

% Get the average and variation
% Suppose the first 10ms(1 frame) are silence parts
eavg = mean(En(1:1));
zcavg = mean(Zc(1:1));
esig = std(En(1:1));
zcsig = std(Zc(1:1));

% Get the IZCT,ITU,ITR and store them in to array
IF = 35;
IZCT = max(IF,zcavg+3*zcsig);
ITU = -20;
ITR = max(ITU-10,eavg+3*esig);
ITU_TRAIN{num,point} = ITU;
ITR_TRAIN{num,point} = ITR;
IZCT_TRAIN{num,point} = IZCT;


% Search the start point and endpoint(5 stages)
% B1
flag = 1;
c = 1;
B1 = 1;
while(flag)
    while(En(c)<=ITR)
        c=c+1;
    end
    B1 = c;
    flag = 0;
    for c=B1+1:B1+3
        if c>length(En)
            break;
        end
        if En(c)<ITU
            flag = 1;
            break;
        end
        if flag
            c = B1+1;
        else 
            break;
        end
    end
    B1 = c;
end
% E1
flag = 1;
c = length(En);
E1 = c;
while(flag)
    while(En(c)<=ITR)
        c=c-1;
    end
    E1 = c;
    flag = 0;
    for c=E1+1:-1:E1-3
        if c>length(length(En))
            break;
        end
        if En(c)<ITU
            flag = 1;
            break;
        end
    end
    if flag
        c=E1-1;
    else
        break;
    end
    E1 = c;
end
% B2
for i=B1:-1:B1-10
    if i<1
        break;
    end
    sumZ = 0;
    ind = [];
    if Zc(i)>IZCT
        sumZ = sumZ+1;
        ind = [i ind];
    end
end
if sumZ>=4
    B2 = ind(i);
else
    B2 = B1;
end
% E2
for i=E1:E1+10
    if i>length(En)
        break;
    end
    sumZ = 0;
    ind = [];
    if Zc(i)>IZCT
        sumZ = sumZ+1;
        ind = [ind i];
    end
end
if sumZ>=4
    E2 = ind(end);
else
    E2 = E1;
end
% Final Check
i = B2-1;
while(1)
    if i<1
        i = 1;
        break
    end
    if En(i)>ITR
        B2 = i;
    else
        break;
    end
    i = i-1;
end
i = E2+1;
while(1)
    if i>length(En)
        i = length(En);
        break;
    end
    if En(i)>ITR
        E2 = i;
    else
        break;
    end
    i = i+1;
end

% Plot the En/Zc of each phoneme
% figure
% plot(En,'b')
% hold on
% plot(1:length(En),ITU,'r');hold on
% plot(1:length(En),ITR,'g');hold on
% ylim=get(gca,'Ylim');
% plot([B1,B1],ylim,'r--','LineWidth',2);hold on
% plot([E1,E1],ylim,'r--','LineWidth',2)
% xlabel('frame')
% ylabel('Er')
% 
% figure
% plot(Zc,'r');hold on
% plot(1:length(Zc),IZCT,'g');hold on
% ylim=get(gca,'Ylim');
% plot([B2,B2],ylim,'b--','LineWidth',2);hold on
% plot([E2,E2],ylim,'b--','LineWidth',2);
% xlabel('frame')
% ylabel('Zr')

% Fix the start/end point of each phoneme
startpoint_v = B1*w_shift;
endpoint_v = E1*w_shift;
startpoint = min(B1,B2)*w_shift;
endpoint = max(E1,E2)*w_shift;
% figure
% plot(data);hold on
% data_v = data(startpoint:endpoint);
% plot(startpoint:endpoint,data_v,'r')
% xlabel('sample')
% ylabel('amplitude')
% legend('Original','V Region')

% Store startpoint and endpoint in THRESHOLD
THRESHOLD{num,point,1} = startpoint;
THRESHOLD{num,point,2} = endpoint;

    end
    % Plot the original audio file
    figure
    plot(data_whole,'b');hold on
    tempdata = data_whole(THRESHOLD{num,1,1}:THRESHOLD{num,1,2});
    plot(THRESHOLD{num,1,1}:THRESHOLD{num,1,2},tempdata,'r');hold on
    for k=1:length(LOCATION)
        tempdata = data_whole(LOCATION(k)+THRESHOLD{num,k+1,1}:LOCATION(k)+THRESHOLD{num,k+1,2});
        plot(LOCATION(k)+THRESHOLD{num,k+1,1}:LOCATION(k)+THRESHOLD{num,k+1,2},tempdata,'r');hold on
    end
    
    xlabel('sample')
    ylabel('amplitude')
    legend('original signal','voiced/uvoiced parts')
    title(filename(num,1:6))
    
end
figure
freqz([1,-0.9],1);

% Get the estimated parameter(ITU,ITR,IZCT) for testing case
ITU_TEST = [];
for m=1:size(filename,1)
    for n = 1:10
        ITU_TEST = [ITU_TEST ITU_TRAIN{m,n}];
    end
end
ITU_TEST = mean(ITU_TEST);

ITR_TEST = [];
for m=1:size(filename,1)
    for n = 1:10
        ITR_TEST = [ITR_TEST ITR_TRAIN{m,n}];
    end
end
ITR_TEST = mean(ITR_TEST);

IZCT_TEST = [];
for m=1:size(filename,1)
    for n = 1:10
        IZCT_TEST = [IZCT_TEST IZCT_TRAIN{m,n}];
    end
end
IZCT_TEST = mean(IZCT_TEST);



















