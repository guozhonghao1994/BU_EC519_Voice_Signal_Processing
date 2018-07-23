function [SNR_e,SNR_u] = LP(start_frame,end_frame,mode)
% Steps to produce the figures:
% 1. derive the residual e(n) for the entire file
% 2. replace e(n) with u(n) (for the entire file)
% 3. create G?u(n) (for the entire file)
% 4. create s¡¯(n) with e(n) as excitation
% 5. create s¡¯(n) with G?u(n) as excitation
% start_frame = 1;
% end frame should be no more than 90
if start_frame<0||end_frame>90
    warning('frame # must be with (1,90)');return
end
p_ORDER = [4;10;18]; % Order
[x,fs] = audioread('s5.wav');
load('pp5.mat');
R = 240;
L = 240;
[UN] = IMPULSE(pp5,R);

% p>p(opt)
x_err_18=[]; err_18=[]; x_imp_18=[]; xtotal=[];  binary_source_18=[]; G_binary_source_18=[];
for k=start_frame:end_frame 
    x_frame=x(1+R*(k-1):R*(k-1)+L);
    xtotal=[xtotal x_frame']; 
    UN_frame=UN(1+R*(k-1):R*(k-1)+L); % Extract the impulse train needed
    binary_source_18=[binary_source_18 UN_frame];
    p = p_ORDER(3); % p=18,p>p(opt)
    % Calculate residual error, sythesized speeches
    [x_syn_err_frame, err_frame,x_syn_imp_frame, al,G] = synthesizing(L, p, x_frame,UN_frame);
    G_binary_source_18 = [G_binary_source_18 G*UN_frame];
    err_18=[err_18 (err_frame')]; % Residual error
    x_err_18=[x_err_18 x_syn_err_frame']; % Synthesized speech from residual error
    x_imp_18=[x_imp_18 x_syn_imp_frame]; % Synthesized speech from impulse train
    x_imp_18=x_imp_18*(sqrt(sum(xtotal.^2))/sqrt(sum((x_imp_18).^2)));
    [h18,w18]=freqz(1,[1 al(2:p+1)]); 
end
[f_18,STFT_18,ww_18,HH_18] = H(start_frame,end_frame,p);

% p=p(opt)
x_err_10=[]; err_10=[]; x_imp_10=[]; xtotal=[];binary_source_10=[]; G_binary_source_10=[];
for k=start_frame:end_frame 
    x_frame=x(1+R*(k-1):R*(k-1)+L);
    xtotal=[xtotal x_frame']; 
    UN_frame=UN(1+R*(k-1):R*(k-1)+L);
    binary_source_10=[binary_source_10 UN_frame];
    p=p_ORDER(2); % p=10,p=p(opt)
    [x_syn_err_frame, err_frame,x_syn_imp_frame,al,G] = synthesizing(L, p, x_frame,UN_frame);
    G_binary_source_10 = [G_binary_source_10 G*UN_frame];
    err_10=[err_10 (err_frame')];
    x_err_10=[x_err_10 x_syn_err_frame'];
    x_imp_10=[x_imp_10 x_syn_imp_frame]; 
    x_imp_10=x_imp_10*(sqrt(sum(xtotal.^2))/sqrt(sum((x_imp_10).^2)));
    [h10,w10]=freqz(1,[1 al(2:p+1)]);
end
[f_10,STFT_10,ww_10,HH_10] = H(start_frame,end_frame,p);

% p<p(opt)
x_err_4=[]; err_4=[]; x_imp_4=[]; xtotal=[];binary_source_4=[]; G_binary_source_4=[];
for k=start_frame:end_frame 
    x_frame=x(1+R*(k-1):R*(k-1)+L);
    xtotal=[xtotal x_frame']; 
    UN_frame=UN(1+R*(k-1):R*(k-1)+L);
    binary_source_4 = [binary_source_4 UN_frame];
    p=p_ORDER(1); % p=4,p<p(opt)
    [x_syn_err_frame, err_frame,x_syn_imp_frame,al,G] = synthesizing(L, p, x_frame,UN_frame);
    G_binary_source_4 = [G_binary_source_4 G*UN_frame];
    err_4 = [err_4 (err_frame')];
    x_err_4 = [x_err_4 x_syn_err_frame'];
    x_imp_4 = [x_imp_4 x_syn_imp_frame]; 
    x_imp_4 = x_imp_4*(sqrt(sum(xtotal.^2))/sqrt(sum((x_imp_4).^2)));
    [h4,w4]=freqz(1,[1 al(2:p+1)]);
end
[f_4,STFT_4,ww_4,HH_4] = H(start_frame,end_frame,p);

% s(n)
subplot(7,3,1)
plot(1:1:(end_frame-start_frame+1)*L,xtotal);
title('s(n)');
subplot(7,3,2)
plot(1:1:(end_frame-start_frame+1)*L,xtotal);
title('s(n)');
subplot(7,3,3)
plot(1:1:(end_frame-start_frame+1)*L,xtotal);
title('s(n)');
% H
if mode == 'w' || mode == 'v'|| mode == 'u'
subplot(7,3,4) 
plot(f_4,STFT_4);hold on
plot(ww_4,HH_4)
title('H(n)')
subplot(7,3,5)
plot(f_10,STFT_10);hold on
plot(ww_10,HH_10)
title('H(n)')
subplot(7,3,6)
plot(f_18,STFT_18);hold on
plot(ww_18,HH_18)
title('H(n)')
end
% e(n)
subplot(7,3,7)
plot(1:1:(end_frame-start_frame+1)*L,err_4);
title('e(n)')
subplot(7,3,8)
plot(1:1:(end_frame-start_frame+1)*L,err_10);
title('e(n)')
subplot(7,3,9)
plot(1:1:(end_frame-start_frame+1)*L,err_18);
title('e(n)')
% u(n)
subplot(7,3,10)
plot(1:1:length(binary_source_4),binary_source_4)
title('u(n)')
subplot(7,3,11)
plot(1:1:length(binary_source_10),binary_source_10)
title('u(n)')
subplot(7,3,12)
plot(1:1:length(binary_source_18),binary_source_18)
title('u(n)')
% G*u(n)
subplot(7,3,13)
plot(1:1:length(G_binary_source_4),G_binary_source_4)
title('G*u(n)')
subplot(7,3,14)
plot(1:1:length(G_binary_source_10),G_binary_source_10)
title('G*u(n)')
subplot(7,3,15)
plot(1:1:length(G_binary_source_18),G_binary_source_18)
title('G*u(n)')
% s'(n) by e(n)
subplot(7,3,16)
plot(1:1:(end_frame-start_frame+1)*L,x_err_4);
title("s'(n) by e(n)")
subplot(7,3,17)
plot(1:1:(end_frame-start_frame+1)*L,x_err_10);
title("s'(n) by e(n)")
subplot(7,3,18)
plot(1:1:(end_frame-start_frame+1)*L,x_err_18);
title("s'(n) by e(n)")
% s'(n) by G*u(n)
subplot(7,3,19)
plot(1:1:(end_frame-start_frame+1)*L,x_imp_4);
title("s'(n) by e(n)")
subplot(7,3,20)
plot(1:1:(end_frame-start_frame+1)*L,x_imp_10);
title("s'(n) by e(n)")
subplot(7,3,21)
plot(1:1:(end_frame-start_frame+1)*L,x_imp_18);
title("s'(n) by e(n)")
if mode == 'w'
    suptitle('Entire File')
else if mode == 'v'
        suptitle('Vowel Window')
    else if mode == 'u'
            suptitle('Unvoiced Window')
        else if mode == 'c'
                suptitle('V-F-V Window')
            end
        end
    end
end

% Write synthesized files to wav
if mode == 'w'
    audiowrite('s5(e(n)).wav',x_err_10,fs);
    audiowrite('s5(u(n)).wav',x_imp_10,fs);
end

% SNR
SNR_e_4 = 10*log10((sum(xtotal.^2))/(sum((xtotal-x_err_4).^2)));
SNR_u_4 = 10*log10((sum(xtotal.^2))/(sum((xtotal-x_imp_4).^2)));

SNR_e_10 = 10*log10((sum(xtotal.^2))/(sum((xtotal-x_err_10).^2)));
SNR_u_10 = 10*log10((sum(xtotal.^2))/(sum((xtotal-x_imp_10).^2)));

SNR_e_18 = 10*log10((sum(xtotal.^2))/(sum((xtotal-x_err_18).^2)));
SNR_u_18 = 10*log10((sum(xtotal.^2))/(sum((xtotal-x_imp_18).^2)));

SNR_e = [SNR_e_4 SNR_e_10 SNR_e_18];
SNR_u = [SNR_u_4 SNR_u_10 SNR_u_18];
end

function [f,STFT,ww,HH] = H(start_frame,end_frame,p)
[x,fs] = audioread('s5.wav');
L =240;
data = x((start_frame-1)*L+1:end_frame*L);
% STFT
X = fft(data);
n = 0:(length(X)-1)/2;
f = fs*n/length(X);
X = X(1:length(X)/2);
X = abs(X);
STFT = 10*log10(X/max(X));
% H
L = length(data);
    r = [];
    for i=0:p
        r = [r; sum(data(1:L-i).* data(1+i:L))];
    end
    RR = toeplitz(r(1:p));
    a = inv(RR)*r(2:p+1);
    al = [1;-a];
    %G1 = sqrt(sum(al.*r));
    %data_pre = filter([0 -al(2:end)'],1,data);
    %err = data - data_pre;
    %data_syn_err = filter(1,[1 al(2:end)'],err);
    [h,w] = freqz(1,[1 al(2:p+1)']);
    HH = 10*log10(abs(h)/abs(max(h)));
    ww = w*fs/pi/2;
end

function [UN] = IMPULSE(pitch,R)
    pprev=0;
    UN=[]; 
    for i=1:length(pitch)
        ppd=pitch(i);
        if(ppd == 0)
            UN=[UN randn(1,R)*0.01];
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
           UN=[UN exc];
       end
       pprev=ppd;
    end
end

%SUPTITLE Puts a title above all subplots.
%	SUPTITLE('text') adds text to the top of the figure
%	above all subplots (a "super title"). Use this function
%	after all subplot commands.

% Drea Thomas 6/15/95 drea@mathworks.com

% Warning: If the figure or axis units are non-default, this
% will break.

function hout=suptitle(varargin)

str = varargin{1};
if length(varargin)>1
	params = varargin(2:2:end);
	vals = varargin(3:2:end);
	if length(params)~=length(vals)
		error('Parameter and value pairs are not balanced');
	end
else
	params = {};
	vals = {};
end

% Parameters used to position the supertitle.

% Amount of the figure window devoted to subplots
plotregion = .92;

% Y position of title in normalized coordinates
titleypos  = .95;

% Fontsize for supertitle
fs = get(gcf,'defaultaxesfontsize')+4;

% Fudge factor to adjust y spacing between subplots
fudge=1;

haold = get(gcf,'currentAxes');
figunits = get(gcf,'units');

% Get the (approximate) difference between full height (plot + title
% + xlabel) and bounding rectangle.

	if (~strcmp(figunits,'pixels')),
		set(gcf,'units','pixels');
		pos = get(gcf,'position');
		set(gcf,'units',figunits);
	else
		pos = get(gcf,'position');
	end
	ff = (fs-4)*1.27*5/pos(4)*fudge;

        % The 5 here reflects about 3 characters of height below
        % an axis and 2 above. 1.27 is pixels per point.

% Determine the bounding rectange for all the plots

% h = findobj('Type','axes');   

% findobj is a 4.2 thing.. if you don't have 4.2 comment out
% the next line and uncomment the following block.
	
h = findobj(gcf,'Type','axes');  % Change suggested by Stacy J. Hills

% If you don't have 4.2, use this code instead
%ch = get(gcf,'children');
%h=[];
%for i=1:length(ch),
%  if strcmp(get(ch(i),'type'),'axes'),
%    h=[h,ch(i)];
%  end
%end

	


max_y=0;
min_y=1;

oldtitle =[];
for i=1:length(h),
	if (~strcmp(get(h(i),'Tag'),'suptitle')),
		pos=get(h(i),'pos');
		if (pos(2) < min_y), min_y=pos(2)-ff/5*3;end;
		if (pos(4)+pos(2) > max_y), max_y=pos(4)+pos(2)+ff/5*2;end;
	else
		oldtitle = h(i);
	end
end

if max_y > plotregion,
	scale = (plotregion-min_y)/(max_y-min_y);
	for i=1:length(h),
		pos = get(h(i),'position');
		pos(2) = (pos(2)-min_y)*scale+min_y;
		pos(4) = pos(4)*scale-(1-scale)*ff/5*3;
		set(h(i),'position',pos);
	end
end

np = get(gcf,'nextplot');
set(gcf,'nextplot','add');
if ~isempty(oldtitle)
	delete(oldtitle);
end
ha=axes('pos',[0 1 1 1],'visible','off','Tag','suptitle');
ht=text(.5,titleypos-1,str);set(ht,'horizontalalignment','center','fontsize',fs,'interpreter','none',params,vals);
set(gcf,'nextplot',np);
set(gcf,'currentAxes',haold);
if nargout,
	hout=ht;
end
end


function [x_syn_err_frame, err_frame,x_syn_imp_frame, al,G] = synthesizing(L, p, x_frame,imptr_frame)
al=lpc(x_frame,p);
xpred = filter([0 -al(2:p+1)],1,x_frame); 
err_frame=(x_frame-xpred); 
clear r; r = zeros(p+1,1);
for i=0:p
   r(i+1) = x_frame(1:L-i)' * x_frame(1+i:L); 
end
G=sqrt(al*r);
x_syn_err_frame=filter(1,[1 al(2:p+1)],err_frame); 
x_syn_imp_frame=filter(1,[1 al(2:p+1)],G*imptr_frame); 
end

