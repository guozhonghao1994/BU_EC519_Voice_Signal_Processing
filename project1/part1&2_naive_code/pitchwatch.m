function [key] = pitchwatch(x,Fs)
% Set constants.
xlen = length(x);
wlen = 400;
wlag = 110;

% Zero-pad signal.
xpad = wlag-rem(xlen-wlen-1,wlag)-1;
x = [x(:); zeros(xpad,1)];
L = (xlen+xpad-wlen)/wlag+1;
L = L-wlen/wlag/2;

% Calculate the AMDF and AMDF-fractional pitch periods.
for k1 = 1:L
    k2 = (k1-1)*wlag;
    for k3 = 1:wlen/2-1, 
        c(k3) = sum(abs(x(k2+(1:wlen))-x(k2+k3+(1:wlen)))); 
    end
    n = findpeaks(-c);
    if ~isempty(n)
        n(find(c(n) > mean(c(n))-sqrt(var(c(n))))) = []; 
    end
    if ~isempty(n) 
        t0 = n(1);
    else           
        t0 = []; 
    end
   
   if ~isempty(t0)
      u1 = x(k2+t0+(1:wlen))-x(k2+(1:wlen));
      u2 = x(k2+t0+(1:wlen))-x(k2+t0+1+(1:wlen));
      t1 = sum(u1.*u2)/sum(u2.*u2);
      y(k1) = Fs/(t0+t1);
   else
      y(k1) = NaN;
   end
end

% Plot pitch keys.
key = 35*log2(y/440)+180;
figure
plot((0.5+(0:L-1))*wlag/Fs,key,'b.');
ylim([70 500])
xlim([0 length(x)/Fs])
xlabel('Time')
ylabel('Pitch Frequency')

% FUNCTIONS

function n = findpeaks(x)
n    = find(diff(diff(x) > 0) < 0);
u    = find(x(n+1) > x(n));
n(u) = n(u)+1;
