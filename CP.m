u1 = 1:8; 
h = [0.4 1 0.4];

figure(1)
clf;
subplot(2,1,1)
stem(u1);
axis([0 10 0 10])
title("Input signal")
subplot(2,1,2)
stem(h);
axis([0 10 0 2])
title("Channel impulse response")

N = length(u1);
yl1 = conv(u1,h);
yc1 = cconv(u1,h,N);
figure(2)
clf;
stem(yl1,"x")
hold on;
stem(yc1,"o")
title(["Convolution Results - N=",int2str(N)])
legend ("Linear","Circular","Location","northwest")

L = length(h);      % Length of channel
N = length(u1);     % Length of input signal
ucp = u1(N-L+1:N);  % Use last samples of input signal as the CP
u2 = [ucp u1];      % Prepend the CP to the input signal
yl2 = conv(u2,h);   % Convolution of input+CP and channel
yl2 = yl2(L+1:end); % Remove CP to compare signals
figure(3);
clf;
stem(yc1,"x")
hold on;
stem(yl2,"o")
title("Convolution Results with Cyclic Prefix")
legend ("Linear","Circular","Location","northwest")

if max(yc1 - yl2(1:N)) < 1e-8
    disp("Linear and circular convolution sequences match.")
else
    disp("Received symbols do not match transmitted symbols.")
end