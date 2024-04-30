close all; clear; clc;

%% FFT/IFFT with Conjugate Symmetry

% Initial Condition
N = 64;                             % Point of the FFT/IFFT   = 64 points
n_in = (N-2)/2;                     % Number of Input Symbols = 32 symbols
raw_input = randn(n_in, 1, "like", 1i);    % Random Inputs complex symbols
conj_input = conj(raw_input);              % Complex Conjugate of Inputs

% Flip the symbols here
tmp = conj_input;                               
for i = 1:n_in
    conj_input(i) = tmp(n_in-i+1);
end

% Concatenate 0 and the conjugate symbols.
input = cat(1, 0, raw_input, 0, conj_input);
in_re = real(input);
in_im = imag(input);

% IFFT block
ifft_output = ifft(input, N);

ifft_re = real(ifft_output);
ifft_im = imag(ifft_output);

% FFT block
output = fft(ifft_output, N);
out_re = real(output);
out_im = imag(output);

data = output(2:32);
%% Graph plotting
figure(1)
plot(in_re, in_im, 'ro')
hold on;
plot(out_re, out_im, 'blackx')
title("input and output")
legend("input","output")
xlabel("real")
ylabel("imaginary")
hold off;
figure(2)
plot(ifft_re, ifft_im, 'b*')
title("ifft output")
xlabel("real")
ylabel("imaginary")
%figure(3)
%plot()