% This MATLAB code was created by Pop

close all; clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 2000;                           % Sampling frequency (Hz)
duration = 2;                       % Signal duration (second)
N = fs*duration;                    % Total number of samples
t = 0:1/fs:duration-1/fs;           % Time vector

% Parameters
a1 = 3;  f1 = 30;  phi1 = 0.6;
a2 = 2;  f2 = 45;  phi2 = -0.8;
a3 = 1;  f3 = 70;  phi3 = 2;
% sinusoidal signals
s1 = a1*cos(2*pi*f1*t + phi1);
s2 = a2*cos(2*pi*f2*t + phi2);
s3 = a3*cos(2*pi*f3*t + phi3);
% Synthesized signal
s = s1 + s2 + s3;       

% QAM Parameters
M = 256;                            % ary-QAM
bitsPerSym = log2(M);
fraction = 4;                       % fixed-point fraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QAM Block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fixedsig = fi(s, 1, bitsPerSym, fraction);
sigbin_tmp = fixedsig.bin;

countbin = 1;
for iterator = 1:(N*bitsPerSym + 3*(N-1))
    if (sigbin_tmp(iterator) ~= ' ')
        tmpbin = sigbin_tmp(iterator);
        sigbin(countbin) = tmpbin;
        sigbindouble(countbin) = str2double(sigbin(countbin));
        countbin = countbin + 1;
    end
end
sigbindouble = sigbindouble.';

txSig = qammod(sigbindouble,M,'bin','InputType','bit','OutputDataType', ...
    numerictype(1,bitsPerSym,fraction),'UnitAveragePower',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IFFT Block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ifftSig = ifft(double(txSig));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Medium Channel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rxSig = awgn(ifftSig, 60);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FFT Block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fftSig = fft(rxSig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% QAM Demodulation Block %%%%%%%%%%%%%%%%%%%%%%%%%%
fixedfftSig = fi(fftSig,1,bitsPerSym,fraction);
demodSig = qamdemod(fixedfftSig, M, 'bin', 'OutputType', 'bit', ...
    'UnitAveragePower', true);
% equal = isequal(sigbindouble,double(demodSig))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% Interpreter Block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uintData = zeros(N,1);
numofbin = 0;
tmpbin = 1;
demodSig = double(demodSig);
for iterator = 1:N
    countbin = 1;
    for count = 1:bitsPerSym 
        temp = demodSig(numofbin+countbin)*(2^(bitsPerSym-countbin));
        uintData(iterator) = uintData(iterator) + temp;
        countbin = countbin + 1;
    end
    numofbin = bitsPerSym*tmpbin;
    tmpbin = tmpbin + 1;
end
uintData = uint8(uintData);
rxData = reinterpretcast(uintData, numerictype(1,bitsPerSym,fraction));
rxData = rxData.';
equal = isequal(double(fixedsig),double(rxData))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots Block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
refQAM = qammod(0:255,256,'UnitAveragePower',true);
cd = comm.ConstellationDiagram('ShowReferenceConstellation', true, ...
    'ReferenceConstellation', refQAM);
cd(fixedfftSig)

figure(1)
subplot(1,2,1)
plot(t,s)
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Time-domain plot')

S = fft(s);
S_oneSide = S(1:N/2);
f = fs*(0:N/2-1)/N;
S_meg = abs(S_oneSide)/(N/2);

subplot(1,2,2)
plot(f,S_meg)
xlabel('Samples');
ylabel('Magnitude');
title('Frequency-domain plot');


figure(2)
subplot(1,2,1)
plot(t,rxData)
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Time-domain plot')

rxS = fft(double(rxData));
rxS_oneSide = rxS(1:N/2);
f = fs*(0:N/2-1)/N;
rxS_meg = abs(rxS_oneSide)/(N/2);

subplot(1,2,2)
plot(f,rxS_meg)
xlabel('Samples');
ylabel('Magnitude');
title('Frequency-domain plot');

clearvars -except bitsPerSym fs t ifftSig fftSig N rxData s cd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
