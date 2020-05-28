function [power, f] = data2power_fft(data, srate)
%data2power_fft Summary of this function goes here
%   Detailed explanation goes here
Fs = srate;         % Sampling frequency                    
% % % T = 1 / Fs;         % Sampling period       
L = length(data);   % Length of signal
% % % t = (0:(L-1)) * T;	% Time vector

n = 2^nextpow2(L);
Y = fft(data, n);

P2 = abs(Y/n);
P1 = P2(1:(n/2+1));
P1(2:(end-1)) = 2 * P1(2:(end-1));

f = Fs * (0:(n/2)) / n;
power = P1;
end

