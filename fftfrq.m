% this code just returns the fft frequencies
% Note N must be divisible by 2

function[frq] = fftfrq(N,Fs)

frq = (Fs/N).*(1:N/2);