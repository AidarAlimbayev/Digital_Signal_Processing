Fs = 1000;          % Sampling frequency
T = 1/Fs;           % Sampling period
L = 1500;           % Length of signal
t = (0:L - 1) * T;  % Time vector

S = 0.7 * sin(2 * pi * 50 * t) + sin(2 * pi * 120 * t);
X = S + 2 * randn(size(t));

plot(1000 * t (1:50), X(1:50))
title('Signal Corrupted with Zero-Mean Randome Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')