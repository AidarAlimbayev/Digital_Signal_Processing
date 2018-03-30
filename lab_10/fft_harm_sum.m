N = 6;
t = -1:0.01:1;
A = 1;
T = 1;
nh = (1:N) * 2 - 1;
% arrays strings - harmony waves
harmonics = cos(2 * pi * nh' * t/T);
Am = 2/pi./nh;
Am(2:2:end) = -Am(2:2:end);
s1 = harmonics.*repmat(Am', 1, length(t));
% lines - harmonics sum
s2 = cumsum(s1);
for k = 1:N
    subplot(ceil(N/2), 2, k), plot(t, s2(k,:)),
end

y = fft(s2,201);
figure
stem(t (1:101), y(1:101))