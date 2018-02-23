t = [0:0.1:2*pi];
%a = sin(t);
a = sinc(t);
b = cos(t);
%figure
%plot(b)
%hold on
%plot(a)
%a = [0 10 20 30 40 50 60 70 80 90];
M = mean(a);
E = sum(a.^2);
P = sum(a.^2)/length(a);
R = xcorr(a, b);
r = xcov(a, b);
figure
plot(a)