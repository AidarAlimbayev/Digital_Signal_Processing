t = [0:0.1:2*pi];
x = sin(t);
% x = [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0];
h = [0, 0, 0.2, 0.5, 0.7, 1, 0.7, 0.5, 0.2, 0, 0];
y1 = conv(x, h)
y2 = conv(h, x)
figure
plot(y1)
figure
plot(y2)
y = conv(x, h, 'same')
y = conv(x, h, 'full')
y = conv(x, h, 'valid')