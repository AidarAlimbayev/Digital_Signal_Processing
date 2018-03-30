%% часть перва€ симул€ци€ сигнала 50 √ц и 120 √ц с шумом
t = 0:0.001:0.6;
x = sin(2 * pi * 50 * t) + sin(2 * pi * 120 * t);
y = x + 2 * randn(size(t));
plot(y(1:50)),grid

%% втора€ часть выделение сигнала с помощью fft
Y = fft(y,512);
Pyy = Y.*conj(Y)/512;
f = 1000 * (0:255)/512;
figure(2), plot(f, Pyy(1:256)), grid

%% часть треть€ симул€ци€ сигнала 50 √ц и 120 √ц с шумом fft2
t = 0:0.001:0.6;
x = sin(2 * pi * 50 * t) + sin(2 * pi * 120 * t);
y1 = x + 2 * randn(size(t));
y2 = x + 2 * randn(size(t));
y = [y1; y2];
plot(y(1:50)),hold on, plot(y(2,1:50)), grid, hold off

%% часть четверта€ выделение сигнала с помощью fft2
Y = fft2(y, 2, 512);
Pyy = Y.*conj(Y)/512;
f = 1000 * (0:255)/512;
figure(2), plot(f, Pyy(1:256)), grid
