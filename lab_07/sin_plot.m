% %t = [0:0.1:2*pi];
% a = cos(t);
% b = sin(t);
% figure
% plot(t,a)
% hold on
% plot(t,b)
% %hold on
% %a = [0 10 20 30 40 50 60 70 80 90 100]
% M = mean(a)
% D = var(a)
% E = sum(a.^2)
% P = sum(a.^2)/length(a)
% [acor,lag] = xcorr(a,b);
% R = xcorr(a,b);
% r = xcov(a);
% figure
% plot(lag,acor);
% %hold on
% %figure
% %plot(R)
% %hold on
% %figure
% %plot(r)
% 
% % y1 = rand(50,1);
% % y2 = rand(50,1);
% % 
% % b = xcorr(y1, y2)
% % plot(b)

x = rand(1, 10);
M = mean(x)
D = var(x)
b = xcorr(x);
plot(b)
%plot(x)

