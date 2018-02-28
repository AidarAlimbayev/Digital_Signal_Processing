function plot_fir(R,b,Fs)
% Вывод графиков ИХ, АЧХ и ФЧХ КИХ-фильтра 
%  
% R – порядок КИХ-фильтра 
% b – вектор коэффициентов КИХ-фильтра (ИХ КИХ-фильтра) 
% Fs – частота дискретизации (Гц)
%  
% a = [1] – коэффициент знаменателя передаточной функции 
% n – вектор дискретного нормированного времени 
% f – сетка частот (Гц) для расчета АЧХ и ФЧХ 
% H – частотная характеристика
% MAG и PHASE – АЧХ и ФЧХ
%  
a = [1];
n = 0:R;
subplot(3,1,1), stem(n,b,'fill','MarkerSize',3) 
xlabel('n'), title('Impulse Response'), grid
f = 0:((Fs/2)/1000):Fs/2;
H = freqz(b,a,f,Fs);
MAG = abs(H);
PHASE = angle(H);
subplot(3,1,2), plot(f,MAG)
xlabel('f (Hz)'), title('MAGNITUDE'), grid
subplot(3,1,3), plot(f,PHASE)
xlabel('f (Hz)'), title('PHASE'), grid


