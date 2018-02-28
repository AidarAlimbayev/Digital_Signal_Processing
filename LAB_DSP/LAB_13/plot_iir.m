function plot_iir(b,a,Fs)
% Вывод графиков АЧХ, ФЧХ, ИХ и карты нулей и полюсов БИХ-фильтра 
%  
% b – вектор коэффициентов числителя передаточной функции 
% a – вектор коэффициентов знаменателя передаточной функции 
% Fs – частота дискретизации (Гц)
%  
% M - длина ИХ БИХ-фильтра, ограниченная до 50-ти отсчетов
% n – вектор дискретного нормированного времени 
% h – вектор отсчетов ИХ 
% f – сетка частот (Гц) для расчета АЧХ и ФЧХ 
% H – частотная характеристика
% MAG и PHASE – АЧХ и ФЧХ
%  
M = 50; 
n = 0:(M-1);
h = impz(b,a,M);
f = 0:((Fs/2)/1000):Fs/2;      
H = freqz(b,a,f,Fs); 
MAG = abs(H);
PHASE = phase(H);
subplot(2,2,1), plot(f,MAG), xlabel('f (Hz)')
title('MAGNITUDE'), grid, ylim([0 1.2])
subplot(2,2,2), zplane(b,a), title('Z-plane zero-pole plot'), grid
subplot(2,2,3), plot(f,PHASE), xlabel('f (Hz)')
title('PHASE'), grid
subplot(2,2,4), stem(n,h,'fill','MarkerSize',3)
xlabel('n'), title('Impulse Response'), grid








