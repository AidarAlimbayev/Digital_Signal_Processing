script
clc
clear
disp('% ËĞ ¹7. ÄÈÑÊĞÅÒÍÛÅ ÑÈÃÍÀËÛ')
disp('%')
disp('%')
disp('% Ââåäèòå ÈÑÕÎÄÍÛÅ ÄÀÍÍÛÅ');
DATA=0;
while DATA==0
Nb = input('Nb = ');            % ÍÎÌÅĞ ÁĞÈÃÀÄÛ
N = input('N = ');              % ÄËÈÍÀ ÏÎÑËÅÄÎÂÀÒÅËÜÍÎÑÒÈ
T = input('T = ');              % ÏÅĞÈÎÄ ÄÈÑÊĞÅÒÈÇÀÖÈÈ
a = input('a = ');              % ÎÑÍÎÂÀÍÈÅ ÄÈÑÊĞÅÒÍÎÉ İÊÑÏÎÍÅÍÒÛ
C = input('C = ');      % ÀÌÏËÈÒÓÄÀ ÄÈÑÊĞÅÒÍÎÃÎ ÃÀĞÌÎÍÈ×ÅÑÊÎÃÎ ÑÈÃÍÀËÀ
w0 = input('w0 = ');    % ×ÀÑÒÎÒÀ ÄÈÑÊĞÅÒÍÎÃÎ ÃÀĞÌÎÍÈ×ÅÑÊÎÃÎ ÑÈÃÍÀËÀ
m = input('m = ');              % ÂÅËÈ×ÈÍÀ ÇÀÄÅĞÆÊÈ
U = input('U = ');              % ÀÌÏËÈÒÓÄÀ ÈÌÏÓËÜÑÀ
n0 = input('n0 = ');            % ÌÎÌÅÍÒ ÍÀ×ÀËÀ ÈÌÏÓËÜÑÀ
n_imp = input('n_imp = ');      % ÄËÈÒÅËÜÍÎÑÒÜ ÈÌÏÓËÜÑÀ
B = input('B = ');              % ÂÅÊÒÎĞ ÀÌÏËÈÒÓÄ
w = input('w = ');              % ÂÅÊÒÎĞ ×ÀÑÒÎÒ 
A = input('A = ');        % ÂÅÊÒÎĞ ÊÎİÔÔÈÖÈÅÍÒÎÂ ËÈÍÅÉÍÎÉ ÊÎÌÁÈÍÀÖÈÈ
Mean = input('Mean = ');  % ÇÀÄÀÍÍÎÅ ÌÀÒÅÌÀÒÈ×ÅÑÊÎÅ ÎÆÈÄÀÍÈÅ ØÓÌÀ
Var = input('Var = ');    % ÇÀÄÀÍÍÀß ÄÈÑÏÅĞÑÈß ØÓÌÀ 
disp('% Ïğîâåğüòå ÏĞÀÂÈËÜÍÎÑÒÜ ââîäà ÈÑÕÎÄÍÛÕ ÄÀÍÍÛÕ')
disp('% Ïğè ÏĞÀÂÈËÜÍÛÕ ÈÑÕÎÄÍÛÕ ÄÀÍÍÛÕ ââåäèòå 1')
disp('% Ïğè ÍÅÏĞÀÂÈËÜÍÛÕ ÈÑÕÎÄÍÛÕ ÄÀÍÍÛÕ ââåäèòå 0 è ÏÎÂÒÎĞÈÒÅ ââîä')
DATA = input('--> '); 
end
disp('%')
disp('%')
disp('% Äëÿ ïğîäîëæåíèÿ íàæìèòå <ENTER>')
pause 
disp('%')
disp('%')
disp('% ï.1. ÖÈÔĞÎÂÎÉ ÅÄÈÍÈ×ÍÛÉ ÈÌÏÓËÜÑ')
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÃĞÀÔÈÊÎÂ öèôğîâîãî åäèíè÷íîãî èìïóëüñà íàæìèòå <ENTER>')
pause 
n = 0:(N-1); nT = T.*n;      % ÄÈÑÊĞÅÒÍÎÅ ÍÎĞÌÈĞÎÂÀÍÍÎÅ È ÍÅÍÎĞÌÈĞÎÂÀÍÍÎÅ ÂĞÅÌß
u0 = [1 zeros(1,(N-1))];     % ÖÈÔĞÎÂÎÉ ÅÄÈÍÈ×ÍÛÉ ÈÌÏÓËÜÑ
figure('Name','Digital Unit Impulse, Unit Step, and Discrete Exponent','NumberTitle', 'off')
subplot(3,2,1),stem(nT,u0,'Linewidth',2), grid
title('Digital Unit Impulse u0(nT)')
subplot(3,2,2),stem(n,u0,'Linewidth',2), grid 
title('Digital Unit Impulse u0(n)')
disp('%')
disp('%')
disp('% Äëÿ ïğîäîëæåíèÿ íàæìèòå <ENTER>')
pause 
disp('%')
disp('%')
disp('% ï.2. ÖÈÔĞÎÂÎÉ ÅÄÈÍÈ×ÍÛÉ ÑÊÀ×ÎÊ');
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÃĞÀÔÈÊÎÂ öèôğîâîãî åäèíè÷íîãî ñêà÷êà íàæìèòå <ENTER>')
pause 
u1 = [1 ones(1,(N-1))];       % ÖÈÔĞÎÂÎÉ ÅÄÈÍÈ×ÍÛÉ ÑÊÀ×ÎÊ
subplot(3,2,3),stem(nT,u1,'Linewidth',2), grid
title('Digital Unit Step u1(nT)'), 
subplot(3,2,4),stem(n,u1,'Linewidth',2), grid
title('Digital Unit Step u1(n)')
disp('%')
disp('%')
disp('% Äëÿ ïğîäîëæåíèÿ íàæìèòå <ENTER>')
pause 
disp('%')
disp('%')
disp('% ï.3. ÄÈÑÊĞÅÒÍÀß İÊÑÏÎÍÅÍÒÀ')
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÃĞÀÔÈÊÎÂ äèñêğåòíîé ıêñïîíåíòû íàæìèòå <ENTER>')
pause
x1 = a.^n;                   % ÄÈÑÊĞÅÒÍÀß İÊÑÏÎÍÅÍÒÀ
subplot(3,2,5),stem(nT,x1,'Linewidth',2), xlabel('nT'), grid
title('Discrete Exponent x1(nT)')
subplot(3,2,6),stem(n, x1,'Linewidth',2), xlabel('n'), grid
title('Discrete Exponent x1(n)'),
disp('%')
disp('%')
disp('% Äëÿ ïğîäîëæåíèÿ íàæìèòå <ENTER>')
pause 
disp('%')
disp('%')
disp('% ï.4. ÄÈÑÊĞÅÒÍÛÉ ÊÎÌÏËÅÊÑÍÛÉ ÃÀĞÌÎÍÈ×ÅÑÊÈÉ ÑÈÃÍÀË')
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÃĞÀÔÈÊÎÂ âåùåñòâåííîé è ìíèìîé ÷àñòåé')
disp('% ãàğìîíè÷åñêîãî ñèãíàëà íàæìèòå <ENTER>')
pause 
x2 = C.*exp(j*w0.*n);  % ÄÈÑÊĞÅÒÍÛÉ ÊÎÌÏËÅÊÑÍÛÉ ÃÀĞÌÎÍÈ×ÅÑÊÈÉ ÑÈÃÍÀË
figure('Name','Discrete Harmonic Signal','NumberTitle', 'off')
subplot(2,1,1),stem(n,real(x2) ,'Linewidth',2), grid
title('Discrete Harmonic Signal: REAL [x2(n)]')
subplot(2,1,2),stem(n,imag(x2) ,'Linewidth',2), xlabel('n'), grid
title(' Discrete Harmonic Signal: IMAG [x2(n)]')
disp('%')
disp('%')
disp('% Äëÿ ïğîäîëæåíèÿ íàæìèòå <ENTER>')
pause 
disp('%')
disp('%')
disp('% ï.5. ÇÀÄÅĞÆÀÍÍÛÅ ÏÎÑËÅÄÎÂÀÒÅËÜÍÎÑÒÈ')
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÃĞÀÔÈÊÎÂ çàäåğæàííûõ ïîñëåäîâàòåëüíîñòåé íàæìèòå <ENTER>')
pause
u0_m = [zeros(1,m) u0(1:(N-m))];    % ÇÀÄÅĞÆÀÍÍÛÉ ÖÈÔĞÎÂÎÉ ÅÄÈÍÈ×ÍÛÉ ÈÌÏÓËÜÑ
u1_m = [zeros(1,m) u1(1:(N-m))];    % ÇÀÄÅĞÆÀÍÍÛÉ ÖÈÔĞÎÂÎÉ ÅÄÈÍÈ×ÍÛÉ ÑÊÀ×ÎÊ
x1_m = [zeros(1,m) x1(1:(N-m))];    % ÇÀÄÅĞÆÀÍÍÀß ÄÈÑÊĞÅÒÍÀß İÊÑÏÎÍÅÍÒÀ
figure('Name','Delayed Discrete Signals','NumberTitle', 'off')
subplot(3,1,1),stem(n,u0_m,'Linewidth',2), grid
title ('Delayed Digital Unit Impulse u0(n-m)')
subplot(3,1,2),stem(n,u1_m,'Linewidth',2), grid
title ('Delayed Digital Unit Step u1(n-m)')
subplot(3,1,3),stem(n,x1_m,'Linewidth',2),xlabel('n'), grid
title ('Delayed Discrete Exponent x1(n-m)')
disp('%')
disp('%')
disp('% Äëÿ ïğîäîëæåíèÿ íàæìèòå <ENTER>')
pause 
disp('%')
disp('%')
disp('% ï.6. ÄÈÑÊĞÅÒÍÛÉ ÏĞßÌÎÓÃÎËÜÍÛÉ ÈÌÏÓËÜÑ')
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÃĞÀÔÈÊÎÂ äèñêğåòíîãî ïğÿìîóãîëüíîãî èìïóëüñà íàæìèòå <ENTER>')
pause
x3_1 = U*rectpuls(n-n0,2*n_imp); x3_1(1:n0) = 0; % ÔÎĞÌÈĞÎÂÀÍÈÅ ÈÌÏÓËÜÑÀ Ñ ÏÎÌÎÙÜŞ ÔÓÍÊÖÈÈ rectpuls 
x3_2 = [zeros(1,n0) U.*u1((n0+1):(n0+n_imp))...
zeros(1,N-(n0+n_imp))];     % ÔÎĞÌÈĞÎÂÀÍÈÅ ÈÌÏÓËÜÑÀ Ñ ÏÎÌÎÙÜŞ ÖÈÔĞÎÂÎÃÎ ÅÄÈÍÈ×ÍÎÃÎ ÑÊÀ×ÊÀ
figure('Name','Discrete Rectangular and Triangular Impulses','NumberTitle', 'off')
subplot(3,1,1),stem(n,x3_1,'Linewidth',2), grid
title('Discrete Rectangular Impulse x3 1(n)')
subplot(3,1,2),stem(n,x3_2,'Linewidth',2), grid
title('Discrete Rectangular Impulse x3 2 (n)')
disp('%')
disp('%')
disp('% Äëÿ ïğîäîëæåíèÿ íàæìèòå <ENTER>')
pause 
disp('%')
disp('%')
disp('% ï.7. ÄÈÑÊĞÅÒÍÛÉ ÒĞÅÓÃÎËÜÍÛÉ ÈÌÏÓËÜÑ')
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÃĞÀÔÈÊÀ äèñêğåòíîãî òğåóãîëüíîãî èìïóëüñà íàæìèòå <ENTER>')
pause
x4 = conv(x3_1,x3_1);           % ÄÈÑÊĞÅÒÍÛÉ ÒĞÅÓÃÎËÜÍÛÉ ÈÌÏÓËÜÑ
L = 2*N-1;                      % ÄËÈÍÀ ÑÂÅĞÒÊÈ
n = 0:(L-1);                    % ÄÈÑÊĞÅÒÍÎÅ ÍÎĞÌÈĞÎÂÀÍÍÎÅ ÂĞÅÌß
subplot(3,1,3),stem(n,x4,'Linewidth',2), xlabel('n'), grid
title('Discrete Triangular Impulse x4(n)')
disp('%')
disp('%')
disp('% Äëÿ ïğîäîëæåíèÿ íàæìèòå <ENTER>')
pause 
disp('%')
disp('%')
disp('% ï.8. ËÈÍÅÉÍÀß ÊÎÌÁÈÍÀÖÈß ÄÈÑÊĞÅÒÍÛÕ ÃÀĞÌÎÍÈ×ÅÑÊÈÕ ÑÈÃÍÀËÎÂ')
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÃĞÀÔÈÊÎÂ ãàğìîíè÷åñêèõ ñèãíàëîâ è èõ ëèíåéíîé êîìáèíàöèè íàæìèòå <ENTER>')
pause
n = 0:(5*N-1);                         % ÄÈÑÊĞÅÒÍÎÅ ÍÎĞÌÈĞÎÂÀÍÍÎÅ ÂĞÅÌß
xi = repmat(B,length(n),1).*sin(n'*w); % ÌÀÒĞÈÖÀ ÄÈÑÊĞÅÒÍÛÕ ÃÀĞÌÎÍÈÊ
ai = repmat(A,length(n),1);            % ÌÀÒĞÈÖÀ ÊÎİÔÔÈÖÈÅÍÒÎÂ
x5 = sum((ai.* xi)');         % ËÈÍÅÉÍÀß ÊÎÌÁÈÍÀÖÈß ÄÈÑÊĞÅÒÍÛÕ ÃÀĞÌÎÍÈÊ
figure('Name','Discrete Harmonic Signals and their Linear Combination','NumberTitle', 'off')
subplot(4,1,1),stem(n, xi(:,1),'Linewidth',2), grid
title('First Discrete Harmonic Signal')
subplot(4,1,2),stem(n, xi(:,2),'Linewidth',2), grid
title('Second Discrete Harmonic Signal')
subplot(4,1,3),stem(n, xi(:,3),'Linewidth',2), grid
title('Third Discrete Harmonic Signal')
subplot(4,1,4),stem(n,x5,'Linewidth',2), xlabel('n'), grid
title('Linear Combination x5(n)') 
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÑĞÅÄÍÅÃÎ ÇÍÀ×ÅÍÈß, İÍÅĞÃÈÈ è ÑĞÅÄÍÅÉ ÌÎÙÍÎÑÒÈ ñèãíàëà x5 íàæìèòå <ENTER>')
pause
mean_x5 = mean(x5);               % ÑĞÅÄÍÅÅ ÇÍÀ×ÅÍÈÅ ÑÈÃÍÀËÀ
E = sum(x5.^2);                   % İÍÅĞÃÈß ÑÈÃÍÀËÀ
P = sum(x5.^2)/length(x5);        % ÑĞÅÄÍßß ÌÎÙÍÎÑÒÜ ÑÈÃÍÀËÀ
disp('%')
disp('%')
disp(['  mean_x5 = ' num2str(mean_x5) '  E = ' num2str(E) '  P = ' num2str(P)])
disp('%')
disp('%')
disp('% Äëÿ ïğîäîëæåíèÿ íàæìèòå <ENTER>')
pause 
disp('%')
disp('%')
disp('% ï.9. ÄÈÑÊĞÅÒÍÛÉ ÃÀĞÌÎÍÈ×ÅÑÊÈÉ ÑÈÃÍÀË Ñ İÊÑÏÎÍÅÍÖÈÀËÜÍÎÉ ÎÃÈÁÀŞÙÅÉ')
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÃĞÀÔÈÊÀ ãàğìîíè÷åñêîãî ñèãíàëà ñ ıêñïîíåíöèàëüíîé îãèáàşùåé íàæìèòå <ENTER>')
pause 
n = 0:(N-1);                       % ÄÈÑÊĞÅÒÍÎÅ ÍÎĞÌÈĞÎÂÀÍÍÎÅ ÂĞÅÌß
x = C.*sin(w0.*n);                 % ÄÈÑÊĞÅÒÍÛÉ ÃÀĞÌÎÍÈ×ÅÑÊÈÉ ÑÈÃÍÀË
x6 = x.*(abs(a).^n);               % ÄÈÑÊĞÅÒÍÛÉ ÃÀĞÌÎÍÈ×ÅÑÊÈÉ ÑÈÃÍÀË Ñ İÊÑÏÎÍÅÍÖÈÀËÜÍÎÉ ÎÃÈÁÀŞÙÅÉ
figure('Name','Harmonic Signal with Exponential Envelope.  Periodic Sequence of Rectangular Impulses','NumberTitle', 'off')
subplot(2,1,1),stem(n,x6,'Linewidth',2), grid
title('Harmonic Signal with Exponential Envelope x6(n)')
disp('%')
disp('%')
disp('% Äëÿ ïğîäîëæåíèÿ íàæìèòå <ENTER>')
pause 
disp('%')
disp('%')
disp('% ï.10. ÏÅĞÈÎÄÈ×ÅÑÊÀß ÏÎÑËÅÄÎÂÀÒÅËÜÍÎÑÒÜ ÄÈÑÊĞÅÒÍÛÕ ÏĞßÌÎÓÃÎËÜÍÛÕ ÈÌÏÓËÜÑÎÂ')
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÃĞÀÔÈÊÀ ïÿòè ïåğèîäîâ ïîñëåäîâàòåëüíîñòè íàæìèòå <ENTER>')
pause
xp = [U.*u1(1:n_imp) zeros(1,n_imp)];    % ÏÅĞÈÎÄ ÏÎÑËÅÄÎÂÀÒÅËÜÍÎÑÒÈ
p = 5;                                   % ×ÈÑËÎ ÏÅĞÈÎÄÎÂ 
x7 =  repmat(xp,1,p);             % ÏÅĞÈÎÄÈ×ÅÑÊÀß ÏÎÑËÅÄÎÂÀÒÅËÜÍÎÑÒÜ
n = 0:(length(x7)-1);             % ÄÈÑÊĞÅÒÍÎÅ ÍÎĞÌÈĞÎÂÀÍÍÎÅ ÂĞÅÌß
subplot(2,1,2), stem(n,x7,'Linewidth',2), xlabel('n'), grid
title('Periodic Sequence of Rectangular Impulses x7(n)') 
disp('%')
disp('%')
disp('% Äëÿ ïğîäîëæåíèÿ íàæìèòå <ENTER>')
pause 
disp('%')
disp('%')
disp('% ï.11. ĞÀÂÍÎÌÅĞÍÛÉ ÁÅËÛÉ ØÓÌ')
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÎÖÅÍÎÊ ÌÀÒÅÌÀÒÈ×ÅÑÊÎÃÎ ÎÆÈÄÀÍÈß è ÄÈÑÏÅĞÑÈÈ ØÓÌÀ íàæìèòå <ENTER>')
pause
r_uniform = rand(1,10000);           % ĞÀÂÍÎÌÅĞÍÛÉ ÁÅËÛÉ ØÓÌ
mean_uniform = mean(r_uniform);      % ÎÖÅÍÊÀ ÌÀÒ. ÎÆÈÄÀÍÈß ØÓÌÀ
var_uniform = var(r_uniform);        % ÎÖÅÍÊÀ ÄÈÑÏÅĞÑÈÈ ØÓÌÀ
disp('%')
disp('%')
disp(['  mean_uniform = ' num2str(mean_uniform) '  var_uniform = ' num2str(var_uniform)]) 
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ãğàôèêà ÀÂÒÎÊÎÂÀĞÈÀÖÈÎÍÍÎÉ ÔÓÍÊÖÈÈ íàæìèòå <ENTER>')
pause
r_r_uniform = (1/length(r_uniform)).*xcov(r_uniform);   % ÎÖÅÍÊÀ ÀÂÒÎÊÎÂÀĞÈÀÖÈÎÍÍÎÉ ÔÓÍÊÖÈÈ ĞÀÂÍÎÌÅĞÍÎÃÎ ÁÅËÎÃÎ ØÓÌÀ
m = -(length(r_uniform)-1):(length(r_uniform)-1);       % ÂÅÊÒÎĞ ÄÈÑÊĞÅÒÍÛÕ ÑÄÂÈÃÎÂ ÄËß ÀÂÒÎÊÎÂÀĞÈÀÖÈÎÍÍÎÉ ÔÓÍÊÖÈÈ 
figure('Name','Autocovariance Function of Uniform White Noise','NumberTitle', 'off')
stem(m,r_r_uniform,'Linewidth',2), xlabel('m'), grid
title('Autocovariance Function of Uniform White Noise')
disp('%')
disp('%')
disp('% Äëÿ ïğîäîëæåíèÿ íàæìèòå <ENTER>')
pause 
disp('%')
disp('%')
disp('% ï.12. ÍÎĞÌÀËÜÍÛÉ ÁÅËÛÉ ØÓÌ')
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÎÖÅÍÎÊ ÌÀÒÅÌÀÒÈ×ÅÑÊÎÃÎ ÎÆÈÄÀÍÈß è ÄÈÑÏÅĞÑÈÈ øóìà íàæìèòå <ENTER>')
pause
r_norm = randn(1,10000);           % ÍÎĞÌÀËÜÍÛÉ ÁÅËÛÉ ØÓÌ
mean_norm = mean(r_norm);          % ÎÖÅÍÊÀ ÌÀÒ. ÎÆÈÄÀÍÈß ØÓÌÀ
var_norm = var(r_norm);            % ÎÖÅÍÊÀ ÄÈÑÏÅĞÑÈÈ ØÓÌÀ
disp('%')
disp('%')
disp(['  mean_norm = ' num2str(mean_norm) '  var_norm = ' num2str(var_norm)]) 
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ãğàôèêà ÀÊÔ íàæìèòå <ENTER>')
pause
R_r_norm = (1/length(r_norm)).*xcorr(r_norm);   % ÎÖÅÍÊÀ ÀÊÔ ÍÎĞÌÀËÜÍÎÃÎ ÁÅËÎÃÎ ØÓÌÀ 
m = -(length(r_norm)-1):(length(r_norm)-1);     % ÂÅÊÒÎĞ ÄÈÑÊĞÅÒÍÛÕ ÑÄÂÈÃÎÂ ÄËß ÀÊÔ 
figure('Name','ACF of White Gaussian Noise','NumberTitle', 'off')
stem(m,R_r_norm,'Linewidth',2), xlabel('m'), grid
title('ACF of White Gaussian Noise')
disp('%')
disp('%')
disp('% Äëÿ ïğîäîëæåíèÿ íàæìèòå <ENTER>')
pause
disp('%')
disp('%')
disp('% ï.13. ÀÄÄÈÒÈÂÍÀß ÑÌÅÑÜ ÄÈÑÊĞÅÒÍÎÃÎ ÃÀĞÌÎÍÈ×ÅÑÊÎÃÎ ÑÈÃÍÀËÀ Ñ ÍÎĞÌÀËÜÍÛÌ ÁÅËÛÌ ØÓÌÎÌ')
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÃĞÀÔÈÊÀ àääèòèâíîé ñìåñè ñèãíàëà ñ øóìîì íàæìèòå <ENTER>')
pause
n = 0:(N-1);                     % ÄÈÑÊĞÅÒÍÎÅ ÍÎĞÌÈĞÎÂÀÍÍÎÅ ÂĞÅÌß
x8 = x+randn(1,N);               % ÀÄÄÈÒÈÂÍÀß ÑÌÅÑÜ ÑÈÃÍÀËÀ Ñ ØÓÌÎÌ
figure('Name','Mixture of Harmonic Signal and White Gaussian Noise and ACF','NumberTitle', 'off')
subplot(2,1,1),stem(n,x8,'Linewidth',2),xlabel('n'), grid
title('Mixture of Harmonic Signal and White Gaussian Noise x8(n)')
disp('%')
disp('%')
disp('% Äëÿ ïğîäîëæåíèÿ íàæìèòå <ENTER>')
pause 
disp('%')
disp('%')
disp('% ï.14. ÀÊÔ ÀÄÄÈÒÈÂÍÎÉ ÑÌÅÑÈ ÄÈÑÊĞÅÒÍÎÃÎ ÃÀĞÌÎÍÈ×ÅÑÊÎÃÎ ÑÈÃÍÀËÀ Ñ ÍÎĞÌÀËÜÍÛÌ ÁÅËÛÌ ØÓÌÎÌ')
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÃĞÀÔÈÊÀ ÀÊÔ íàæìèòå <ENTER>')
pause 
R = (1/N).*xcorr(x8);            % ÎÖÅÍÊÀ ÀÊÔ 
m = -(N-1):(N-1);                % ÂÅÊÒÎĞ ÄÈÑÊĞÅÒÍÛÕ ÑÄÂÈÃÎÂ ÄËß ÀÊÔ 
subplot(2,1,2),stem((m),R,'Linewidth',2),xlabel('m'), grid
title('ACF R(m)')
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÄÈÑÏÅĞÑÈÈ àääèòèâíîé ñìåñè ñèãíàëà ñ øóìîì è ÀÊÔ R(N) íàæìèòå <ENTER>') 
pause 
disp('%')
disp('%')
disp(['  var_x8 = ' num2str(var(x8))])
disp(['  R(N) = ' num2str(R(N))])
disp('%')
disp('%')
disp('% Äëÿ ïğîäîëæåíèÿ íàæìèòå <ENTER>')
pause 
disp('%')
disp('%')
disp('% ï.15. ÍÎĞÌÀËÜÍÛÉ ÁÅËÛÉ ØÓÌ Ñ ÇÀÄÀÍÍÛÌÈ ÑÒÀÒÈÑÒÈ×ÅÑÊÈÌÈ ÕÀĞÀÊÒÅĞÈÑÒÈÊÀÌÈ')
r_normMean = randn(1,10000)+Mean;      % ÍÎĞÌÀËÜÍÛÉ ÁÅËÛÉ ØÓÌ Ñ ÇÀÄÀÍÍÛÌ ÌÀÒÅÌÀÒÈ×ÅÑÊÈÌ ÎÆÈÄÀÍÈÅÌ
r_normVar = sqrt(Var).*randn(1,10000); % ÍÎĞÌÀËÜÍÛÉ ÁÅËÛÉ ØÓÌ Ñ ÇÀÄÀÍÍÎÉ ÄÈÑÏÅĞÑÈÅÉ
r_normMeanVar = sqrt(Var).*randn(1,10000)+ Mean; % ÍÎĞÌÀËÜÍÛÉ ÁÅËÛÉ ØÓÌ Ñ ÇÀÄÀÍÍÛÌÈ ÌÀÒÅÌÀÒÈ×ÅÑÊÈÌ ÎÆÈÄÀÍÈÅÌ È ÄÈÑÏÅĞÑÈÅÉ
MAX = max([r_norm r_normMean r_normVar r_normMeanVar]); 
% ÌÀÊÑÈÌÀËÜÍÎÅ ÇÍÀ×ÅÍÈÅ ØÓÌÀ ÑĞÅÄÈ ×ÅÒÛĞÅÕ ÅÃÎ ĞÀÇÍÎÂÈÄÍÎÑÒÅÉ
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÃĞÀÔÈÊÎÂ íîğìàëüíîãî áåëîãî øóìà íàæìèòå <ENTER>')
pause
figure('Name','White Gaussian Noises with different statistics','NumberTitle', 'off')
subplot(4,1,1), plot(r_norm), grid, ylim([-MAX MAX])
title(strcat([' Mean value = ',num2str(mean(r_norm)),'   Variance = ',num2str(var(r_norm))]))
subplot(4,1,2), plot(r_normMean), grid, ylim([-MAX MAX])
title(strcat([' Mean value = ',num2str(mean(r_normMean)),'   Variance = ',num2str(var(r_normMean))]))
subplot(4,1,3), plot(r_normVar), grid, ylim([-MAX MAX])
title(strcat([' Mean value = ',num2str(mean(r_normVar)),'   Variance = ',num2str(var(r_normVar))]))
subplot(4,1,4), plot(r_normMeanVar), xlabel('n'), grid, ylim([-MAX MAX])
title(strcat([' Mean value = ',num2str(mean(r_normMeanVar)),'   Variance = ',num2str(var(r_normMeanVar))]))
disp('%')
disp('%')
disp('% Äëÿ âûâîäà ÃÈÑÒÎÃĞÀÌÌ íîğìàëüíîãî áåëîãî øóìà íàæìèòå <ENTER>')
pause
figure('Name','Histograms with different statistics','NumberTitle', 'off')
subplot(4,1,1), hist(r_norm), grid, xlim([-MAX MAX]) 
title(strcat([' Mean value = ',num2str(mean(r_norm)),'   Variance = ',num2str(var(r_norm))]))
subplot(4,1,2), hist(r_normMean), grid, xlim([-MAX MAX])
title(strcat([' Mean value =  ',num2str(mean(r_normMean)),'   Variance = ',num2str(var(r_normMean))]))
subplot(4,1,3), hist(r_normVar), grid, xlim([-MAX MAX])
title(strcat([' Mean value = ',num2str(mean(r_normVar)),'   Variance = ',num2str(var(r_normVar))]))
subplot(4,1,4),hist(r_normMeanVar), grid, xlim([-MAX MAX])
title(strcat([' Mean value = ',num2str(mean(r_normMeanVar)),'   Variance = ',num2str(var(r_normMeanVar))]))
disp('%')
disp('%')
disp('% ĞÀÁÎÒÀ ÇÀÂÅĞØÅÍÀ')















































































