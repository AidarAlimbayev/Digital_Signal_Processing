script
clc
disp('% кп ╧20. лндекхпнбюмхе яхярелш ндмнйпюрмни дежхлюжхх я онкхтюгмни ярпсйрспни')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше');
DATA=0;
while DATA==0
Nb = input('Nb = ');  % мнлеп апхцюдш
A1 = input('A1 = ');  % юлокхрсдш дхяйпермшу цюплнмхй бундмнцн яхцмюкю
A2 = input('A2 = ');
f1 = input('f1 = ');  % вюярнрш (цЖ) дхяйпермшу цюплнмхй бундмнцн яхцмюкю
f2 = input('f2 = ');
Nd = input('Nd = ');        % оепхнд бундмнцн яхцмюкю
Fs_d = input('Fs_d = ');    % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
M = input('M = ');          % йнщттхжхемр дежхлюжхх
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
disp('%')
disp('% дКЪ БШБНДЮ ГМЮВЕМХИ Nd, Fs_d Х M МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp(['	N = Nd = ',num2str(Nd),'      Fs = Fs_d = ',num2str(Fs_d)])
disp('%')
disp(['      	M = ',num2str(M)])
N = Nd;           % оепхнд бундмнцн яхцмюкю
Fs = Fs_d;        % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
n = 0:(N-1);      % дхяйпермне мнплхпнбюммне бпелъ бундмнцн яхцмюкю
x = A1*sin((2*pi*f1/Fs).*n)+A2*sin((2*pi*f2/Fs).*n);% бундмни яхцмюк
y = filter(Hd,x);               % бшундмни яхцмюк
nd_start = ceil(length(y)/2+1); % мювюкн сярюмнбхбьецняъ пефхлю
yd = y(nd_start:end);           % бшундмни яхцмюк б сярюмнбхбьеляъ пефхле
nd = nd_start:(nd_start+length(yd)-1);  % дхяйпермне мнплхпнбюммне бпелъ дкъ сярюмнбхбьецняъ пефхлю
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб бундмнцн х бшундмнцн ЯХЦМЮКНБ МЮФЛХРЕ <ENTER>')
pause
figure('Name','Input and Output Signals in Polyphase Structure Decimation System','NumberTitle', 'off')
subplot(3,1,1), stem(n,x), grid, xlabel('n')
title(strcat(['Input Signal x(n)      N = ',num2str(N)]))
subplot(3,1,2), stem(0:length(y)-1,y), grid
title(strcat(['M = ',num2str(M),'  Output Signal y(n)   length(y) = ',num2str(length(y))]))
subplot(3,1,3), stem(nd,yd), grid, xlabel('n')
title(strcat(['Output Signal y(n)   Starting point n = ',num2str(nd_start),'   length(y) = ',num2str(length(yd))]))
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юлокхрсдмшу яоейрпнб Х юву йху-ТХКЭРПЮ МЮФЛХРЕ <ENTER>')
pause
X = fft(x);                % дот бундмнцн яхцмюкю
MODX = (2/N)*abs(X);       % юлокхрсдмши яоейрп бундмнцн яхцмюкю
MODX(1) = (1/N)*abs(X(1));
Yd = fft(yd);              % дот бшундмнцн яхцмюкю б сярюмнбхбьеляъ пефхле
MODYd = (2/length(Yd))*abs(Yd);    % юлокхрсдмши яоейрп бшундмнцн яхцмюкю
MODYd(1) = (1/length(Yd))*abs(Yd(1));
k = 0:(N-1);           % дхяйпермше мнплхпнбюммше вюярнрш бундмнцн яхцмюкю
kd = 0:(length(Yd)-1); % дхяйпермше мнплхпнбюммше вюярнрш бшундмнцн яхцмюкю
f = 0:Fs/1000:Fs;                   % яерйю вюярнр дкъ юву йху-тхкэрпю
hd = Hd.Numerator;                  % ху йху-тхкэрпю
MAG = abs(freqz(hd,1,f,Fs));        % юву йху-тхкэрпю
figure('Name','Amplitude Spectrums and Magnitude Response in Polyphase Structure Decimation System','NumberTitle', 'off')
subplot(2,1,1), stem(k*(Fs/N),MODX,'MarkerSize',3,'Linewidth',2)
grid, xlabel('f (Hz)'), xlim([0 Fs])
title(strcat('Amplitude Spectrum x(n) and FIR Magnitude Response'))
hold on, plot(f,MAG,'r','Linewidth',2)
subplot(2,1,2),stem(kd*((Fs/M)/length(Yd)),MODYd,'MarkerSize',3, 'Linewidth',2)
grid, xlabel('f (Hz)'), xlim([0 Fs])
title(strcat(['M = ',num2str(M),'   Amplitude Spectrum y(n)']))
disp('%')
disp('%')
disp('% лндекхпнбюмхе яхярелш ндмнйпюрмни дежхлюжхх я онкхтюгмни ярпсйрспни гюбепьемн')


























































































