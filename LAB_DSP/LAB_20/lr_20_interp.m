script
clc
disp('% кп ╧20. лндекхпнбюмхе яхярелш ндмнйпюрмни хмрепонкъжхх я онкхтюгмни ярпсйрспни')
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
Ni = input('Ni = ');        % оепхнд бундмнцн яхцмюкю
Fs_i = input('Fs_i = ');    % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
L = input('L = ');          % йнщттхжхемр хмрепонкъжхх
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
disp('%')
disp('% дКЪ БШБНДЮ ГМЮВЕМХИ Ni, Fs_i Х L МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp(['	N = Ni = ',num2str(Ni),'      Fs = Fs_i = ',num2str(Fs_i)])
disp('%')
disp(['      	L = ',num2str(L)])
N = Ni;                  % оепхнд бундмнцн яхцмюкю
Fs = Fs_i;               % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
n = 0:(N-1);             % дхяйпермне мнплхпнбюммне бпелъ бундмнцн яхцмюкю
x = A1*sin((2*pi*f1/Fs).*n)+A2*sin((2*pi*f2/Fs).*n); % бундмни яхцмюк
y = L*filter(Hi,x);      % бшундмни яхцмюк, слмнфеммши мю йнщттхжхемр хмрепонкъжхх
ni_start = ceil(length(y)/2+1);  % мювюкн сярюмнбхбьецняъ пефхлю
yi = y(ni_start:end);            % бшундмни яхцмюк б сярюмнбхбьеляъ пефхле
ni = ni_start:(ni_start+length(yi)-1);  % дхяйпермне мнплхпнбюммне бпелъ дкъ сярюмнбхбьецняъ пефхлю
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб бундмнцн х бшундмнцн ЯХЦМЮКНБ МЮФЛХРЕ <ENTER>')
pause
figure('Name','Input and Output Signals in Polyphase Structure Interpolation System','NumberTitle', 'off')
subplot(3,1,1), stem(n,x), grid, xlabel('n')
title(strcat(['Input Signal x(n)      N = ',num2str(N)]))
subplot(3,1,2), stem(0:length(y)-1,y), grid
title(strcat(['L = ',num2str(L),'   Output Signal y(n)   length(y) = ',num2str(length(y))]))
subplot(3,1,3), stem(ni,yi), grid, xlabel('n')
title(strcat(['Output Signal y(n)   Starting point n = ',num2str(ni_start),'   length(y) = ',num2str(length(yi))]))
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юлокхрсдмшу яоейрпнб Х юву йху-ТХКЭРПЮ МЮФЛХРЕ <ENTER>')
pause
X = fft(x);                      % дот бундмнцн яхцмюкю
MODX = (2/N)*abs(X);             % юлокхрсдмши яоейрп бундмнцн яхцмюкю
MODX(1) = (1/N)*abs(X(1));
Yi = fft(yi);                    % дот бшундмнцн яхцмюкю
MODYi = (2/length(Yi))*abs(Yi);  % юлокхрсдмши яоейрп бшундмнцн яхцмюкю
MODYi(1) = (1/length(Yi))*abs(Yi(1));
k = 0:(N-1);           % дхяйпермше мнплхпнбюммше вюярнрш бундмнцн яхцмюкю
ki = 0:(length(Yi)-1); % дхяйпермше мнплхпнбюммше вюярнрш бшундмнцн яхцмюкю
f = 0:(Fs*L)/1000:Fs*L;               % яерйю вюярнр дкъ юву йху-тхкэрпю
h = Hi.Numerator;                     % ху йху-тхкэрпю
MAG = abs(freqz(h,1,f,Fs*L));         % юву йху-тхкэрпю
figure('Name','Amplitude Spectrums and Magnitude Response in Polyphase Structure Interpolation System','NumberTitle', 'off')
subplot(2,1,1), stem(k*(Fs/N),MODX,'MarkerSize',3,'Linewidth',2)
grid, xlabel('f (Hz)'), xlim([0 L*Fs])
title(strcat(['Amplitude Spectrum x(n)   N = ',num2str(N)]))
subplot(2,1,2), stem(ki*(Fs*L/length(Yi)),MODYi,'MarkerSize',3, 'Linewidth',2)
grid, xlabel('f (Hz)'), xlim([0 L*Fs])
title(strcat(['L = ',num2str(L),'   Amplitude Spectrum y(n) and FIR Magnitude Response']))
hold on, plot(f,MAG,'r','Linewidth',2)
disp('%')
disp('%')
disp('% лндекхпнбюмхе яхярелш ндмнйпюрмни хмрепонкъжхх я онкхтюгмни ярпсйрспни гюбепьемн')































































































