script
clc
disp('% кп ╧20. лндекхпнбюмхе яхярелш ндмнйпюрмни оепедхяйперхгюжхх я онкхтюгмни ярпсйрспни')
disp('% опх онбшьемхх вюярнрш дхяйперхгюжхх')
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
Nri = input('Nri = ');        % оепхнд бундмнцн яхцмюкю
Fs_ri = input('Fs_ri = ');    % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
L = input('L = ');            % йнщттхжхемр хмрепонкъжхх
M = input('M = ');            % йнщттхжхемр дежхлюжхх
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ГМЮВЕМХИ Nri, Fs_ri Х L/M МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp(['      N = Nri = ',num2str(Nri),'      Fs = Fs_ri = ',num2str(Fs_ri)])
disp('%')
disp(['      	L/M = ',num2str(L),'/',num2str(M)])
N = Nri;           % оепхнд бундмнцн яхцмюкю
Fs = Fs_ri;        % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
n = 0:(N-1);       % дхяйпермне мнплхпнбюммне бпелъ бундмнцн яхцмюкю
x = A1*sin((2*pi*f1/Fs).*n)+A2*sin((2*pi*f2/Fs).*n); % бундмни яхцмюк
y = L*filter(Hri,x);             % бшундмни яхцмюк, слмнфеммши мю йнщттхжхемр хмрепонкъжхх
nri_start = ceil(length(y)/2+1); % мювюкн сярюмнбхбьецняъ пефхлю
yri = y(nri_start:end);          % бшундмни яхцмюк б сярюмнбхбьеляъ пефхле
nri = nri_start:(nri_start+length(yri)-1);  % дхяйпермне мнплхпнбюммне бпелъ дкъ сярюмнбхбьецняъ пефхлю
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб бундмнцн х бшундмнцн ЯХЦМЮКНБ МЮФЛХРЕ <ENTER>')
pause
figure('Name','Input and Output Signals in Polyphase Structure Resampling System','NumberTitle', 'off')
subplot(3,1,1), stem(n,x), grid, xlabel('n')
title(strcat(['Input Signal x(n)      N = ',num2str(N)]))
subplot(3,1,2), stem(0:length(y)-1,y), grid
title(strcat(['L = ',num2str(L),' M = ',num2str(M),'  Output Signal y(n)   length(y) = ',num2str(length(y))]))
subplot(3,1,3), stem(nri,yri), grid, xlabel('n')
title(strcat(['Output Signal y(n)   Starting point n = ',num2str(nri_start),'   length(y) = ',num2str(length(yri))]))
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юлокхрсдмшу яоейрпнб Х юву йху-ТХКЭРПЮ МЮФЛХРЕ <ENTER>')
pause
X = fft(x);               % дот бундмнцн яхцмюкю
MODX = (2/N)*abs(X);      % юлокхрсдмши яоейрп бундмнцн яхцмюкю
MODX(1) = (1/N)*abs(X(1));
Yri = fft(yri);           % дот бшундмнцн яхцмюкю б сярюмнбхбьеляъ пефхле
MODYri = (2/length(Yri))*abs(Yri);    % юлокхрсдмши яоейрп бшундмнцн яхцмюкю
MODYri(1) = (1/length(Yri))*abs(Yri(1));
k = 0:(N-1);           % дхяйпермше мнплхпнбюммше вюярнрш бундмнцн яхцмюкю
kri = 0:(length(Yri)-1); % дхяйпермше мнплхпнбюммше вюярнрш бшундмнцн яхцмюкю
f = 0:(Fs*(L/M))/1000:Fs*(L/M);        % яерйю вюярнр дкъ юву йху-тхкэрпю
h = Hri.Numerator;                     % ху йху-тхкэрпю
MAG = abs(freqz(h,1,f,Fs*(L/M)));      % юву йху-тхкэрпю
figure('Name','Amplitude Spectrums and Magnitude Response in Polyphase Structure Resampling System','NumberTitle', 'off')
subplot(2,1,1), stem(k*(Fs/N),MODX,'MarkerSize',3,'Linewidth',2)
grid, xlabel('f (Hz)'), xlim([0 (L/M)*Fs])
title(strcat(['Amplitude Spectrum x(n)   N = ',num2str(N)]))
subplot(2,1,2)
stem(kri*(Fs*(L/M)/length(Yri)),MODYri,'MarkerSize',3, 'Linewidth',2)
grid, xlabel('f (Hz)'), xlim([0 (L/M)*Fs])
title(strcat(['L/M = ',num2str(L),'/',num2str(M),'  Amplitude Spectrum y(n) and FIR Magnitude Response']))
hold on, plot(f,MAG,'r','Linewidth',2)
disp('%')
disp('%')
disp('% лндекхпнбюмхе яхярелш ндмнйпюрмни оепедхяйперхгюжхх я онкхтюгмни ярпсйрспни гюбепьемн')




























































































