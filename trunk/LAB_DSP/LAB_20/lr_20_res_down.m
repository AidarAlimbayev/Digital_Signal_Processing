script
clc
disp('% кп ╧20. лндекхпнбюмхе яхярелш ндмнйпюрмни оепедхяйперхгюжхх я онкхтюгмни ярпсйрспни')
disp('% опх онмхфемхх вюярнрш дхяйперхгюжхх')
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
Nrd = input('Nrd = ');        % оепхнд бундмнцн яхцмюкю
Fs_rd = input('Fs_rd = ');    % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
L = input('L = ');            % йнщттхжхемр хмрепонкъжхх
M = input('M = ');            % йнщттхжхемр дежхлюжхх
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ГМЮВЕМХИ Nrd, Fs_rd Х L/M МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp(['      N = Nrd = ',num2str(Nrd),'      Fs = Fs_rd = ',num2str(Fs_rd)])
disp('%')
disp(['      	L/M = ',num2str(L),'/',num2str(M)])
pause
N = Nrd;           % оепхнд бундмнцн яхцмюкю
Fs = Fs_rd;        % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
n = 0:(N-1);       % дхяйпермне мнплхпнбюммне бпелъ бундмнцн яхцмюкю
x = A1*sin((2*pi*f1/Fs).*n)+A2*sin((2*pi*f2/Fs).*n);% бундмни яхцмюк
y = L*filter(Hrd,x);       % бшундмни яхцмюк, слмнфеммши мю йнщттхжхемр хмрепонкъжхх
nrd_start = ceil(length(y)/2+1); % мювюкн сярюмнбхбьецняъ пефхлю
yrd = y(nrd_start:end);          % бшундмни яхцмюк б сярюмнбхбьеляъ пефхле
nrd = nrd_start:(nrd_start+length(yrd)-1);  % дхяйпермне мнплхпнбюммне бпелъ дкъ сярюмнбхбьецняъ пефхлю
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб бундмнцн х бшундмнцн ЯХЦМЮКНБ МЮФЛХРЕ <ENTER>')
figure('Name','Input and Output Signals in Polyphase Structure Resampling System','NumberTitle', 'off')
subplot(3,1,1), stem(n,x), grid, xlabel('n')
title(strcat(['Input Signal x(n)      N = ',num2str(N)]))
subplot(3,1,2), stem(0:length(y)-1,y), grid
title(strcat(['L = ',num2str(L),' M = ',num2str(M),'  Output Signal y(n)   length(y) = ',num2str(length(y))]))
subplot(3,1,3), stem(nrd,yrd), grid, xlabel('n')
title(strcat(['Output Signal y(n)   Starting point n = ',num2str(nrd_start),'   length(y) = ',num2str(length(yrd))]))
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юлокхрсдмшу яоейрпнб Х юву йху-ТХКЭРПЮ МЮФЛХРЕ <ENTER>')
pause
X = fft(x);                % дот бундмнцн яхцмюкю
MODX = (2/N)*abs(X);       % юлокхрсдмши яоейрп бундмнцн яхцмюкю
MODX(1) = (1/N)*abs(X(1));
Yrd = fft(yrd);            % дот бшундмнцн яхцмюкю б сярюмнбхбьеляъ пефхле
MODYrd = (2/length(Yrd))*abs(Yrd);    % юлокхрсдмши яоейрп бшундмнцн яхцмюкю
MODYrd(1) = (1/length(Yrd))*abs(Yrd(1));
k = 0:(N-1);           % дхяйпермше мнплхпнбюммше вюярнрш бундмнцн яхцмюкю
krd = 0:(length(Yrd)-1); % дхяйпермше мнплхпнбюммше вюярнрш бшундмнцн яхцмюкю
f = 0:Fs/1000:Fs;                   % яерйю вюярнр дкъ юву йху-тхкэрпю
h = Hrd.Numerator;                  % ху йху-тхкэрпю
MAG = abs(freqz(h,1,f,Fs));         % юву йху-тхкэрпю
figure('Name','Amplitude Spectrums and Magnitude Response in Polyphase Structure Resampling System','NumberTitle', 'off')
subplot(2,1,1), stem(k*(Fs/N),MODX,'MarkerSize',3,'Linewidth',2)
grid, xlabel('f (Hz)'), xlim([0 Fs])
title(strcat('Amplitude Spectrum x(n) and FIR Magnitude Response'))
hold on, plot(f,MAG,'r','Linewidth',2)
subplot(2,1,2)
stem(krd*(Fs*(L/M)/length(Yrd)),MODYrd,'MarkerSize',3, 'Linewidth',2)
grid, xlabel('f (Hz)'), xlim([0 Fs])
title(strcat(['L/M = ',num2str(L),'/',num2str(M),'  Amplitude Spectrum y(n)']))
disp('%')
disp('%')
disp('% лндекхпнбюмхе яхярелш ндмнйпюрмни оепедхяйперхгюжхх я онкхтюгмни ярпсйрспни гюбепьемн')
































































































