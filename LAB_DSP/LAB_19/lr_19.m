script
clc
clear
disp('% кп ╧19. лмнцняйнпнярмше яхярелш жня')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше');
DATA=0;
while DATA==0
Nb = input('Nb = ');   % мнлеп апхцюдш
A1 = input('A1 = ');   % юлокхрсдш дхяйпермшу цюплнмхй бундмнцн яхцмюкю
A2 = input('A2 = ');
f1 = input('f1 = ');   % вюярнрш (цЖ) дхяйпермшу цюплнмхй бундмнцн яхцмюкю
f2 = input('f2 = ');
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA=input('--> ');
end
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.1. лндекхпнбюмхе яхярелш ндмнйпюрмни хмрепонкъжхх')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше');
DATA=0;
while DATA==0
Ni = input('Ni = ');        % оепхнд бундмнцн яхцмюкю
Fs_i = input('Fs_i = ');    % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
L = input('L = ');          % бейрнп йнщттхжхемрнб хмрепонкъжхх
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
disp('%')
disp('% дКЪ БШБНДЮ ГМЮВЕМХИ Ni, Fs_i Х L МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp(['      N = Ni = ',num2str(Ni),'      Fs = Fs_i = ',num2str(Fs_i)])
disp('%')
disp(['      L = [',num2str(L) ']'])
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб бундмнцн Х бшундмшу ЯХЦМЮКНБ МЮФЛХРЕ <ENTER>')
pause
N = Ni;           % оепхнд бундмнцн яхцмюкю
Fs = Fs_i;        % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
n = 0:(N-1);      % дхяйпермне мнплхпнбюммне бпелъ бундмнцн яхцмюкю
x = A1*sin((2*pi*f1/Fs).*n)+A2*sin((2*pi*f2/Fs).*n);% бундмни яхцмюк
figure('Name','Input and Output Signals in Interpolation System','NumberTitle', 'off')
subplot(4,1,1), stem(n,x), grid, xlabel('n')
title(strcat(['Input Signal x(n)   N = ',num2str(N)]))
for i = 1:length(L)                % хмдейяш бейрнпю L
[y{i},h{i}] = interp(x,L(i));      % бшундмни яхцмюк х ху йху-тхкэрпю (cell array)
subplot(4,1,i+1), stem(0:length(y{i})-1,y{i}), grid, xlabel('n')
title(strcat(['L = ',num2str(L(i)),'   Output Signal y(n)   length(y) = ',num2str(length(y{i}))]))
end
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб хлоскэямшу уюпюйрепхярхй йху-ТХКЭРПНБ МЮФЛХРЕ <ENTER>')
pause
figure('Name','Impulse Responses FIR in Interpolation System','NumberTitle', 'off')
for i = 1:length(L)                 % хмдейяш бейрнпю L
subplot(3,1,i), stem(0:length(h{i})-1,   h{i},'MarkerSize',3,'Linewidth',2), grid, xlabel('n')
title(strcat(['L = ',num2str(L(i)),'   Impulse Response h(n)   length(h) = ',num2str(length(h{i}))]))
end
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.2. бшвхякемхе юлокхрсдмшу яоейрпнб яхцмюкнб х юву йху-тхкэрпнб')
disp('% яхярелш ндмнйпюрмни хмрепонкъжхх')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юлокхрсдмшу яоейрпнб Х юву йху-ТХКЭРПНБ')
disp('% Б "мнбни" ЬЙЮКЕ ВЮЯРНР МЮФЛХРЕ <ENTER>')
pause
X = fft(x);                    % дот бундмнцн яхцмюкю
MODX = (2/N)*abs(X);           % юлокхрсдмши яоейрп бундмнцн яхцмюкю
MODX(1) = (1/N)*abs(X(1));
k = 0:N-1;       % дхяйпермше мнплхпнбюммше вюярнрш бундмнцн яхцмюкю
figure('Name','Amplitude Spectrums and Magnitude Responses ("new" frequencies)','NumberTitle', 'off')
subplot(4,1,1), stem(k*(Fs/N),MODX,'MarkerSize',3,'Linewidth',2)
grid, xlabel('f (Hz)'), title(strcat(['Amplitude Spectrum x(n)   N = ',num2str(N)]))
f = 0:Fs/1000:Fs;    % яерйю вюярнр дкъ юву йху-тхкэрпю б "мнбни" ьйюке
for i = 1:length(L)          % хмдейяш бейрнпю L
Y{i} = fft(y{i});            % дот бшундмнцн яхцмюкю(cell array)
Y_D = Y{i};                  % дот бшундмнцн яхцмюкю(double)
MODY_D = (2/length(Y_D))*abs(Y_D); % юлокхрсдмши яоейрп бшундмнцн яхцмюкю(double) 
MODY_D(1) = (1/length(Y_D))*abs(Y_D(1));
MODY{i} = MODY_D;  % юлокхрсдмши яоейрп бшундмнцн яхцмюкю(cell array)
MAG = abs(freqz(h{i},1,f,Fs));   % юву йху-тхкэрпю
k_out = 0:length(Y{i})-1;        % дхяйпермше мнплхпнбюммше вюярнрш бшундмнцн яхцмюкю
subplot(4,1,i+1), stem(k_out*(Fs/length(Y{i})),MODY{i},'MarkerSize',3,'Linewidth',2)
grid, xlabel('f ▓ (Hz) ("new" frequencies)')
title(strcat(['L = ',num2str(L(i)),'   Amplitude Spectrum y(n) and FIR Magnitude Response']))
hold on, plot(f,MAG,'r','Linewidth',2)
end
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юлокхрсдмшу яоейрпнб Х юву йху-ТХКЭРПНБ')
disp('% Б "ярюпни" ЬЙЮКЕ ВЮЯРНР МЮФЛХРЕ <ENTER>')
pause
figure('Name','Amplitude Spectrums and Magnitude Responses ("old" frequencies)','NumberTitle', 'off')
subplot(4,1,1), stem(k*(Fs/N),MODX,'MarkerSize',3,'Linewidth',2)
grid, xlabel('f (Hz)'), xlim([0 max(L)*Fs])
title(strcat(['Amplitude Spectrum x(n)    N=',num2str(N)]))
for i = 1:length(L)                   % хмдейяш бейрнпю L
f = 0:(Fs*L(i))/1000:(Fs*L(i));       % яерйю вюярнр дкъ юву йху-тхкэрпю б "ярюпни" ьйюке
MAG = abs(freqz(h{i},1,f,(Fs*L(i)))); % юву йху-тхкэрпю
k_out = 0:length(Y{i})-1;             % дхяйпермше мнплхпнбюммше вюярнрш бшундмнцн яхцмюкю
subplot(4,1,i+1),stem(k_out*(Fs*L(i)/length(Y{i})),MODY{i},'MarkerSize',3,'Linewidth',2)
grid, xlabel('f (Hz) ("old" frequencies)'), xlim([0 max(L)*Fs])
title(strcat(['L = ',num2str(L(i)),'   Amplitude Spectrum y(n) and FIR Magnitude Response']))
hold on, plot(f,MAG,'r','Linewidth',2)
end
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.3. лндекхпнбюмхе яхярелш ндмнйпюрмни дежхлюжхх')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше');
DATA=0;
while DATA==0
Nd = input('Nd = ');        % оепхнд бундмнцн яхцмюкю
Fs_d = input('Fs_d = ');    % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
M = input('M = ');          % бейрнп йнщттхжхемрнб дежхлюжхх
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
disp('%')
disp('% дКЪ БШБНДЮ ГМЮВЕМХИ Nd, Fs_d Х M МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp(['      	N = Nd = ',num2str(Nd),'      Fs = Fs_d = ',num2str(Fs_d)])
disp('%')
disp(['      	M = [',num2str(M) ']'])
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб бундмнцн Х бшундмшу ЯХЦМЮКНБ МЮФЛХРЕ <ENTER>')
pause
N = Nd;            % оепхнд бундмнцн яхцмюкю
Fs = Fs_d;         % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
n = 0:(N-1);       % дхяйпермне мнплхпнбюммне бпелъ бундмнцн яхцмюкю
x = A1*sin((2*pi*f1/Fs).*n)+A2*sin((2*pi*f2/Fs).*n);% бундмни яхцмюк
figure('Name','Input and Output Signals in Decimation System','NumberTitle', 'off')
subplot(4,1,1), stem(n,x), grid, xlabel('n')
title(strcat(['Input Signal x(n)      N = ',num2str(N)]))
for i = 1:length(M)              % хмдейяш бейрнпю M
y{i} = decimate(x,M(i),'fir');   % бшундмни яхцмюк(cell array)
subplot(4,1,i+1), stem(0:length(y{i})-1,y{i}), grid, xlabel('n')
title(strcat(['M = ',num2str(M(i)),'   Output Signal y(n)   length(y) = ',num2str(length(y{i}))]))
end
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.4. бшвхякемхе юлокхрсдмшу яоейрпнб яхцмюкнб')
disp('% яхярелш ндмнйпюрмни дежхлюжхх')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юлокхрсдмшу яоейрпнб МЮФЛХРЕ <ENTER>')
pause
X = fft(x);                     % дот бундмнцн яхцмюкю
MODX = (2/N)*abs(X);            % юлокхрсдмши яоейрп бундмнцн яхцмюкю
MODX(1) = (1/N)*abs(X(1));
k = 0:N-1;         % дхяйпермше мнплхпнбюммше вюярнрш бундмнцн яхцмюкю
figure('Name','Amplitude Spectrums in Decimation System','NumberTitle', 'off')
subplot(4,1,1), stem(k*(Fs/N),MODX,'MarkerSize',3,'Linewidth',2)
grid, xlabel('f (Hz)'),title(strcat(['Amplitude Spectrum x(n)   N = ',num2str(N)])), xlim([0 Fs])
for i = 1:length(M)                   % хмдейяш бейрнпю M 
Y{i} = fft(y{i});                     % дот бшундмнцн яхцмюкю
Y_D = Y{i};                           % дот бшундмнцн яхцмюкю(double)
MODY_D = (2/length(Y_D))*abs(Y_D);    % юлокхрсдмши яоейрп бшундмнцн яхцмюкю(double)
MODY_D(1) = (1/length(Y_D))*abs(Y_D(1));
MODY{i} = MODY_D; % юлокхрсдмши яоейрп бшундмнцн яхцмюкю(cell array)
k_out = 0:length(Y{i})-1;             % дхяйпермше мнплхпнбюммше вюярнрш бшундмнцн яхцмюкю
subplot(4,1,i+1),stem(k_out*(Fs/M(i)/length(Y{i})),MODY{i},'MarkerSize',3, 'Linewidth',2)
grid, xlabel('f (Hz)'), xlim([0 Fs])
title(strcat(['M = ',num2str(M(i)),'   Amplitude Spectrum y(n)']))
end
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.5. лндекхпнбюмхе яхярелш ндмнйпюрмни оепедхяйперхгюжхх')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше');
DATA=0;
while DATA==0
Nr = input('Nr = ');        % оепхнд бундмнцн яхцмюкю
Fs_r = input('Fs_r = ');    % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
L = input('L = ');          % йнщттхжхемр хмрепонкъжхх
M = input('M = ');          % йнщттхжхемр дежхлюжхх
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
disp('%')
disp('% дКЪ БШБНДЮ ГМЮВЕМХИ Nr, Fs_r Х L/M МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp(['      N = Nr = ',num2str(Nr),'      Fs = Fs_r = ',num2str(Fs_r)])
disp('%')
disp(['      	L/M = ',num2str(L),'/',num2str(M)])
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб бундмнцн Х бшундмшу ЯХЦМЮКНБ МЮФЛХРЕ <ENTER>')
pause
N = Nr;              % оепхнд бундмнцн яхцмюкю
Fs = Fs_r;           % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
n = 0:(N-1);         % дхяйпермне мнплхпнбюммне бпелъ бундмнцн яхцмюкю
x = A1*sin((2*pi*f1/Fs).*n)+A2*sin((2*pi*f2/Fs).*n);  % бундмни яхцмюк
figure('Name','Input and Output Signals in Resampling System','NumberTitle', 'off')
subplot(2,1,1), stem(n,x), grid, xlabel('n')
title(strcat(['Input Signal x(n)      N = ',num2str(N)]))
[y h] = resample(x,L,M);            % бшундмни яхцмюк х ху йху-тхкэрпю
subplot(2,1,2), stem(0:length(y)-1,y), grid, xlabel('n')
title(strcat(['L/M = ',num2str(L),'/',num2str(M),'  Output Signal y(n)   length(y) = ',num2str(length(y))]))
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.6. бшвхякемхе юлокхрсдмшу яоейрпнб яхцмюкнб х юву йху-тхкэрпю')
disp('% яхярелш ндмнйпюрмни оепедхяйперхгюжхх')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юлокхрсдмшу яоейрпнб Х юву йху-ТХКЭРПЮ МЮФЛХРЕ <ENTER>')
pause
figure('Name','Amplitude Spectrums and Magnitude Responses in Resampling System','NumberTitle', 'off')
X = fft(x);                     % дот бундмнцн яхцмюкю
MODX = (2/N)*abs(X);            % юлокхрсдмши яоейрп бундмнцн яхцмюкю
MODX(1) = (1/N)*abs(X(1));
k = 0:N-1;            % дхяйпермше мнплхпнбюммше вюярнрш бундмнцн яхцмюкю
subplot(2,1,1),stem(k*(Fs/N),MODX,'MarkerSize',3,'Linewidth',2)
grid, xlabel('f (Hz)'), title(strcat(['Amplitude Spectrum x(n)   N = ',num2str(N)])), xlim([0 (L/M)*Fs])
Y = fft(y);                         % дот бшундмнцн яхцмюкю
MODY = (2/length(Y))*abs(Y);        % юлокхрсдмши яоейрп бшундмнцн яхцмюкю
MODY(1) = (1/length(Y))*abs(Y(1));
k_out = 0:length(Y)-1;              % дхяйпермше мнплхпнбюммше вюярнрш бшундмнцн яхцмюкю
f = 0:Fs*(L/M)/1000:(Fs*(L/M));     % яерйю вюярнр дкъ юву йху-тхкэрпю
MAG = abs(freqz(h,1,f,(Fs*(L/M)))); % юву йху-тхкэрпю
subplot(2,1,2), stem(k_out*(Fs*(L/M)/length(Y)),MODY,'MarkerSize',3, 'Linewidth',2)
grid, xlabel('f (Hz)')
title(strcat(['L/M = ',num2str(L),'/',num2str(M),'  Amplitude Spectrum y(n) and FIR Magnitude Response']))
hold on, plot(f,MAG,'r','Linewidth',2), xlim([0 (L/M)*Fs])
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.7. лндекхпнбюмхе яхярелш ндмнйпюрмни хмрепонкъжхх я онкхтюгмни ярпсйрспни')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше');
DATA=0;
while DATA==0
Nip = input('Nip = ');       % оепхнд бундмнцн яхцмюкю
Fs_ip = input('Fs_ip = ');   % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
L = input('L = ');           % йнщттхжхемр хмрепонкъжхх
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
disp('%')
disp('% дКЪ БШБНДЮ ГМЮВЕМХИ Ni, Fs_ip Х L МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp(['      	N = Nip = ',num2str(Nip),'      Fs = Fs_ip = ',num2str(Fs_ip)])
disp('%')
disp(['      	L = ',num2str(L)])
N = Nip;                % оепхнд бундмнцн яхцмюкю
Fs = Fs_ip;             % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
n = 0:(N-1);            % дхяйпермне мнплхпнбюммне бпелъ бундмнцн яхцмюкю
x = A1*sin((2*pi*f1/Fs).*n)+A2*sin((2*pi*f2/Fs).*n);% бундмни яхцмюк
disp('%')
disp('% дКЪ БШБНДЮ ябниярб назейрю Hi МЮФЛХРЕ <ENTER>')
pause
Hi = mfilt.firinterp(L)          % нохяюмхе онкхтюгмни ярсйрспш яхярелш ндмнйпюрмни хмрепонкъжхх
y = filter(Hi,x);                % бшундмни яхцмюк
ni_start = ceil(length(y)/2+1);  % мювюкн сярюмнбхбьецняъ пефхлю
yi = y(ni_start:end);            % бшундмни яхцмюк б сярюмнбхбьеляъ пефхле
ni = ni_start:(ni_start+length(yi)-1);  % дхяйпермне мнплхпнбюммне бпелъ дкъ сярюмнбхбьецняъ пефхлю
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб бундмнцн Х бшундмнцн ЯХЦМЮКНБ МЮФЛХРЕ <ENTER>')
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
X = fft(x);                    % дот бундмнцн яхцмюкю
MODX = (2/N)*abs(X);           % юлокхрсдмши яоейрп бундмнцн яхцмюкю
MODX(1) = (1/N)*abs(X(1));
Yi = fft(yi);                    % дот бшундмнцн яхцмюкю
MODYi = (2/length(Yi))*abs(Yi);  % юлокхрсдмши яоейрп бшундмнцн яхцмюкю
MODYi(1) = (1/length(Yi))*abs(Yi(1));
k = 0:(N-1);           % дхяйпермше мнплхпнбюммше вюярнрш бундмнцн яхцмюкю
ki = 0:(length(Yi)-1); % дхяйпермше мнплхпнбюммше вюярнрш бшундмнцн яхцмюкю
f = 0:(Fs*L)/1000:Fs*L;           % яерйю вюярнр дкъ юву йху-тхкэрпю
h = Hi.Numerator;                 % ху йху-тхкэрпю
MAG = abs(freqz(h,1,f,Fs*L));     % юву йху-тхкэрпю
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
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.8. лндекхпнбюмхе яхярелш ндмнйпюрмни дежхлюжхх я онкхтюгмни ярпсйрспни')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше');
DATA=0;
while DATA==0
Ndp = input('Ndp = ');          % оепхнд бундмнцн яхцмюкю
Fs_dp = input('Fs_dp = ');      % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
M = input('M = ');              % йнщттхжхемр дежхлюжхх
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
disp('%')
disp('% дКЪ БШБНДЮ ГМЮВЕМХИ Ndp, Fs_dp Х M МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp(['      	N = Ndp = ',num2str(Ndp),'      Fs = Fs_dp = ',num2str(Fs_dp)])
disp('%')
disp(['      	M = ',num2str(M)])
N = Ndp;               % оепхнд бундмнцн яхцмюкю
Fs = Fs_dp;            % вюярнрю дхяйперхгюжхх бундмнцн яхцмюкю
n = 0:(N-1);           % дхяйпермне мнплхпнбюммне бпелъ бундмнцн яхцмюкю
x = A1*sin((2*pi*f1/Fs).*n)+A2*sin((2*pi*f2/Fs).*n);    % бундмни яхцмюк
disp('%')
disp('% дКЪ БШБНДЮ ябниярб назейрю Hd МЮФЛХРЕ <ENTER>')
pause
Hd = mfilt.firdecim(M) % нохяюмхе онкхтюгмни ярсйрспш яхярелш ндмнйпюрмни дежхлюжхх
y = filter(Hd,x);                % бшундмни яхцмюк
nd_start = ceil(length(y)/2+1);  % мювюкн сярюмнбхбьецняъ пефхлю
yd = y(nd_start:end);            % бшундмни яхцмюк б сярюмнбхбьеляъ пефхле
nd = nd_start:(nd_start+length(yd)-1);  % дхяйпермне мнплхпнбюммне бпелъ дкъ сярюмнбхбьецняъ пефхлю
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб бундмнцн Х бшундмнцн ЯХЦМЮКНБ МЮФЛХРЕ <ENTER>')
pause
figure('Name','Input and Output Signals in Polyphase Structure Decimation System','NumberTitle', 'off')
subplot(3,1,1), stem(n,x), grid, xlabel('n')
title(strcat(['Input Signal x(n)   N = ',num2str(N)]))
subplot(3,1,2), stem(0:length(y)-1,y), grid
title(strcat(['M = ',num2str(M),'  Output Signal y(n)   length(y) = ',num2str(length(y))]))
subplot(3,1,3), stem(nd,yd), grid, xlabel('n')
title(strcat(['Output Signal y(n)   Starting point n = ',num2str(nd_start),'   length(y) = ',num2str(length(yd))]))
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юлокхрсдмшу яоейрпнб Х юву йху-ТХКЭРПЮ МЮФЛХРЕ <ENTER>')
pause
X = fft(x);                         % дот бундмнцн яхцмюкю
MODX = (2/N)*abs(X);                % юлокхрсдмши яоейрп бундмнцн яхцмюкю
MODX(1) = (1/N)*abs(X(1));
Yd = fft(yd);              % дот бшундмнцн яхцмюкю б сярюмнбхбьеляъ пефхле
MODYd = (2/length(Yd))*abs(Yd);     % юлокхрсдмши яоейрп бшундмнцн яхцмюкю
MODYd(1) = (1/length(Yd))*abs(Yd(1));
k = 0:(N-1);           % дхяйпермше мнплхпнбюммше вюярнрш бундмнцн яхцмюкю
kd = 0:(length(Yd)-1); % дхяйпермше мнплхпнбюммше вюярнрш бшундмнцн яхцмюкю
f = 0:Fs/1000:Fs;      % яерйю вюярнр дкъ дкъ юву йху-тхкэрпю
hd = Hd.Numerator;     % ху йху-тхкэрпю
MAG = abs(freqz(hd,1,f,Fs)); % юву йху-тхкэрпю
figure('Name','Amplitude Spectrums and Magnitude Response in Polyphase Structure Decimation System','NumberTitle', 'off')
subplot(2,1,1), stem(k*(Fs/N),MODX,'MarkerSize',3,'Linewidth',2)
grid, xlabel('f (Hz)'), xlim([0 Fs])
title(strcat('Amplitude Spectrum x(n) and FIR Magnitude Response'))
hold on, plot(f,MAG,'r','Linewidth',2)
subplot(2,1,2), stem(kd*((Fs/M)/length(Yd)),MODYd,'MarkerSize',3, 'Linewidth',2)
grid, xlabel('f (Hz)'), xlim([0 Fs])
title(strcat(['M = ',num2str(M),'   Amplitude Spectrum y(n)']))
disp('%')
disp('%')
disp('% пюанрю гюбепьемю')
























    

















































































