clc
clear
disp('% кп ╧21. юдюорхбмше тхкэрпш')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше')
DATA=0;
while DATA==0
Nb = input('Nb = ');      % мнлеп апхцюдш
N = input('N = ');        % дкхмю йху-тхкэрпю б янярюбе ют
L = input('L = ');        % дкхмю бундмнцн яхцмюкю ют
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.1. лндекхпнбюмхе мнплюкэмнцн аекнцн ьслю')
n = 0:(L-1);              %дхяйпермне мнплхпнбюммне бпелъ
r_gauss = randn(1,L);     % мнплюкэмши аекши ьсл
disp('%')
disp('%')
disp('% яНГДЮММЮЪ ЛНДЕКЭ мнплюкэмнцн аекнцн ьслю АСДЕР ХЯОНКЭГНБЮМЮ Б ДЮКЭМЕИЬХУ ХЯЯКЕДНБЮМХЪУ')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.2. лндекхпнбюмхе ярпсйрспш ют я юкцнпхрлнл LMS')
x = r_gauss;                % бундмни яхцмюк ют
Px = var(x);                % япедмхи йбюдпюр бундмнцн яхцмюкю ют
mu_max = 2/(N*Px);          % люйяхлюкэмши ьюц юдюорюжхх
mu = 0.5*mu_max;            % гюдюммши ьюц юдюорюжхх
Hlms = adaptfilt.lms(N,mu)  % ярпсйрспю ют я юкцнпхрлнл LMS
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.3. лндекхпнбюмхе ярпсйрспш ют я юкцнпхрлнл NLMS')
epsilon = 1e-6; % йнмярюмрю, нопедекъчыюъ люйяхлюкэмши ьюц юдюорюжхх
Hnlms = adaptfilt.nlms(N,1,1,epsilon)    % ярпсйрспю ют я юкцнпхрлнл NLMS
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.4. лндекхпнбюмхе ярпсйрспш ют я юкцнпхрлнл RLS')
Hrls = adaptfilt.rls(N)                 % ярпсйрспю ют я юкцнпхрлнл RLS
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.5. нжемйю хлоскэямни уюпюйрепхярхйх мехгбеярмни кдя')
disp('%')
disp('%')
disp('% 5.1. нжемйю хлоскэямни уюпюйрепхярхйх мехгбеярмни кдя - йху-тхкэрпю тмв (FIR)')
disp('%')
disp('%')
x = r_gauss;         % бундмни яхцмюк мехгбеярмни кдя
R1 = round(N/2);     % онпъднй йху-тхкэрпю тмв
wc = 0.5;            % мнплхпнбюммюъ вюярнрю пюгпшбю йху-тхкэрпю тмв
b = fir1(R1,wc);     % йнщттхжхемрш йху-тхкэрпю тмв
d = filter(b,1,x);   % бшундмни яхцмюк мехгбеярмни кдя (FIR)
h = b;               % хярхммюъ хлоскэямюъ уюпюйрепхярхйю мехгбеярмни кдя (FIR) дкхмш (R1+1)
[y_lms,e_lms] = filter(Hlms,x,d); % бшундмни яхцмюк х яхцмюк ньхайх ют я юкцнпхрлнл LMS
[y_rls,e_rls] = filter(Hrls,x,d); % бшундмни яхцмюк х яхцмюк ньхайх ют я юкцнпхрлнл RLS
h_lms = Hlms.coefficients;        % нжемйю хлоскэямни уюпюйрепхярхйх мехгбеярмни кдя (FIR) - оюпюлерпш ют я юкцнпхрлнл LMS
h_rls = Hrls.coefficients;        % нжемйю хлоскэямни уюпюйрепхярхйх мехгбеярмни кдя (FIR) - оюпюлерпш ют я юкцнпхрлнл RLS
disp('% дКЪ БШБНДЮ цпютхйнб яхцмюкю ньхайх ют МЮФЛХРЕ <ENTER>')
pause
n = 0:length(x)-1;                % дхяйпермне мнплхпнбюммне бпелъ дкъ яхцмюкю ньхайх ют
figure('Name','Error signal of AF for LMS and RLS','NumberTitle', 'off')
subplot(2,1,1), plot(n,e_lms), grid, xlabel('n'), title('Error signal for LMS')
subplot(2,1,2), plot(n,e_rls), grid, xlabel('n'), title('Error signal for RLS')
n1 = 0:N-1;            % хмрепбюк дхяйпермнцн мнплхпнбюммнцн бпелемх дкъ нжемнй хлоскэямни уюпюйрепхярхйх
if length(h)<N
    h = [h zeros(1,(N-length(h)))];  % хярхммюъ хлоскэямюъ уюпюйрепхярхйю мехгбеярмни кдя (FIR), днонкмеммюъ мскълх дн дкхмш N
end
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб хярхммни хлоскэямни уюпюйрепхярхйх (FIR) Х ЕЕ нжемнй МЮФЛХРЕ <ENTER>')
pause
figure('Name','True Impulse response FIR and its Estimates','NumberTitle', 'off')
subplot(3,1,1), stem(n1,h), grid, xlabel('n1')
title('True Impulse response FIR - h(n) with length N')
subplot(3,1,2), stem(n1,h_lms), grid, xlabel('n1')
title('LMS Impulse response FIR - h lms')
subplot(3,1,3), stem(n1,h_rls), grid, xlabel('n1')
title('RLS Impulse response FIR - h rls')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ япедмецн юаянкчрмнцн нрйкнмемхъ')
disp('% нрявернб хлоскэямни уюпюйрепхярхйх НР ЕЕ нжемнй (FIR) МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp(['         norm1_lms = ',num2str((1/N)*norm(h-h_lms,1))])
disp(['         norm1_rls = ',num2str((1/N)*norm(h-h_rls,1))])
disp('%')
disp('%')
disp('% 5.2. нжемйю хлоскэямни уюпюйрепхярхйх мехгбеярмни кдя - аху-тхкэрпю тбв (IIR)')
disp('%')
disp('%')
R1 = round(N/2);       % онпъднй аху-тхкэрпю тбв
WDn = 0.3;             % мнплхпнбюммюъ вюярнрю япегю аху-тхкэрпю тбв
[b,a] = butter(R1,WDn,'high');     % йнщттхжхемрш аху-тхкэрпю тбв
d = filter(b,a,x);     % бшундмни яхцмюк мехгбеярмни кдя (IIR)
h = impz(b,a,N);       % хярхммюъ хлоскэямюъ уюпюйрепхярхйю мехгбеярмни кдя (IIR) дкхмш N (бейрнп-ярнкаеж)
[y_lms,e_lms] = filter(Hlms,x,d);  % бшундмни яхцмюк х яхцмюк ньхайх ют я юкцнпхрлнл LMS
[y_rls,e_rls] = filter(Hrls,x,d);  % бшундмни яхцмюк х яхцмюк ньхайх ют я юкцнпхрлнл RLS
h_lms = Hlms.coefficients;         % оюпюлерпш ют я юкцнпхрлнл LMS
h_rls = Hrls.coefficients;         % оюпюлерпш ют я юкцнпхрлнл RLS
disp('% дКЪ БШБНДЮ цпютхйнб яхцмюкю ньхайх ют МЮФЛХРЕ <ENTER>')
pause
figure('Name','Error signal of AF for LMS and RLS','NumberTitle', 'off')
subplot(2,1,1), plot(n,e_lms), grid, xlabel('n'), title('Error signal for LMS')
subplot(2,1,2), plot(n,e_rls), grid, xlabel('n'), title('Error signal for RLS')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб хярхммни хлоскэямни уюпюйрепхярхйх (IIR) Х ЕЕ нжемнй МЮФЛХРЕ <ENTER>')
pause
figure('Name','True Impulse response IIR and its Estimates','NumberTitle', 'off')
subplot(3,1,1), stem(n1,h), grid, xlabel('n1')
title('True Impulse response IIR - h(n) with length N')
subplot(3,1,2), stem(n1,h_lms), grid, xlabel('n1')
title('LMS Impulse response IIR - h lms(n)')
subplot(3,1,3), stem(n1,h_rls), grid, xlabel('n1')
title('RLS Impulse response IIR - h rls(n)')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ япедмецн юаянкчрмнцн нрйкнмемхъ')
disp('% нрявернб хлоскэямни уюпюйрепхярхйх НР ЕЕ нжемнй (IIR) МЮФЛХРЕ <ENTER>')
disp('%')
disp(['         norm1_lms = ',num2str((1/N)*norm(h'-h_lms,1))])
disp(['         norm1_rls = ',num2str((1/N)*norm(h'-h_rls,1))])
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.6. нвхярйю яхцмюкю нр ьслю')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше')
DATA=0;
while DATA==0
Fs = input('Fs = ');      % вюярнрю дхяйперхгюжхх (цЖ)
A1 = input('A1 = ');      % юлокхрсдш дхяйпермшу цюплнмхй
A2 = input('A2 = ');
f1 = input('f1 = ');      % вюярнрш (цЖ) дхяйпермшу цюплнмхй
f2 = input('f2 = ');
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
x_noise = r_gauss;   % бундмни яхцмюк ют - мнплюкэмши аекши ьсл
b = fir1(R1,wc);     % йнщттхжхемрш йху-тхкэрпю тмв, хяйюфючыецн ьсл
x_noiseNEW = filter(b,1,x_noise);   % хяйюфеммши ьсл мю бшунде йху-тхкэрпю тмв
w1 = 2*pi*f1/Fs;  w2 = 2*pi*f2/Fs;  % мнплхпнбюммше вюярнрш дхяйпермшу цюплнмхй (пюд)
n = 0:(L-1);                        % дхяйпермне мнплхпнбюммне бпелъ
s = A1*cos(w1*n)+A2*cos(w2*n);      % онкегмши яхцмюк
d = s+x_noiseNEW;                   % напюгжнбши яхцмюк ют
[y_lms,e_lms] = filter(Hlms,x_noise,d); % бшундмни яхцмюк х яхцмюк ньхайх ют я юкцнпхрлнл LMS
[y_rls,e_rls] = filter(Hrls,x_noise,d); % бшундмни яхцмюк х яхцмюк ньхайх ют я юкцнпхрлнл RLS
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб онкегмнцн яхцмюкю, ецн юддхрхбмни ялеях я ьслнл')
disp('% Х нжемнй онкегмнцн яхцмюкю (LMS) х (RLS) МЮФЛХРЕ <ENTER>')
pause
figure('Name','Harmonic Signal, Mixture of Signal and Noise, Estimates of Harmonic Signal (LMS and RLS)','NumberTitle', 'off')
subplot(4,1,1), plot(n,s), grid, title('Harmonic Signal - s(n)')
subplot(4,1,2), plot(n,d), grid, title('Mixture of Harmonic Signal and Noise - d(n)')
subplot(4,1,3), plot(n,e_lms), grid, title('Estimate Harmonic Signal (LMS)')
subplot(4,1,4), plot(n,e_rls), grid, xlabel('n'), title('Estimate Harmonic Signal (RLS)')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юлокхрсдмшу яоейрпнб онкегмнцн яхцмюкю, ецн юддхрхбмни ялеях я ьслнл')
disp('% Х нжемнй онкегмнцн яхцмюкю (LMS) х (RLS)')
disp('% Б сярюмнбхбьеляъ пефхле МЮФЛХРЕ <ENTER>')
pause
k = 0:(L-1);                  % дхяйпермюъ мнплхпнбюммюъ вюярнрю
f = (k/L)*Fs;                 % вюярнрю (цЖ)
MOD_s = (2/L)*abs(fft(s));    % юлокхрсдмши яоейрп онкегмнцн яхцмюкю
MOD_s(1) = (1/L)*abs(fft(s(1)));
MOD_d = (2/L)*abs(fft(d));    % юлокхрсдмши яоейрп ялеях онкегмнцн яхцмюкю я ьслнл
MOD_d(1) = (1/L)*abs(fft(d(1)));
n_start = round(0.05*L);      % мювюкн сярюмнбхбьецняъ пефхлю
Ls = (L- n_start)+1;          % дкхмю нжемнй яхцмюкю х ьслю б сярюмнбхбьеляъ пефхле
e1_lms = e_lms(n_start:end);  % нжемйю онкегмнцн яхцмюкю б сярюмнбхбьеляъ пефхле (LMS)
e1_rls = e_rls(n_start:end);  % нжемйю онкегмнцн яхцмюкю б сярюмнбхбьеляъ пефхле (RLS)
ks = 0:(Ls-1);                % дхяйпермюъ мнплхпнбюммюъ вюярнрю дкъ сярюмнбхбьецняъ пефхлю
fs = (ks/Ls)*Fs;              % вюярнрю (цЖ)
MOD_lms = (2/Ls)*abs(fft(e1_lms));     % юлокхрсдмши яоейрп нжемйх онкегмнцн яхцмюкю б сярюмнбхбьеляъ пефхле (LMS)
MOD_lms(1) = (1/Ls)*abs(fft(e1_lms(1)));
MOD_rls = (2/Ls)*abs(fft(e1_rls));     % юлокхрсдмши яоейрп нжемйх онкегмнцн яхцмюкю б сярюмнбхбьеляъ пефхле (RLS)
MOD_rls(1) = (1/Ls)*abs(fft(e1_rls(1)));
figure('Name','Amplitude Spectrums of Harmonic Signal, Mixture of Signal and Noise, and Estimates of Harmonic Signal (LMS and RLS)','NumberTitle', 'off')
subplot(4,1,1), stem(f,MOD_s), grid, title('Amplitude Spectrum of Harmonic Signal')
subplot(4,1,2), stem(f,MOD_d), grid, title('Amplitude Spectrum of Harmonic Signal and Noise - d(n)')
subplot(4,1,3), stem(fs,MOD_lms), grid, title ('Amplitude Spectrum of Estimate Harmonic Signal (LMS) ')
subplot(4,1,4), stem(fs,MOD_rls), grid, xlabel('f'), title('Amplitude Spectrum of Estimate Harmonic (RLS)')
ns = 0:(Ls-1);                 % дхяйпермне мнплхпнбюммне бпелъ дкъ сярюмнбхбьецняъ пефхлю
s1 = s(n_start:end);           % онкегмши яхцмюк мю хмрепбюке бпелемх дкъ ецн нжемнй
RMSE_lms = sqrt((1/Ls).*sum((s1-e1_lms).^2));     % RMSE дкъ нжемйх онкегмнцн яхцмюкю (LMS)
RMSE_rls = sqrt((1/Ls).*sum((s1-e1_rls).^2));     % RMSE дкъ нжемйх онкегмнцн яхцмюкю (RLS)
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ гмювемхи RMSE МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp(['         RMSE_lms =     ',num2str(RMSE_lms)])
disp(['         RMSE_rls =     ',num2str(RMSE_rls)])
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.7. бшпюбмхбюмхе вюярнрмни уюпюйрепхярхйх мехгбеярмни кдя')
f2 = 3*f1;                  % хглемеммюъ вюярнрю (цЖ) цюплнмхйх
w2 = 2*pi*(3*f1)/Fs;        % мнплхпнбюммюъ дхяйпермюъ вюярнрю (пюд)
n = 0:(L-1);                % дхяйпермне мнплхпнбюммне бпелъ
s = A1*cos(w1*n)+A2*cos(w2*n)+3*max(A1,A2)*r_gauss; % бундмни яхцмюк мехгбеярмни кдя - напюгжнбши яхцмюк ют
k = 0:(L-1);                % дхяйпермюъ мнплхпнбюммюъ вюярнрю
f = k*(Fs/L);               % вюярнрю (цЖ)
S = fft(s);                 % дот бундмнцн яхцмюкю мехгбеярмни кдя
MODS = (2/L)*abs(S);        % юлокхрсдмши яоейрп бундмнцн яхцмюкю мехгбеярмни кдя
MODS(1) = (1/L)*abs(S(1));
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб бундмнцн яхцмюкю мехгбеярмни кдя Х ЕЦН юлокхрсдмнцн яоейрпю   МЮФЛХРЕ <ENTER>')
pause
figure('Name','Input Signal of Unknown LDS and its Amplitude Spectrum','NumberTitle', 'off')
subplot(2,1,1), plot(n,s), grid, xlabel('n'), title('Input Signal of Unknown LDS')
subplot(2,1,2), stem(f,MODS), grid, xlabel('f'), title('Amplitude Spectrum')
d = s;       % напюгжнбши яхцмюк ют - бундмни яхцмюк мехгбеярмни кдя
D = floor(N/2);                % бекхвхмю гюдепфйх напюгжнцн яхцмюкю
D_zeros = zeros(1,round(N/2)); % тнплхпнбюмхе мювюкэмшу мскебшу гмювемхи гюдепфюммнцн напюгжнбнцн яхцмюкю
d_delay = [D_zeros d(1:(length(d)-length(D_zeros)))];  % гюдепфюммши напюгжнбши яхцмюк
R2 = round(N/7);     % онпъднй йху-тхкэрпю тмв
wc = 0.5;            % мнплхпнбюммюъ вюярнрю пюгпшбю йху-тхкэрпю тмв
b = fir1(R2,wc);     % хлоскэямюъ уюпюйрепхярхйю мехгбеярмни кдя
x = filter(b,1,d);   % бшундмни яхцмюк мехгбеярмни кдя - бундмни яхцмюк ют
X = fft(x);          % дот бшундмнцн яхцмюкю мехгбеярмни кдя
MODX = (2/L)*abs(X); % юлокхрсдмши яоейрп бшундмнцн яхцмюкю мехгбеярмни кдя
MODX(1) = (1/L)*abs(X(1));
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб бшундмнцн яхцмюкю мехгбеярмни кдя Х ЕЦН юлокхрсдмнцн яоейрпю МЮФЛХРЕ <ENTER>')
pause
figure('Name','Output Signal of Unknown LDS and its Amplitude Spectrum','NumberTitle', 'off')
subplot(2,1,1), plot(n,x), grid, xlabel('n'), title('Output Signal of Unknown LDS')
subplot(2,1,2), stem(f,MODX), grid, xlabel('f'), title('Amplitude Spectrum')
Hrls = adaptfilt.rls(N);         % ярпсйрспю ют я юкцнпхрлнл RLS
[y_rls,e_rls] = filter(Hrls,x,d_delay);   % бшундмни яхцмюк х яхцмюк ньхайх ют я юкцнпхрлнл RLS
h_rls = Hrls.coefficients;       % оюпюлерпш ют (хлоскэямюъ уюпюйрепхярхйю ецн йху-тхкэрпю)
h_conv = conv(h_rls,b);          % хлоскэямюъ уюпюйрепхярхйю йюяйюдмнцн янедхмемхъ мехгбеярмни кдя х ют
n2 = 0:length(b)-1;              % дхяйпермне мнплхпнбюммне бпелъ дкъ хлоскэямни уюпюйрепхярхйх мехгбеярмни кдя (йху-тхкэрпю тмв)
n3 = 0:length(h_conv)-1;         % дхяйпермне мнплхпнбюммне бпелъ дкъ ябепрйх хлоскэямни уюпюйрепхярхйх мехгбеярмни кдя х ют
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб хлоскэямшу уюпюйрепхярхй мехгбеярмни кдя, йху-тхкэрпю б янярюбе ют')
disp('% Х ху йюяйюдмнцн янедхмемхъ  МЮФЛХРЕ <ENTER>')
pause
figure('Name','Impulse Response of AF, Unknown LDS, and Cascade Connection','NumberTitle', 'off')
subplot(3,1,1), stem(n2,b), grid, title('Impulse Response of Unknown LDS')
subplot(3,1,2), stem(n1,h_rls), grid, title('Impulse Response of AF')
subplot(3,1,3), stem(n3,h_conv), grid, xlabel('n'), title('Impulse Response of Cascade Connection')
fa = 0:Fs/500:(Fs/2);                 % бейрнп вюярнр дкъ юву
MAG_US = abs(freqz(b,1,fa,Fs));       % юву мехгбеярмни кдя
MAG_AF = abs(freqz(h_rls,1,fa,Fs));   % юву йху-тхкэрпю б янярюбе ют
MAG = MAG_US.*MAG_AF;                 % юву йюяйюдмнцн янедхмемхъ мехгбеярмни кдя х ют
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юву мехгбеярмни кдя, йху-тхкэрпю б янярюбе ют,')
disp('% Х ху йюяйюдмнцн янедхмемхъ МЮФЛХРЕ <ENTER>')
pause
figure('Name','Magnitude Response of Unknown LDS, AF, and Cascade Connection','NumberTitle', 'off')
subplot(3,1,1), plot(fa,MAG_US), grid, title('Magnitude Response of Unknown LDS')
subplot(3,1,2), plot(fa,MAG_AF), grid, title('Magnitude Response of AF')
subplot(3,1,3), plot(fa,MAG), grid, xlabel('f'), title('Magnitude Response of Cascade Connection')
PH_US = angle(freqz(b,1,fa,Fs));      % тву мехгбеярмни кдя
PH_AF = angle(freqz(h_rls,1,fa,Fs));  % тву йху-тхкэрпю б янярюбе ют
PH = PH_US+PH_AF;                     % тву йюяйюдмнцн янедхмемхъ мехгбеярмни кдя х ют
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб тву мехгбеярмни кдя, йху-тхкэрпю б янярюбе ют')
disp('% Х ху йюяйюдмнцн янедхмемхъ МЮФЛХРЕ <ENTER>')
pause
figure('Name','Phase Response of Unknown LDS, AF, and Cascade Connection','NumberTitle', 'off')
subplot(3,1,1), plot(fa,PH_US), grid, title('Phase Response of Unknown LDS')
subplot(3,1,2), plot(fa,PH_AF), grid, title('Phase Response of AF')
subplot(3,1,3), plot(fa,PH), grid, xlabel('f'), title('Phase Response of Cascade Connection')
Y = fft(y_rls);              % дот бшундмнцн яхцмюкю ют
MODY = (2/L)*abs(Y);         % юлокхрсдмши яоейрп бшундмнцн яхцмюкю ют
MODY(1) = (1/L)*abs(Y(1));
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб бшундмнцн яхцмюкю ют Х ЕЦН юлокхрсдмнцн яоейрпю МЮФЛХРЕ <ENTER>')
pause
figure('Name','Output Signal and Amplitude Spectrum of AF','NumberTitle', 'off')
subplot(2,1,1), plot(n,y_rls), grid, xlabel('n'), title('Output Signal of AF')
subplot(2,1,2), stem(f,MODY), grid, xlabel('f'), title('Amplitude Spectrum')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.8. бшвхякемхе нжемнй оюпюлерпнб юп-лндекх х нжемнй оюпюлерпнб кхмеимнцн опедяйюгюмхъ')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше')
DATA=0;
while DATA==0
a = input('a = ');              % бейрнп оюпюлерпнб юп-лндекх
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
e_AR = r_gauss;                 % бундмни яхцмюк юп-лндекх
y_AR = filter(1,a,e_AR);        % бшундмни яхцмюк юп-лндекх
x = y_AR;                       % юмюкхгхпселши яхцмюк онкюцюеряъ пюбмшл лндекхпселнлс яхцмюкс
p = length(a)-1;                % онпъднй юп-лндекх
[aYW DYW] = aryule(x,p);        % нжемйю оюпюлерпнб юп-лндекх он лерндс чкю-снкйепю
d = x;                          % напюгжнбши яхцмюк ют
x_delay = [0 x(1:(L-1))];       % бундмни яхцмюк ют - гюдепфюммши напюгжнбши яхцмюк
[y_rls,e_rls] = filter(Hrls,x_delay,d);  % бшундмни яхцмюк х яхцмюк ньхайх ют я юкцнпхрлнл RLS
h_rls = -Hrls.coefficients;      % бейрнп нжемнй оюпюлерпнб кхмеимнцн опедяйюгюмхъ
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб бейрнпю гюдюммшу оюпюлерпнб юп-лндекх, ХУ нжемнй')
disp('% Х бейрнпю нжемнй оюпюлерпнб кхмеимнцн опедяйюгюмхъ МЮФЛХРЕ <ENTER>')
pause
figure('Name','True AR parameters, Estimated AR parameters, and Estimated Linear Prediction parameters','NumberTitle', 'off')
YMAX = max([max(abs(a)) max(abs(aYW)) max(abs(h_rls))]); % люйяхлюкэмши он лндскч щкелемр бейрнпнб
subplot(3,1,1), stem(a(2:end)), grid, xlim([0 N-1]), ylim([-YMAX YMAX])
title('True parameters - a')
subplot(3,1,2), stem(aYW(2:end)), grid, xlim([0 N-1]), ylim([-YMAX YMAX])
title('Estimated parameters AR - aYW')
subplot(3,1,3), stem(h_rls), grid, xlim([0 N-1]), ylim([-YMAX YMAX])
title('Estimated Linear Prediction parameters - h rls')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ гмювемхи MAE МЮФЛХРЕ <ENTER>')
pause
a_N = [a(2:end) zeros(1,(N-length(a(2:end))))];       % бейрнп хярхммшу оюпюлерпнб юп-лндекх, днонкмеммши мскълх дн дкхмш N
aYW_N = [aYW(2:end) zeros(1,(N-length(aYW(2:end))))]; % бейрнп нжемнй оюпюлерпнб юп-лндекх, днонкмеммши мскълх дн дкхмш N
MAE_AR = (1/N).*sum(abs(a_N-aYW_N));      % MAE нжемнй оюпюлерпнб юп-лндекх
MAE_LP = (1/N).*sum(abs(h_rls-a_N));      % MAE нжемнй оюпюлерпнб кхмеимнцн опедяйюгюмхъ
disp('%')
disp(['         MAE_AR =     ',num2str(MAE_AR)])
disp(['         MAE_LP =     ',num2str(MAE_LP)])
disp('%')
disp('%')
disp('% пюанрю гюбепьемю')




























































































































