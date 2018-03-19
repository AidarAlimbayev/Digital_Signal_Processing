script
clc
clear
disp('% кп ╧10. дхяйпермне опенапюгнбюмхе тспэе (ВЮЯРЭ 2)')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше')
DATA=0;
while DATA==0
Nb = input('Nb = ');      % мнлеп апхцюдш
N = input('N = ');        % дкхмю (оепхнд) онякеднбюрекэмнярх
Fs = input('Fs = ');      % вюярнрю дхяйперхгюжхх (цЖ)
A1 = input('A1 = ');      % юлокхрсдш дхяйпермшу цюплнмхй
A2 = input('A2 = ');
f1 = input('f1 = ');      % вюярнрш дхяйпермшу цюплнмхй (цЖ)
f2 = input('f2 = ');
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
disp('% О.1. опнбепйю пюбемярбю оюпяебюкъ')
n = 0:(N-1);              % дхяйпермне мнплхпнбюммне бпелъ
k = 0:(N-1);              % дхяйпермюъ мнплхпнбюммюъ вюярнрю
w1 = 2*pi*f1/Fs; w2 = 2*pi*f2/Fs;  % мнплхпнбюммше вюярнрш дхяйпермшу цюплнмхй (пюд)
x = A1*cos(w1*n)+A2*cos(w2*n);     % онякеднбюрекэмнярэ(оепхнд N)
X = fft(x);               % дот онякеднбюрекэмнярх
E1 = sum(x.^2);           % щмепцхъ онякеднбюрекэмнярх, бшвхякеммюъ он ее нряверюл
E2 = (1/N)*sum(abs(X).^2); % щмепцхъ онякеднбюрекэмнярх, бшвхякеммюъ он нряверюл дот
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ КЕБНИ (E1) Х ОПЮБНИ (E2) ВЮЯРЕИ пюбемярбю оюпяебюкъ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp(['      E1 = ',num2str(E1),'      E2 = ' num2str(E2)])
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.2. хяякеднбюмхе щттейрю пюярейюмхъ яоейрпю дкъ ндмни дхяйпермни цюплнмхйх')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше')
DATA=0;
while DATA==0
M = input('M = ');    % оепхнд онякеднбюрекэмнярх M
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
n = 0:(N-1);          % дхяйпермне мнплхпнбюммне бпелъ (оепхнд N)
k = 0:(N-1);          % дхяйпермюъ мнплхпнбюммюъ вюярнрю (оепхнд N)
w1 = 2*pi*f1/Fs;      % мнплхпнбюммюъ вюярнрю (пюд)
x_N = A1*cos(w1*n);   % онякеднбюрекэмнярэ (оепхнд N)
X_N = fft(x_N);       % дот онякеднбюрекэмнярх (оепхнд N)
MOD_N = (2/N)*abs(X_N);  % юлокхрсдмши яоейрп онякеднбюрекэмнярх (оепхнд N)
MOD_N(1) = (1/N)*abs(X_N(1));
n1 = 0:(M-1);         % дхяйпермне мнплхпнбюммне бпелъ (оепхнд M)
k1 = 0:(M-1);         % дхяйпермюъ мнплхпнбюммюъ вюярнрю (оепхнд M)
x_M = A1*cos(w1*n1);  % онякеднбюрекэмнярэ (оепхнд M)
X_M = fft(x_M);          % дот онякеднбюрекэмнярх (оепхнд M)
MOD_M = (2/M)*abs(X_M);  % юлокхрсдмши яоейрп онякеднбюрекэмнярх (оепхнд M)
MOD_M(1) = (1/M)*abs(X_M(1));
P_N = N*f1/Fs;        % вхякн оепхнднб дхяйпермни цюплнмхйх я вюярнрни f1 мю оепхнде онякеднбюрекэмнярх N
P_M = M*f1/Fs;        % вхякн оепхнднб дхяйпермни цюплнмхйх я вюярнрни f1 мю оепхнде онякеднбюрекэмнярх M
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ вхякю оепхнднб ДХЯЙПЕРМНИ ЦЮПЛНМХЙХ я вюярнрни f1 МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp(['N = ',num2str(N),'  -->  P_N = ' num2str(P_N)])
disp(['M = ',num2str(M),'  -->  P_M = ' num2str(P_M)])
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юлокхрсдмшу яоейрпнб МЮФЛХРЕ <ENTER>')
pause
figure('Name','Amplitude Spectrum','NumberTitle', 'off')
subplot(2,1,1), stem(k,MOD_N,'MarkerSize',3), grid, xlabel('k')
title(strcat(['Amplitude Spectrum of the Periodic Sequence N = ',num2str(N)]))
subplot(2,1,2), stem(k1,MOD_M,'MarkerSize',3), grid, xlabel('k')
title(strcat(['Amplitude Spectrum of the Periodic Sequence M = ',num2str(M)]))
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.3. хяякеднбюмхе бнглнфмнярх слемэьемхъ пюярейюмхъ яоейрпю я онлныэч нймю')
win_M  = hamming(M)';    % нймн ущллхмцю ≈ бейрнп-ярнкаеж дкхмш M
xw_M = x_M.*win_M;       % онякеднбюрекэмнярэ, бгбеьеммюъ нймнл
XW_M = fft(xw_M);        % дот бгбеьеммни онякеднбюрекэмнярх
MODW_M =(2/M)*abs(XW_M); % юлокхрсдмши яоейрп бгбеьеммни онякеднбюрекэмнярх
MODW_M(1) =(1/M)*abs(XW_M(1));
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юлокхрсдмшу яоейрпнб дн Х оняке ОПХЛЕМЕМХЪ нймю МЮФЛХРЕ <ENTER>')
pause
figure('Name','Reducing Spectrum Leakage with the help of Window Functions','NumberTitle', 'off')
subplot(2,1,1), stem(k1,MOD_M,'MarkerSize',3), grid, xlabel('k')
title(strcat(['Amplitude spectrum without windowing M = ',num2str(M)]))
subplot(2,1,2), stem(k1,MODW_M,'MarkerSize',3), grid, xlabel('k')
title(strcat(['Amplitude spectrum with Hamming Window M = ',num2str(M)]))
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.4. хяякеднбюмхе щттейрю пюярейюмхъ яоейрпю дкъ ясллш дбсу дхяйпермшу цюплнмхй')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше')
DATA=0;
while DATA==0
f1_1 = input('f1_1 = ');   % вюярнрш дхяйпермшу цюплнмхй (цЖ)
f2_1 = input('f2_1 = ');
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
n = 0:(N-1);               % дхяйпермне мнплхпнбюммне бпелъ
k = 0:(N-1);               % дхяйпермюъ мнплхпнбюммюъ вюярнрю
w1_1 = 2*pi*f1_1/Fs; w2_1 = 2*pi*f2_1/Fs; % мнплхпнбюммше вюярнрш дхяйпермшу цюплнмхй (пюд)
x1 = A1*cos(w1_1*n)+A2*cos(w2_1*n);  % онякеднбюрекэмнярэ (оепхнд N)
X1 = fft(x1);           % дот онякеднбюрекэмнярх (оепхнд N)
MOD1 = (2/N)*abs(X1);   % юлокхрсдмши яоейрп онякеднбюрекэмнярх
MOD1(1) = (1/N)*abs(X1(1));
P1_1 = N*f1_1/Fs;        % вхякн оепхнднб дхяйпермни цюплнмхйх я вюярнрни f1_1 мю оепхнде онякеднбюрекэмнярх N
P2_1 = N*f2_1/Fs;        % вхякн оепхнднб дхяйпермни цюплнмхйх я вюярнрни f2_1 мю оепхнде онякеднбюрекэмнярх N
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ вхякю оепхнднб ДХЯЙПЕРМШУ ЦЮПЛНМХЙ я вюярнрюлх f1_1 Х f2_1 МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp(['      f1_1 = ',num2str(f1_1),'  -->  P1_1 = ' num2str(P1_1)])
disp(['      f2_1 = ',num2str(f2_1),'  -->  P2_1 = ' num2str(P2_1)])
win_N  = hamming(N)';     % нймн ущллхмцю ≈ бейрнп-ярнкаеж дкхмш N
xw1 = x1.*win_N;   % онякеднбюрекэмнярэ, бгбеьеммюъ нймнл (оепхнд N)
XW1 = fft(xw1);    % дот бгбеьеммни онякеднбюрекэмнярх (оепхнд N)
MODW1 =(2/N)*abs(XW1);   % юлокхрсдмши яоейрп онякеднбюрекэмнярх
MODW1(1) =(1/M)*abs(XW1(1));
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юлокхрсдмшу яоейрпнб дн Х оняке ОПХЛЕМЕМХЪ нймю МЮФЛХРЕ <ENTER>')
pause
figure('Name','Reducing Spectrum Leakage with the help of Window Functions','NumberTitle', 'off')
subplot(2,1,1), stem(k,MOD1,'MarkerSize',3), grid, xlabel('k')
title(strcat(['Amplitude spectrum without windowing N = ',num2str(N)]))
subplot(2,1,2), stem(k,MODW1,'MarkerSize',3), grid, xlabel('k')
title(strcat(['Amplitude spectrum with Hamming Window N = ',num2str(N)]))
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.5. сксвьемхе пюгкхвемхъ дхяйпермшу цюплнмхй я акхгйн пюяонкнфеммшлх вюярнрюлх')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше')
DATA=0;
while DATA==0
f1_2 = input('f1_2 = ');        % вюярнрш дхяйпермшу цюплнмхй (цЖ)
f2_2 = input('f2_2 = ');
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ оепхндю онякеднбюрекэмнярх Х')
disp('% вюярнр цюплнмхй МЮФЛХРЕ <ENTER>')
disp('%')
disp('%')
disp(['      N = ',num2str(N)])
disp(['      f1_2 = ',num2str(f1_2),'       f2_2 = ' num2str(f2_2)])
Delta_N = Fs/N;                   % пюгпеьемхе он вюярнре
Delta_f = abs(f1_2-f2_2);         % пюяярнъмхе лефдс вюярнрюлх
L = ceil(Fs/(Delta_f-Delta_N));   % бшапюммюъ дкхмю L
Delta_L = Fs/L;        % оепхнд дхяйперхгюжхх он вюярнре опх дкхме L
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ пюгпеьемхъ он вюярнре Delta_N,')
disp('% пюяярнъмхъ ЛЕФДС вюярнрюлх Delta_f,')
disp('% дкхмш L ОНЯКЕДНБЮРЕКЭМНЯРХ')
disp('% Х оепхндю дхяйперхгюжхх он вюярнре Delta_L МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp(['      Delta_N = ',num2str(Delta_N)])
disp(['      Delta_f = ',num2str(Delta_f)])
disp(['      L = ',num2str(L)])
disp(['      Delta_L = ',num2str(Delta_L)])
disp('%')
disp('%')
n = 0:(N-1);                              % дхяйпермне мнплхпнбюммне бпелъ
w1_2 = 2*pi*f1_2/Fs; w2_2 = 2*pi*f2_2/Fs; % мнплхпнбюммше вюярнрш
x2 = A1*cos(w1_2*n)+A2*cos(w2_2*n);       % йнмевмюъ онякеднбюрекэмнярэ
X2 = fft(x2);              % дот йнмевмни онякеднбюрекэмнярх дкхмш N
MOD2 = abs(X2);            % лндскэ дот
X2_L = fft(x2,L);          % дот йнмевмни онякеднбюрекэмнярх, днонкмеммни мскълх дн дкхмш L
MOD2_L = abs(X2_L);        % лндскэ дот
disp('% дКЪ БШБНДЮ цпютхйнб N-рнвевмнцн дот Х лндскъ яоейрпюкэмни')
disp('% окнрмнярх, бняярюмнбкеммни он L рнвйюл, МЮФЛХРЕ <ENTER>')
pause
k = 0:(N-1);          % дхяйпермюъ мнплхпнбюммюъ вюярнрю опх дкхме N
k1 = 0:(L-1);         % дхяйпермюъ мнплхпнбюммюъ вюярнрю опх дкхме L
figure('Name','Discrete Harmonic Signal with Close Frequencies','NumberTitle', 'off')
subplot(2,1,1), stem(k,MOD2), grid, xlabel('k')
title(strcat(['DFT Modulus N = ',num2str(N)]))
subplot(2,1,2), plot(k1,MOD2_L,'r','MarkerSize',3, 'Linewidth',2)
grid, hold on, stem(k1,MOD2_L,':'), xlabel('k')
title(strcat(['Spectral Density Modulus L = ',num2str(L)]))
L_2 = ceil(L/2);                     % нямнбмюъ онкняю вюярнр L/2
[MODm m]= max(MOD2_L(1:(L_2)));      % люйяхлсл MODm х хмдейя m бейрнпю MOD2_L (оепбши охй)
k_1 = (m-1); f_1 = k_1*Delta_L;      % дхяйпермюъ мнплхпнбюммюъ х юаянкчрмюъ (цЖ) вюярнрш оепбнцн охйю
K = ceil(L/N);   % йнкхвеярбн нрявернб мю оепхнде дхяйперхгюжхх Fs/N
K1 = m+K; K2 = m+2*K-1;              % мхфмъъ K1 Х бепумъъ K2 цпюмхжш хмрепбюкю опх онхяйе брнпнцн охйю яопюбю
[MODm1 m1]= max(MOD2_L(K1:K2));      % люйяхлсл MODm1 х хмдейя m1 лндскъ дот MOD2_L мю хмрепбюке [K1 K2]
K3 = m-(2*K-1); K4 = m-K;            % мхфмъъ K3 Х бепумъъ K4 цпюмхжш хмрепбюкю опх онхяйе брнпнцн охйю якебю
[MODm2 m2]= max(MOD2_L(K3:K4));      % люйяхлсл MODm2 х хмдейя m2 лндскъ дот MOD2_L мю хмрепбюке [K3 K4]
if (MODm1>MODm2)
    k_2 = (K1+m1-1)-1; f_2 = k_2*Delta_L; % дхяйпермюъ мнплхпнбюммюъ х юаянкчрмюъ (цЖ) вюярнрш брнпнцн охйю, еякх нм яопюбю нр оепбнцн
    else
    k_2 = (K3+m2-1)-1; f_2 = k_2*Delta_L; % дхяйпермюъ мнплхпнбюммюъ х юаянкчрмюъ (цЖ) вюярнрш брнпнцн охйю, еякх нм якебю нр оепбнцн
end
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ вюярнр цюплнмхй МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp(['      k_1 = ',num2str(k_1),'      f_1 = ' num2str(f_1)])
disp(['      k_2 = ',num2str(k_2),'      f_2 = ' num2str(f_2)])
disp('%')
disp('%')
disp('% нОПЕДЕКХРЕ вюярнрш цюплнмхй ОН цпютхйс')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.6. бшвхякемхе йпсцнбни ябепрйх')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше')
DATA=0;
while DATA==0
x3 = input('x3 = ');       % оепбюъ онякеднбюрекэмнярэ
x4 = input('x4 = ');       % брнпюъ онякеднбюрекэмнярэ
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
y34 = ifft(fft(x3).*fft(x4));% йпсцнбюъ ябепрйю онякеднбюрекэмняреи
L34 = length(y34);           % оепхнд йпсцнбни ябепрйх
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ЦПЮТХЙНБ онякеднбюрекэмняреи Х йпсцнбни ЯБЕПРЙХ (3 ОЕПХНДЮ) МЮФЛХРЕ <ENTER>')
pause
figure('Name','Sequences x3, x4, y34','NumberTitle', 'off')
subplot(3,1,1), stem((0:3*L34-1),...
repmat(x3,1,3),'fill','Linewidth',2,'MarkerSize',3), grid
xlabel('n'), title('Periodic Sequence x3(n)')
subplot(3,1,2), stem((0:3*L34-1), repmat(x4,1,3),'fill', 'Linewidth',2,'MarkerSize',3), grid
xlabel('n'), title('Periodic Sequence x4(n)')
subplot(3,1,3), stem((0:3*L34-1), repmat(y34,1,3),'fill', 'Linewidth',2,'MarkerSize',3), grid, xlabel('n')
title('Periodic Sequence y34(n) ≈ Convolution with FFT and IFFT')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.7. бшвхякемхе кхмеимни ябепрйх')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше')
DATA=0;
while DATA==0
x5 = input('x5 = ');       % оепбюъ онякеднбюрекэмнярэ
x6 = input('x6 = ');       % брнпюъ онякеднбюрекэмнярэ
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
y56_1 = conv(x5,x6);       % кхмеимюъ ябепрйю,бшвхякеммюъ я онлныэч тсмйжхх conv
y56_2 = fftfilt(x5,x6);    % кхмеимюъ ябепрйю, бшвхякеммюъ я онлныэч тсмйжхх fftfilt
MAX = max([length(y56_1) length(y56_2)]); % люйяхлюкэмюъ дкхмю ябепрйх
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб онякеднбюрекэмняреи Х кхмеимни ЯБЕПРЙХ МЮФЛХРЕ <ENTER>')
pause
figure('Name','Sequences x5, x6, y56_1, y56_2','NumberTitle', 'off')
subplot(4,1,1)
stem((0:length(x5)-1),x5,'fill','Linewidth',2,'MarkerSize',3)
grid, xlabel('n'), title('Sequence x5(n)'), xlim([0 MAX-1])
subplot(4,1,2)
stem((0:length(x6)-1),x6,'fill','Linewidth',2,'MarkerSize',3)
grid, xlabel('n'), title('Sequence x6(n)'), xlim([0 MAX-1])
subplot(4,1,3)
stem((0:length(y56_1)-1),y56_1,'fill','Linewidth',2,'MarkerSize',3)
grid, xlabel('n'), title('Sequence y56(n) ≈ Convolution'), xlim([0 MAX-1])
subplot(4,1,4)
stem((0:length(y56_2)-1),y56_2,'fill','Linewidth',2,'MarkerSize',3)
grid, xlabel('n'), title('Sequence y56(n) ≈ Convolution with FFT and IFFT'), xlim([0 MAX-1])
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.8. бшвхякемхе пеюйжхх кдя он тнплске ябепрйх')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше')
DATA=0;
while DATA==0
b = input('b = ');        % йнщттхжхемрш вхякхрекъ оепедюрнвмни тсмйжхх
a = input('a = ');        % йнщттхжхемрш гмюлемюрекъ оепедюрнвмни тсмйжхх
N1 = input('N1 = ');      % дкхмю хлоскэямни уюпюйрепхярхйх
N2 = input('N2 = ');      % дкхмю бнгдеиярбхъ
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
h = impz(b,a,N1)';       % хлоскэямюъ уюпюйрепхярхйю
x7 = input_1(N2);        % бнгдеиярбхе
y7_1 = conv(x7,h);       % пеюйжхъ, бшвхякеммюъ я онлныэч тсмйжхх conv
y7_2 = fftfilt(h,x7);    % пеюйжхъ, бшвхякеммюъ я онлныэч тсмйжхх fftfilt
L=N1+N2-1;            % дкхмю ябепрйх, бшвхякеммни я онлныэч тсмйжхх conv
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ЦПЮТХЙНБ ху, бнгдеиярбхъ Х пеюйжхх МЮФЛХРЕ <ENTER>')
pause
figure('Name','Impulse Response, Input and Output Signals','NumberTitle', 'off')
subplot(4,1,1)
stem(0:length(h)-1,h,'Linewidth',2,'MarkerSize',3), grid,
xlabel('n'), title('Impulse Response  h(n)'), xlim([0 L-1])
subplot(4,1,2)
stem(0:length(x7)-1,x7,'Linewidth',2,'MarkerSize',3), grid
xlabel('n'), title('Input Signal x7(n)'), xlim([0 L-1])
subplot(4,1,3)
stem(0:length(y7_1)-1,y7_1,'Linewidth',2,'MarkerSize',3),
grid
xlabel('n'), title('Output Signal y7(n) ≈ Convolution'), xlim([0 L-1])
subplot(4,1,4)
stem(0:length(y7_2)-1,y7_2,'Linewidth',2,'MarkerSize',3), grid
xlabel('n'), title('Output Signal y7(n) ≈ Convolution with FFT and IFFT'), xlim([0 L-1])
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.9. бшвхякемхе пеюйжхх кдя лернднл оепейпшрхъ я мюйнокемхел')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше')
DATA=0;
while DATA==0
N3 = input('N3 = ');   % дкхмю бнгдеиярбхъ
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
x8 = input_1(N3);        % бнгдеиярбхе
y8_1 = fftfilt(h,x8);    % пеюйжхъ, бшвхякеммюъ я онлныэч тсмйжхх fftfilt
y8_2 = fftfilt(h,x8,N1); % пеюйжхъ, бшвхякеммюъ я онлныэч тсмйжхх fftfilt лернднл мюйнокемхъ я оепейпшрхел
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб бнгдеиярбхъ Х пеюйжхх МЮФЛХРЕ <ENTER>')
pause
figure('Name','Impulse Response, Input and Output Signals ≈ Overlap-add method','NumberTitle', 'off')
subplot(4,1,1)
stem(0:length(h)-1,h,'MarkerSize',3), grid
xlabel('n'), title('Impulse Response h(n)'), xlim([0 N3-1])
subplot(4,1,2), stem(0:length(x8)-1,x8,'MarkerSize',3), grid
xlabel('n'), title('Input Signal x8(n)')
subplot(4,1,3),stem(0:length(y8_1)-1,y8_1,'MarkerSize',3), grid
xlabel('n')
title('Output Signal y8(n) ≈ Convolution with FFT and IFFT')
subplot(4,1,4), stem(0:length(y8_2)-1,y8_2,'MarkerSize',3), grid
xlabel('n')
title('Output Signal y8(n) ≈ Convolution with Overlap-add method')
disp('%')
disp('%')
disp('% пюанрю гюбепьемю')






























































































