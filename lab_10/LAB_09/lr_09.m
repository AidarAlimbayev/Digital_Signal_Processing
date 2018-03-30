script
clc
clear
disp('% кп ╧9. дхяйпермне опенапюгнбюмхе тспэе (ВЮЯРЭ 1)')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше');
DATA=0;
while DATA==0
Nb = input('Nb = ');           % мнлеп апхцюдш
N = input('N = ');             % дкхмю (оепхнд) онякеднбюрекэмнярх
Fs = input('Fs = ');           % вюярнрю дхяйперхгюжхх
T = input('T = ');             % оепхнд дхяйперхгюжхх 1/Fs
A1 = input('A1 = ');           % юлокхрсдш дхяйпермшу цюплнмхй
A2 = input('A2 = ');
f1 = input('f1 = ');           % вюярнрш (цЖ) дхяйпермшу цюплнмхй
f2 = input('f2 = ');
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
disp('%')
disp('% дКЪ БШБНДЮ хяундмшу юлокхрсд Х вюярнр дхяйпермшу цюплнмхй МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp(['      A1 = ' num2str(A1) '      A2 = ' num2str(A2)])
disp(['      f1 = ' num2str(f1) '      f2 = ' num2str(f2)])
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.1. бшвхякемхе юлокхрсдмнцн х тюгнбнцн яоейрпнб оепхндхвеяйни онякеднбюрекэмнярх')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб ОЕПХНДХВЕЯЙНИ ОНЯКЕДНБЮРЕКЭМНЯРХ МЮФЛХРЕ <ENTER>')
pause
n = 0:(N-1);                 	 % дхяйпермне мнплхпнбюммне бпелъ
k = 0:(N-1);                      % дхяйпермюъ мнплхпнбюммюъ вюярнрю
w1 = 2*pi*f1/Fs; w2 = 2*pi*f2/Fs;  % мнплхпнбюммше вюярнрш дхяйпермшу цюплнмхй (пюд)
x = A1*cos(w1*n+pi/4)+A2*cos(w2*n+pi/8);   % оепхндхвеяйюъ онякеднбюрекэмнярэ
X = fft(x);                        % дот оепхндхвеяйни онякеднбюрекэмнярх
MOD = (2/N)*abs(X);                % юлокхрсдмши яоейрп оепхндхвеяйни онякеднбюрекэмнярх
MOD(1) = (1/N)*abs(X(1));
PHASE = angle(X);                  % тюгнбши яоейрп оепхндхвеяйни онякеднбюрекэмнярх
for i = 1:N
    if (abs(X(i)) < 1e-4)
       PHASE(i)=0;
    end
end
figure('Name','Periodic Sequence','NumberTitle','off')
subplot(3,1,1), stem(n,x, 'MarkerSize',3,'Linewidth',2)
grid, xlabel('n')
ylabel('x(n)'), title(strcat(['Periodic Sequence x(n)  N = ',num2str(N)]))
subplot(3,1,2), stem(n/Fs,x,'MarkerSize',3,'Linewidth',2)
grid, xlabel('nT')
ylabel('x(nT)'), title(strcat(['Periodic Sequence x(nT)  N = ',num2str(N)]))
x = ifft(X);                     % оепхндхвеяйюъ онякеднбюрекэмнярэ, бшвхякеммюъ я онлныэч ндот
subplot(3,1,3), stem(n,x,'MarkerSize',3,'Linewidth',2)
grid, xlabel('n')
ylabel('x(n)'), title(strcat(['Periodic Sequence x = ifft(X)  N = ',num2str(N)]))
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юлокхрсдмнцн яоейрпю ОЕПХНДХВЕЯЙНИ ОНЯКЕДНБЮРЕКЭМНЯРХ МЮФЛХРЕ <ENTER>')
pause
figure('Name','Amplitude Spectrum','NumberTitle', 'off')
subplot(2,1,1), stem(k,MOD,'MarkerSize',3,'Linewidth',2), grid
xlabel('k'), ylabel('1/N|X(k)|')
title(strcat(['Amplitude Spectrum of the Periodic Sequence N = ',num2str(N)]))
subplot(2,1,2), stem(k*(Fs/N),MOD,'MarkerSize',3,'Linewidth',2),grid
xlabel('f (Hz)'), ylabel('1/N|X(f)|')
title(strcat(['Amplitude Spectrum of the Periodic Sequence N = ',num2str(N)]))
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб тюгнбнцн яоейрпю ОЕПХНДХВЕЯЙНИ ОНЯКЕДНБЮРЕКЭМНЯРХ МЮФЛХРЕ <ENTER>')
pause
figure('Name','Phase Spectrum','NumberTitle', 'off')
subplot(2,1,1), stem(k, PHASE,'MarkerSize',3,'Linewidth',2), grid
xlabel('k'), ylabel('arg{X(k)} (rad)')
title(strcat(['Phase Spectrum of the Periodic Sequence N = ',num2str(N)]))
subplot(2,1,2), stem(k*(Fs/N),PHASE,'MarkerSize',3,'Linewidth',2)
grid, xlabel('f (Hz)'), ylabel('arg{X(f)} (rad)')
title(strcat(['Phase Spectrum of the Periodic Sequence N = ',num2str(N)]))
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.2. бшвхякемхе дот йнмевмни онякеднбюрекэмнярх')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб лндскъ дот ЙНМЕВМНИ ОНЯКЕДНБЮРЕКЭМНЯРХ Х юлокхрсдмнцн яоейрпю')
disp('% ОЕПХНДХВЕЯЙНИ ОНЯКЕДНБЮРЕКЭМНЯРХ МЮФЛХРЕ <ENTER>')
pause
MOD_K = abs(fft(x));        % лндскэ дот йнмевмни онякеднбюрекэмнярх
figure('Name','DFT Modulus and Amplitude Spectrum', 'NumberTitle','off')
subplot(2,1,1), stem(k,MOD_K,'MarkerSize',3,'Linewidth',2), grid
xlabel('k'), ylabel('|X(k)|')
title('DFT Modulus of the Finite Sequence')
subplot(2,1,2), stem(k,MOD,'MarkerSize',3,'Linewidth',2), grid
xlabel('k'), ylabel('1/N |X(k)|')
title('Amplitude Spectrum of the Periodic Sequence')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.3. нопедекемхе юлокхрсд х вюярнр дхяйпермшу цюплнмхй')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ бшундмшу оюпюлерпнб тсмйжхх fft_e1 МЮФЛХРЕ <ENTER>')
pause
e1 = 1e-7;                    % гмювемхе онпнцю дкъ оепбнцн йпхрепхъ
[MODm,m] = fft_e1(MOD,e1)     % бмеьмъъ тсмйжхъ дкъ бшдекемхъ юлокхрсд х вюярнр цюплнмхй онкегмнцн яхцмюкю он оепбнлс йпхрепхч
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ юлокхрсд Х вюярнр дхяйпермшу цюплнмхй МЮФЛХРЕ <ENTER>')
pause
A1 = MODm(1); A2 = MODm(2);    % юлокхрсдш дхяйпермшу цюплнмхй
k1 = m(1);  k2 = m(2);         % дхяйпермше мнплхпнбюммше вюярнрш
f1 = k1*Fs/N; f2 = k2*Fs/N;    % вюярнрш (цЖ) дхяйпермшу цюплнмхй
disp('%')
disp('%')
disp(['      A1 = ' num2str(A1) '      A2 = ' num2str(A2)])
disp(['      k1 = ' num2str(k1) '      k2 = ' num2str(k2)])
disp(['      f1 = ' num2str(f1) '      f2 = ' num2str(f2)])
disp('%')
disp('%')
disp('% япюбмхре Я бшундмшлх оюпюлерпюлх ТСМЙЖХХ fft_e1 Х ХЯУНДМШЛХ ДЮММШЛХ')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.4. цпюмхвмше гмювемхъ онпнцнб дкъ оепбнцн х брнпнцн йпхрепхеб бшдекемхъ онкегмнцн яхцмюкю')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ЦПЮМХВМШУ ГМЮВЕМХИ ОНПНЦЮ ДКЪ оепбнцн йпхрепхъ МЮФЛХРЕ <ENTER>')
pause
noise = randn(1,N);             % мнплюкэмши аекши ьсл
s = x+noise;                    % юддхрхбмюъ ялеяэ яхцмюкю я ьслнл
S = fft(s);                     % дот ялеях яхцмюкю я ьслнл
MODS = (2/N)*abs(S);            % юлокхрсдмши яоейрп ялеях яхцмюкю я ьслнл
MODS(1) = (1/N)*abs(S(1));
NOISE = fft(noise);             % дот ьслю
MODNOISE = (2/N)*abs(NOISE);    % юлокхрсдмши яоейрп ьслю
MODNOISE(1) = (1/N)*abs(NOISE(1));
MAX_NOISE = max(MODNOISE);     % люйяхлсл юлокхрсдмнцн яоейрпю ьслю
MAXS = max(MODS);              % люйяхлсл юлокхрсдмнцн яоейрпю ялеях яхцмюкю я ьслнл
e1_low = MAX_NOISE/MAXS; % мхфмъъ цпюмхжю онпнцю дкъ оепбнцн йпхрепхъ
e1_up = 1;               % бепумъъ цпюмхжю онпнцю дкъ оепбнцн йпхрепхъ
P = (1/N)*sum(MODS.^2);  % япедмъъ лнымнярэ ялеях яхцмюкю я ьслнл
MAXS2 = MAXS.^2;         % йбюдпюр люйяхлслю юлокхрсдмнцн яоейрпю ялеях яхцмюкю я ьслнл
MAX_NOISE2 = MAX_NOISE.^2;  % йбюдпюр люйяхлслю юлокхрсдмнцн яоейрпю ьслю
disp('%')
disp('%')
disp(['  e1_low = ' num2str(e1_low) '  e1_up = ' num2str(e1_up)])
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ЦПЮМХВМШУ ГМЮВЕМХИ ОНПНЦЮ ДКЪ брнпнцн йпхрепхъ МЮФЛХРЕ <ENTER>')
pause
e2_low = MAX_NOISE2/P;     % мхфмъъ цпюмхжю онпнцю дкъ брнпнцн йпхрепхъ
e2_up = MAXS2/P;           % бепумъъ цпюмхжю онпнцю дкъ брнпнцн йпхрепхъ
disp('%')
disp('%')
disp(['   e2_low = ' num2str(e2_low) '   e2_up = ' num2str(e2_up)])
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.5. бшдекемхе онкегмнцн яхцмюкю он оепбнлс йпхрепхч')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйю ЮДДХРХБМНИ ЯЛЕЯХ ЯХЦМЮКЮ Я ЬСЛНЛ МЮФЛХРЕ <ENTER>')
pause
figure('Name','Mixture of Signal and Noise','NumberTitle', 'off')
stem(n,s,'MarkerSize',3,'Linewidth',2), grid
xlabel('n'), ylabel('s(n)')
title(strcat(['Mixture of Signal and Noise  N = ',num2str(N)]))
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб ЮЛОКХРСДМНЦН Х мнплхпнбюммнцн ЮЛОКХРСДМНЦН ЯОЕЙРПНБ')
disp('% ЮДДХРХБМНИ ЯЛЕЯХ ЯХЦМЮКЮ Я ЬСЛНЛ МЮФЛХРЕ <ENTER>')
pause
figure('Name','Amplitude Spectrum and Normalized Amplitude Spectrum','NumberTitle', 'off')
subplot(2,1,1), stem(k,MODS,'MarkerSize',3,'Linewidth',2), grid
xlabel('k'), ylabel('|S(k)|')
title(strcat(['Amplitude Spectrum N = ',num2str(N)]))
subplot(2,1,2), stem(k, MODS/MAXS,'MarkerSize',3,'Linewidth',2)
grid, xlabel('k'), ylabel('|S(k)|/max|S(k)|')
title(strcat(['Normalized Amplitude Spectrum N = ',num2str(N)]))
disp('%')
disp('%')
disp('% бБЕДХРЕ БШАПЮММНЕ ГМЮВЕМХЕ ОНПНЦЮ e1 ДКЪ оепбнцн йпхрепхъ')
disp('%')
e1 = input('    e1 = ');  % бшапюммне гмювемхе онпнцю дкъ оепбнцн йпхрепхъ
disp('%')
disp('% дКЪ БШБНДЮ бшундмшу оюпюлерпнб тсмйжхх fft_e1 МЮФЛХРЕ <ENTER>')
pause
[MODm,m] = fft_e1(MODS,e1) % бмеьмъъ тсмйжхъ дкъ бшдекемхъ юлокхрсд х вюярнр цюплнмхй онкегмнцн яхцмюкю он оепбнлс йпхрепхч
disp('%')
disp('%')
disp('% япюбмхре ГМЮВЕМХЪ бшдекеммшу он оепбнлс йпхрепхч юлокхрсд х вюярнр')
disp('% Я ХЯУНДМШЛХ ДЮММШЛХ')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.6. бшдекемхе онкегмнцн яхцмюкю он брнпнлс йпхрепхч')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб ЮЛОКХРСДМНЦН ЯОЕЙРПЮ Х йбюдпюрю ЮЛОКХРСДМНЦН')
disp('% ЯОЕЙРПЮ, мнплхпнбюммнцн Й БЕКХВХМЕ ЯПЕДМЕИ ЛНЫМНЯРХ')
disp('% ЮДДХРХБМНИ ЯЛЕЯХ ЯХЦМЮКЮ Я ЬСЛНЛ, МЮФЛХРЕ <ENTER>')
pause
figure('Name','Amplitude Spectrum and Normalized Amplitude Spectrum Squire','NumberTitle', 'off')
subplot(2,1,1), stem(k,MODS,'MarkerSize',3,'Linewidth',2), grid
xlabel('k'), ylabel('|S(k)|')
title(strcat(['Amplitude Spectrum N = ',num2str(N)]))
subplot(2,1,2), stem(k,(MODS.^2)/P,'MarkerSize',3,'Linewidth',2)
grid, xlabel('k'), ylabel('|S(k)|^2/P')
title(strcat(['Normalized Amplitude Spectrum Squire N = ',num2str(N)]))
disp('%')
disp('%')
disp('% бБЕДХРЕ БШАПЮММНЕ ГМЮВЕМХЕ ОНПНЦЮ e2 ДКЪ брнпнцн йпхрепхъ')
disp('%')
e2 = input('    e2 = ');  % бшапюммне гмювемхе онпнцю дкъ брнпнцн йпхрепхъ
disp('%')
disp('% дКЪ БШБНДЮ бшундмшу оюпюлерпнб тсмйжхх fft_e2 МЮФЛХРЕ <ENTER> ')
pause
[MODm,m] = fft_e2(MODS,e2)% бмеьмъъ тсмйжхъ дкъ бшдекемхъ юлокхрсд х вюярнр цюплнмхй онкегмнцн яхцмюкю он брнпнлс йпхрепхч
disp('%')
disp('%')
disp('% япюбмхре ГМЮВЕМХЪ бшдекеммшу он брнпнлс йпхрепхч юлокхрсд х вюярнр')
disp('% Я ХЯУНДМШЛХ ДЮММШЛХ')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.7. бняярюмнбкемхе юмюкнцнбнцн яхцмюкю')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб онякеднбюрекэмнярх Х лндскъ ЕЕ дот,')
disp('% бняярюмнбкеммнцн юмюкнцнбнцн яхцмюкю Х ЕЦН яоейрпю')
disp('% Х хяундмнцн юмюкнцнбнцн яхцмюкю МЮФЛХРЕ <ENTER>')
pause
Xa  = [X(N/2+1:N),X(1:N/2)]; % яоейрп юмюкнцнбнцн яхцмюкю (я рнвмнярэч дн онярнъммнцн лмнфхрекъ)
i = 1;                       % явервхй гмювемхи юмюкнцнбнцн яхцмюкю
for t = 0:0.25*T:(N-1)*T     % гмювемхъ меопепшбмнцн бпелемх
    s = 0;        
    for k = -N/2:N/2-1       % дхяйпермюъ мнплхпнбюммюъ вюярнрю
        s = s + Xa(k+N/2+1)*exp(j*2*pi*k*t/(N*T)); % бняярюмнбкемхе юмюкнцнбнцн яхцмюкю
    end
    xa(i) = (1/N).*s;        % гмювемхъ бняярюмнбкеммнцн юмюкнцнбнцн яхцмюкю
    i = i+1;
end
t = 0:0.25*T:(N-1)*T;
xt = A1*cos(2*pi*f1*t+pi/4)+A2*cos(2*pi*f2*t+pi/8);    % гмювемхъ хяундмнцн юмюкнцнбнцн яхцмюкю
k = 0:N-1;                    % дхяйпермюъ мнплхпнбюммюъ вюярнрю
MODa = (2/N)*abs(Xa);         % юлокхрсдмши яоейрп бняярюмнбкеммнцн юмюкнцнбнцн яхцмюкю
MODa(1) = (1/N)*abs(Xa(1));
figure('Name','Original Periodic Sequence & FFT, Reconstructed Analog Signal & Spectrum, Original Analog Signal','NumberTitle', 'off')
subplot(3,2,1), stem(n,x,'MarkerSize',3), grid
xlabel('n'), ylabel('x(n)')
title(strcat(['Original Periodic Sequence  N = ',num2str(N)]))
subplot(3,2,2), stem(k,abs(X),'MarkerSize',3,'Linewidth',2), grid
xlabel('k'), ylabel('|X(k)|')
title(strcat(['DFT of Original Periodic Sequence  N = ',num2str(N)]))
subplot(3,2,3), plot(t,real(xa)), grid, xlabel('t')
ylabel('x(t)'),title('Reconstructed Analog Signal')
k = -N/2:N/2-1;
subplot(3,2,4), stem(k,MODa,'MarkerSize',3,'Linewidth',2), grid
xlabel('k'), ylabel('|Xa(k)|')
title('Amplitude Spectrum of Reconstructed Analog Signal')
subplot(3,2,5), plot(t,xt), grid, xlabel('t')
ylabel('x(t)'), title('Original Analog Signal')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.8. бняярюмнбкемхе яоейрпюкэмни окнрмнярх йнмевмни онякеднбюрекэмнярх')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб дот Х яоейрпюкэмни окнрмнярх ЙНМЕВМНИ')
disp('% ОНЯКЕДНБЮРЕКЭМНЯРХ, БШВХЯКЕММНИ дбслъ ЯОНЯНАЮЛХ, МЮФЛХРЕ <ENTER>')
pause
L = 2*N;          % йнкхвеярбн нрявернб яоейрпюкэмни окнрмнярх мю оепхнде
l = 0;           
for l = 0:(L-1)         % дхяйпермюъ мнплхпнбюммюъ вюярнрю
    S = 0;    
    for n = 0:(N-1)     % дхяйпермне мнплхпнбюммне бпелъ
        S = S + x(n+1)*exp(-j*2*pi*l*n/L);     % бняярюмнбкемхе яоейрпюкэмни окнрмнярх
    end
    XW(l+1) = S;        % гмювемхъ бняярюмнбкеммни яоейрпюкэмни окнрмнярх
    l = l+1;
end
xz = [x zeros(1,(L-N))];      % онякеднбюрекэмнярэ, днонкмеммюъ мскълх дн дкхмш L
XZ = fft(xz);                 % дот онякеднбюрекэмнярх, днонкмеммни мскълх
k = 0:(N-1);                  % дхяйпермюъ мнплхпнбюммюъ вюярнрю
w = 0:2*pi/L:2*pi-2*pi/L;     % мнплхпнбюммюъ вюярнрю
l = 0:(L-1);                  % дхяйпермюъ мнплхпнбюммюъ вюярнрю
figure('Name','DFT and Spectral Density','NumberTitle', 'off')
subplot(3,1,1), stem(k,abs(X),'MarkerSize',3,'Linewidth',2)
grid, xlabel('k'), ylabel('|X(k)')
title(strcat(['DFT Modulus N = ',num2str(N)]))
subplot(3,1,2), plot(w,abs(XW),'MarkerSize',3,'Linewidth',2)
grid, xlabel('w'), ylabel('|X(w)|')
title(strcat(['Spectral Density Modulus (option 1) L = ',num2str(L)]))
subplot(3,1,3), plot(w,abs(XZ),'MarkerSize',3,'Linewidth',2)
grid, xlabel('w'), ylabel('|X(w)|')
title(strcat(['Spectral Density Modulus (option 2) L = ',num2str(L)]))
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.9. слемэьемхе оепхндю дхяйперхгюжхх он вюярнре опх бшвхякемхх дот')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб йнмевмшу онякеднбюрекэмняреи,')
disp('% дот Х яоейрпюкэмшу окнрмняреи МЮФЛХРЕ <ENTER>')
pause
figure('Name','Finite Sequences, DFT and Spectral Densities','NumberTitle', 'off')
L = [N 2*N 4*N];
for i = 1:length(L)
    xz = [x zeros(1,(L(i)-N))];   % онякеднбюрекэмнярэ, днонкмеммюъ мскълх дн дкхмш L(i)
    XZ = fft(xz);
    Delta_f(i) = Fs/L(i);
    n = 0:length(xz)-1;           % дхяйпермне мнплхпнбюммне бпелъ
    k = 0:length(XZ)-1;           % дхяйпермюъ мнплхпнбюммюъ вюярнрю
subplot(3,2,2*i-1), stem(n,xz,'MarkerSize',3), xlabel('n'), grid
title(strcat(['Finite Sequence x(n) L = ',num2str(L(i))]))
subplot(3,2,2*i), plot(k,abs(XZ), 'r','MarkerSize',3, 'Linewidth',2), grid, hold on, stem(k,abs(XZ),':'), xlabel('k')
title(strcat(['DFT and Spectral Density Modulus L = ',num2str(L(i))]))
end
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ оепхнднб дот Х оепхнднб дхяйперхгюжхх он вюярнре МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp(['      L = [',num2str(L) ']'])
disp('%')
disp(['      Delta_f = [',num2str(Delta_f) ']'])
disp('%')
disp('%')
disp('% пюанрю гюбепьемю')






























































































