script
clc
clear
disp('% кп ╧16. меоюпюлерпхвеяйхи яоейрпюкэмши юмюкхг')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше')
DATA=0;
while DATA==0
Nb = input('Nb = ');      % мнлеп апхцюдш
N = input('N = ');        % дкхмю онякеднбюрекэмнярх
Fs = input('Fs = ');      % вюярнрю дхяйперхгюжхх (цЖ)
A1 = input('A1 = ');      % юлокхрсдш дхяйпермшу цюплнмхй
A2 = input('A2 = ');   
f1 = input('f1 = ');      % вюярнрш дхяйпермшу цюплнмхй (цЖ)
f2 = input('f2 = ');    
sigma = input('sigma = ');% бейрнп гмювемхи яйн ьслю
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
disp('% О.1. опнбепйю хмтнплюрхбмнярх оепхнднцпюллш б гюбхяхлнярх нр спнбмъ ьслю')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб онякеднбюрекэмняреи Х оепхнднцпюлл')
disp('% Я ПЮГКХВМШЛХ яйн ЬСЛЮ МЮФЛХРЕ <ENTER>')
pause
w1 = 2*pi*f1/Fs; w2 = 2*pi*f2/Fs;  % мнплхпнбюммше вюярнрш дхяйпермшу цюплнмхй (пюд)
n = 0:(N-1);                       % дхяйпермне мнплхпнбюммне бпелъ
x = A1*cos(w1*n')+A2*cos(w2*n');   % онякеднбюрекэмнярэ ≈ ясллю дбсу цюплнмхй (бейрнп-ярнкаеж)
figure('Name',' Harmonic Signals Embedded in White Gaussian Noise and Periodograms','NumberTitle', 'off')
for i = 1:length(sigma)         % хмдейя щкелемнб бейрнпю sigma
xe = x'+sigma(i).*randn(1,N);   % онякеднбюрекэмнярх я гюдюммшлх яйн ьслю (x'-бейрнп-ярпнйю)
subplot(4,2,2*i-1), plot(n,xe,'Linewidth',2), grid
xlabel('n'), ylabel(strcat('xe',num2str(i),'(n)'))
title(strcat(['Sequence ',num2str(i), '   STD = ',num2str(sigma(i))]))
[Se,f] = periodogram(xe,[],N,Fs,'twosided'); % оепхнднцпюллю онякеднбюрекэмнярх я гюдюммшл яйн ьслю
subplot(4,2,2*i), plot(f,Se,'Linewidth',2), grid
xlabel('f (Hz)'), ylabel(strcat('S',num2str(i),'(f)'))
title(strcat(['Periodogram ',num2str(i),' (Vt/Hz)   STD = ',num2str(sigma(i))]))
end
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.2. опнбепйю хмтнплюрхбмнярх оепхнднцпюллш б гюбхяхлнярх')
disp('% нр оепхндю дхяйперхгюжхх он вюярнре')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб оепхнднцпюлл ОПХ ПЮГКХВМНИ ПЮГЛЕПМНЯРХ дот МЮФЛХРЕ <ENTER>')
pause
figure('Name','Periodograms of Harmonic Signal with defferent DFT size','NumberTitle', 'off')
xe = x'+randn(1,N);         % онякеднбюрекэмнярэ дкхмш N
M = [N/16 N/8 N 8*N];       % бейрнп пюглепмняреи дот
for i = 1:length(M)
[Se,f] = periodogram(xe,[],M(i),Fs,'onesided');   % оепхнднцпюллю б нямнбмни онкняе вюярнр опх гюдюммни пюглепмнярх дот  
subplot(4,1,i), plot(f, Se,'Linewidth',2),grid
xlabel('f (Hz)'), ylabel(strcat('S',num2str(i),'(f)'))
title(strcat(['Periodogram ',num2str(i),' (Vt/Hz)   length DFT=', num2str(M(i))]))
end
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.3. опнбепйю нжемйх яол мю юяхлорнрхвеяйсч меялеыеммнярэ х янярнърекэмнярэ')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб ГЮБХЯХЛНЯРХ ялеыемхъ Х япедмецн йбюдпюрю')
disp('% нрйкнмемхъ хярхммни яол нр ее нжемйх МЮФЛХРЕ <ENTER>')
pause
SWN = var(randn(1,100000))/Fs;        % хярхммюъ яол мнплюкэмнцн аекнцн ьслю
N_WN = 1000:1000:100000;              % бейрнп дкхм ьслю
for i = 1:length(N_WN)                % хмдейя щкелемнб бейрнпю N_WN
e = randn(1,N_WN(i));                 % ьсл гюдюммни дкхмш
SWN_estimate = periodogram(e,[],N_WN(i),Fs,'twosided');   % нжемйю яол ьслю
beta(i) = mean(SWN-SWN_estimate');    % ялеыемххе нжемйх яол
mean_square (i) = var(SWN_estimate)+beta(i)^2;            % япедмхи йбюдпюр нрйкнмемхъ хярхммни яол нр ее нжемйх
end
figure('Name','Bias and Mean square Deviation of true PSD from its Estimate','NumberTitle', 'off')
subplot(2,1,1), plot(N_WN, beta,'LineWidth',2), grid, xlabel('N')
ylabel('beta')
title('Bias')
subplot(2,1,2), plot(N_WN, mean_square,'LineWidth',2), grid, xlabel('N')
ylabel('meanerr')
title('Mean square Deviation of true PSD from its Estimate')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.4. тнплхпнбюмхе яксвюимни онякеднбюрекэмнярх я рпеаселни юйт')
N = 1000;                        % дкхмю ьслю
xw = randn(1,N);                 % мнплюкэмши аекши ьсл дкхмш N
N02 = var(xw)+(mean(xw)).^2;     % йнмярюмрю N0/2
m = -(N-1):(N-1);                % дхяйпермне мнплхпнбюммне бпелъ дкъ юйт, жемрпхпнбюммни нрмняхрекэмн m = 0
Ry_required = 0.25.*0.95.^abs(m); % рпеаселюъ юйт
L = 2*N-1;                        % дкхмю юйт
m = 0:L-1;                        % дхяйпермне мнплхпнбюммне бпелъ дкъ юйт, жемрпхпнбюммни нрмняхрекэмн m = N
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйю рпеаселни юйт МЮФЛХРЕ <ENTER>')
pause
figure('Name','Required ACF','NumberTitle', 'off')
stem(m,Ry_required), grid, xlabel('m'), title('Required ACF Ry')
disp('%')
disp('%')
disp('% я ОНЛНЫЭЧ ЙМНОЙХ Zoom in НОПЕДЕКХРЕ Х ББЕДХРЕ вермши онпъднй йху-тхкэрпю')
DATA=0;
while DATA==0
disp('%')
R = input('    R = ');            % вермши онпъднй йху-тхкэрпю
disp('%')
disp('% оПНБЕПЭРЕ опюбхкэмнярэ ББНДЮ хяундмшу дюммшу')
disp('% оПХ опюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 1')
disp('% оПХ меопюбхкэмшу хяундмшу дюммшу ББЕДХРЕ 0 Х онбрнпхре ББНД')
DATA = input('--> ');
end
Sy = 2*real(fft(Ry_required(N:L),L))-Ry_required (N); % яол б L рнвйюу, бшвхякеммюъ он рпеаселни юйт
A = sqrt(real(Sy)./N02);          % юву йху-тхкэрпю б L рнвйюу
k = 0:L-1;                        % дхяйпермюъ мнплхпнбюммюъ вюярнрю
F = -k*pi*R/L;                    % ктву йху-тхкэрпю б L рнвйюу
j = sqrt(-1);                  
H = A.*exp(j*F);                  % ву йху-тхкэрпю б L рнвйюу
h1 = real(ifft(H));               % ху йху-тхкэрпю мю оепхнде L йпсцнбни ябепрйх
y_ACF = fftfilt(h1,xw);           % яксвюимюъ онякеднбюрекэмнярэ дкхмш N я рпеаселни юйт мю бшунде йху-тхкэрпю
h = h1(1:R+1);                    % ху дкхмш R+1
Ry_estimate = xcorr(y_ACF)./N;    % нжемйю юйт пеюйжхх йху-тхкэрпю, жемрпхпнбюммюъ нрмняхрекэмн N
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйю рпеаселни юйт х нжемйх юйт пеюйжхх МЮФЛХРЕ <ENTER>')
pause
figure('Name','Required ACF and ACF Estimate of Output Signal','NumberTitle', 'off')
subplot(2,1,1), stem(m,Ry_required), grid, title('Required ACF Ry')
subplot(2,1,2), stem(m,Ry_estimate), grid, xlabel('m')
title(' ACF Estimate of Output Signal - Ry estimate')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб ху, мнплюкэмнцн аекнцн ьслю х пеюйжхх йху-ТХКЭРПЮ МЮФЛХРЕ <ENTER>')
pause
n = 0:(N-1);  % дхяйпермне мнплхпнбюммне бпелъ дкъ бнгдеиярбхъ х пеюйжхх йху-тхкэрпю
figure('Name','Impulse Response, Input and Output Signals','NumberTitle', 'off')
subplot(3,1,1), stem(0:R,h), grid, title('Impulse Response h(n)')
subplot(3,1,2), plot(n,xw), grid
title('Input Signal - White Gaussian Noise')
subplot(3,1,3), plot(n,y_ACF), grid, xlabel('n')
title('Output Signal with Required ACF - y ACF')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.5. тхкэрпюжхъ яксвюимни онякеднбюрекэмнярх я рпеаселни юйт')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб хяундмни Х опенапюгнбюммни ОНЯКЕДНБЮРЕКЭМНЯРЕИ МЮФЛХРЕ <ENTER>')
pause
[b,a] = butter(3,0.3,'high');      % йнщттхжхемрш аху-тхкэрпю тбв аюррепбнпрю
y = filter(b,a,y_ACF);             % онякеднбюрекэмнярэ мю бшунде аху-тхкэрпю оняке сдюкемхъ рпемдю
figure('Name','Input and Output Signals of Butterworth filter','NumberTitle', 'off')
subplot(2,1,1), plot(n,y_ACF), grid, title('Input Signal of Butterworth filter ≈ y ACF')
subplot(2,1,2), plot(n,y), grid, xlabel('n'), ylabel('y(n)')
title('Output Signal of Butterworth filter - y')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.6. пюявер оепхнднцпюллш')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйю оепхнднцпюллш МЮФЛХРЕ <ENTER>')
pause
[S,f] = periodogram(y,[],N,Fs,'twosided');         % оепхнднцпюллю яксвюимни онякеднбюрекэмнярх
figure('Name','Periodogram of the Non-white Gaussian Noise','NumberTitle', 'off')
subplot(2,1,1), plot(f,S,'Linewidth',2), grid
xlabel('f (Hz)'), ylabel('S(f) (Vt/Hz)')
title('Periodogram of the Non-white Gaussian noise')
subplot(2,1,2), periodogram(y,[],N,Fs,'twosided')
xlabel('f (Hz)'), ylabel('S(f) (dB/Hz)')
title('Periodogram of the Non-white Gaussian noise')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.7. пюявер оепхнднцпюллш дюмэеккю')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб оепхнднцпюлл дюммэеккю МЮФЛХРЕ <ENTER>')
pause
K = [5 10 20];                   % бейрнп йнкхвеярбю сяпедмъелшу вюярнр
figure('Name','Daniell Periodograms for the Different Number of Frequency Intervals','NumberTitle', 'off')
for i = 1:3                      % хмдейя щкелемнб бейрнпю K
S1 = [S(N-K(i)+1:N); S; S(1:K(i))];  % оепхнднцпюллю, оепхндхвеяйх опнднкфеммюъ якебю х яопюбю мю K нрявернб
S2 = smooth(S1,K(i));            % пегскэрюр бшвхякемхъ яйнкэгъыецн япедмецн (бейрнп дкхмш N+2K(i))
SD(:,i) = S2(K(i)+1:N+K(i));     % оепхнднцпюллю дюммэеккю (бейрнп дкхмш N) дкъ йнкхвеярбю сяпедмъелшу вюярнр K(i)
subplot(4,1,i+1), plot(f, SD(:,i),'Linewidth',2), grid
xlabel('f (Hz)'), ylabel(strcat('SD',num2str(i),'(f)'))
title(strcat(['Daniell Periodogram ',num2str(i),'  Frequency Interval 2K+1, K=',num2str(K(i))]))
end
subplot(4,1,1)
plot(f,S,'Linewidth',2), grid, xlabel('f (Hz)'), ylabel('S(f)')
title('Original non-modified periodogram (Vt/Hz)')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.8. пюявер оепхнднцпюллш аюпркеррю')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб оепхнднцпюллш аюпркеррю МЮФЛХРЕ <ENTER>')
pause
L = [10 20 40];                         % бейрнп дкхм тпюцлемрнб
figure('Name',' Bartlett Periodograms for Different Fragment Lengths','NumberTitle', 'off')
for i = 1:3                             % хмдейя щкелемнб бейрнпю L
SB(:,i) = pwelch(y,rectwin(L(i)),0,N,Fs,'twosided'); % оепхнднцпюллю аюпркеррю дкъ тпюцлемрю дкхмш L(i)
subplot(4,1,i+1), plot(f, SB(:,i),'Linewidth',2), grid
xlabel('f (Hz)'), ylabel(strcat('SB',num2str(i),'(f)'))
title(strcat([' Bartlett Periodogram ',num2str(i),'  L =',num2str(L(i))]))
end
subplot(4,1,1)
plot(f,S,'Linewidth',2), grid, xlabel('f (Hz)'), ylabel('S(f)')
title('Original non-modified periodogram (Vt/Hz)')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.9. пюявер оепхнднцпюллш сщквю')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб оепхнднцпюллш сщквю МЮФЛХРЕ <ENTER>')
pause
num = 1;                      % онпъдйнбши мнлеп оепхнднцпюллш сщквю
L = ceil([0.1 0.05 0.025].*length(y));      % бейрнп дкхм тпюцлемрнб
figure('Name','Welch Periodograms for Different Fragment Lengths and Overlapping','NumberTitle', 'off')
for i = 1:3                              % хмдейя щкелемнб бейрнпю L
Q = ceil(0.0125*length(y));              % бекхвхмю оепейпшрхъ
SW(:,num) = pwelch(y,L(i),Q,N,Fs,'twosided'); % оепхнднцпюллю сщквю опх дкхме тпюцлемрю L(i) х бекхвхме оепейпшрхъ Q
subplot(4,1,i+1), plot(f,SW(:,num),'Linewidth',2), grid
xlabel('f (Hz)'), ylabel('S(f)')
title(strcat(['Welch Periodogram Fragment length L = ',num2str(L(i)), '  Overlapping Q = ',num2str(Q)]))
num = num+1;
subplot(4,1,1), plot(f,S,'Linewidth',2), grid, xlabel('f (Hz)')
ylabel('S(f)'), title('Original non-modified periodogram (Vt/Hz)')
end
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.10. пюявер нжемйх яол он лерндс акщйлюмю-рэчйх')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб нжемйх яол, БШВХЯКЕММШУ ОН ЛЕРНДС акщйлюмю-рэчйх, МЮФЛХРЕ <ENTER>')
pause
N1 = ceil(N/10);               % люйяхлюкэмши ядбхц он бпелемх дкъ нжемйх юйт нртхкэрпнбюммнцн ьслю
R1 = (1/N).*xcorr(y,N1);       % нжемйю юйт нртхкэрпнбюммнцн ьслю дкхмш 2N1+1, жемрпхпнбюммюъ нрмняхрекэмн N1+1
R = R1(2:length(R1)-1);        % нжемйю юйт нртхкэрпнбюммнцн ьслю дкхмш 2N1-1, жемрпхпнбюммюъ нрмняхрекэмн N1
L1 = length(R);                % дкхмю нжемйх юйт
Rw(:,1) = R'.*rectwin(L1);     % юйт, бгбеьеммюъ опълнсцнкэмшл нймнл
Rw(:,2) = R'.*hamming(L1);     % юйт, бгбеьеммюъ нймнл ущллхмцю
Rw(:,3) = R'.*chebwin(L1);     % юйт, бгбеьеммюъ нймнл веашьебю
name(1).win = 'Rectangular Window';   % хлемю нйнм (люяяхб гюохяеи name я ндмхл онкел win)
name(2).win = 'Hamming Window';
name(3).win = 'Chebyshev Window';
f = 0:Fs/N:Fs-Fs/N;             % бейрнп вюярнр (цЖ)
figure('Name','PSD estimates by the Blackman-Tukey method','NumberTitle', 'off')
for i = 1:3                     % мнлепю нйнм
SBT(:,i) = (1/Fs)*(2*real(fft(Rw(N1:L1,i),N))- Rw(N1,i));  % нжемйю яол дкхмш N, бшвхякеммюъ он юйт, бгбеьеммни нймнл
subplot(4,1,i+1), plot(f,SBT(:,i),'Linewidth',2), grid
xlabel('f (Hz)'), ylabel(strcat('SBT',num2str(i),'(f)'))
title(['PSD Estimate  -  ',strcat(name(i).win)])
end
subplot(4,1,1)
plot(f,S,'Linewidth',2), grid, xlabel('f (Hz)'), ylabel('S(f)')
title('Original non-modified periodogram (Vt/Hz)')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.11. нопедекемхе онйюгюрекеи йювеярбю нжемнй яол')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ гмювемхи яйн МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('   яйн ОЕПХНДНЦПЮЛЛШ')
format long                   % дкхммши тнплюр дкъ бшбндю яйн
STD_S = std(S)                % яйн оепхнднцпюллш
disp('   яйн ОЕПХНДНЦПЮЛЛ дюмэеккю ОПХ ПЮГКХВМНЛ ЙНКХВЕЯРБЕ СЯПЕДМЪЕЛШУ ВЮЯРНР K')
STD_SD = [K' std(SD)']        % йнкхвеярбн сяпедмъелшу вюярнр х яйн оепхнднцпюлл дюмэеккю
disp('   яйн ОЕПХНДНЦПЮЛЛ аюпркеррю ОПХ ПЮГКХВМНИ ДКХМЕ ТПЮЦЛЕМРЮ L')
STD_SB = [L' std(SB)']        % дкхмш тпюцлемрнб х яйн оепхнднцпюлл аюпркеррю
disp('   яйн ОЕПХНДНЦПЮЛЛШ сщквю ОПХ ПЮГКХВМНИ ДКХМЕ ТПЮЦЛЕМРЮ L Х БЕКХВХМЕ ОЕПЕЙПШРХЪ Q')
LL = [L(1) L(2) L(3)];        % дкхмш тпюцлемрнб
Q = [Q Q Q];                  % бекхвхмю оепейпшрхъ
STD_SW = [LL' Q' std(SW)']    % дкхмш тпюцлемрнб, бекхвхмю оепейпшрхъ х яйн оепхнднцпюлл сщквю
disp('   яйн НЖЕМНЙ яол ОН ЛЕРНДС акщйлюмю-рэчйх ОПХ ПЮГКХВМШУ НЙМЮУ')
WINDOW = {name.win};          % хлемю нйнм (люяяхб ъвеей ≈ cell array)
STD_SBT = std(SBT);           % яйн нжемнй яол он лерндс акщйлюмю-рэчйх
STD_SBT = [WINDOW(1)' STD_SBT(1)'; WINDOW(2)' STD_SBT(2)'; WINDOW(3)' STD_SBT(3)']    % хлемю нйнм х яйн нжемнй яол он лерндс акщйлюмю-рэчйх
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ гмювемхи днапнрмнярх МЮФЛХРЕ <ENTER>')
pause
disp('%')
format          % бнгпюр й хяундмнлс тнплюрс дкъ бшбндю днапнрмняреи
disp('   днапнрмнярэ ОЕПХНДНЦПЮЛЛШ')
Q_S = mean(S).^2/var(S)               % днапнрмнярэ оепхнднцпюллш
disp('   днапнрмнярэ ОЕПХНДНЦПЮЛЛ дюмэеккю ОПХ ПЮГКХВМНЛ ЙНКХВЕЯРБЕ СЯПЕДМЪЕЛШУ ВЮЯРНР K')
Q_SD = mean(SD).^2./var(SD);          % днапнрмнярэ оепхнднцпюлл дюмэеккю
Q_SD = [K' Q_SD']
disp('   днапнрмнярэ ОЕПХНДНЦПЮЛЛ аюпркеррю ОПХ ПЮГКХВМНИ ДКХМЕ ТПЮЦЛЕМРЮ L')
Q_SB = mean(SB).^2./var(SB);          % днапнрмнярэ оепхнднцпюлл аюпркеррю
Q_SB = [L' Q_SB']
disp('   днапнрмнярэ ОЕПХНДНЦПЮЛЛ сщквю ОПХ ПЮГКХВМНИ ДКХМЕ ТПЮЦЛЕМРЮ L Х БЕКХВХМЕ ОЕПЕЙПШРХЪ Q')
Q_SW = mean(SW).^2./var(SW);          % днапнрмнярэ оепхнднцпюлл сщквю
Q_SW = [LL' Q' Q_SW']
disp('   днапнрмнярэ НЖЕМНЙ яол ОН ЛЕРНДС акщйлюмю-рэчйх ОПХ ПЮГКХВМШУ НЙМЮУ')
Q_SBT = real(mean(SBT).^2./var(SBT)); % днапнрмнярэ нжемнй яол он лерндс акщйлюмю-рэчйх
Q_SBT = [WINDOW(1)' Q_SBT(1)'; WINDOW(2)' Q_SBT(2)'; WINDOW(3)' Q_SBT(3)']
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.12. онярпнемхе яоейрпнцпюллш')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйю яоейрпнцпюллш ДХЯЙПЕРМНЦН ЦЮПЛНМХВЕЯЙНЦН ЯХЦМЮКЮ МЮФЛХРЕ <ENTER>')
pause
N = 4000;                          % дкхмю онякеднбюрекэмнярх
n = 0:N-1;                         % дхяйпермне мнплхпнбюммне бпелъ
x = A1*cos(w1*n')+A2*cos(w2*n');   % онякеднбюрекэмнярэ ≈ ясллю дбсу цюплнмхй
figure('Name','Harmonic Signal Spectrogram','NumberTitle', 'off')
spectrogram(x,128,120,128,Fs,'yaxis')
colorbar
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Harmonic Signal Spectrogram')
disp('%')
disp('%')
disp('% пюанрю гюбепьемю')








































































































