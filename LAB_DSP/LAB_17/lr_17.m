script
clc
clear
disp('% кп ╧17. оюпюлерпхвеяйхи яоейрпюкэмши юмюкхг')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше')
DATA=0;
while DATA==0
Nb = input('Nb = ');      % мнлеп апхцюдш
L = input('L = ');        % дкхмю онякеднбюрекэмнярх
Fs = input('Fs = ');      % вюярнрю дхяйперхгюжхх (цЖ)
a = input('a = ');        % бейрнп оюпюлерпнб юп-лндекх
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
disp('% О.1. лндекхпнбюмхе яксвюимни онякеднбюрекэмнярх мю нямнбе юп-лндекх')
n = 0:(L-1);           % дхяйпермне мнплхпнбюммне бпелъ
e = randn(1,L);        % бнгдеиярбхе - мнплюкэмши аекши ьсл дкхмш L
y = filter(1,a,e);     % лндекхпселюъ яксвюимюъ онякеднбюрекэмнярэ дкхмш L
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйю лндекхпселни ОНЯКЕДНБЮРЕКЭМНЯРХ МЮФЛХРЕ <ENTER>')
pause
figure('Name','AR-sequence and True PSD','NumberTitle', 'off')
subplot(2,1,1), plot(n,y), grid, xlabel('n'), ylabel('y(n)')
title('AR-sequence')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.2. бшвхякемхе хярхммни яол лндекхпселни онякеднбюрекэмнярх')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйю хярхммни яол ЛНДЕКХПСЕЛНИ ОНЯКЕДНБЮРЕКЭМНЯРХ МЮФЛХРЕ <ENTER>')
pause
f = 0:Fs/(L-1):Fs;          % бейрнп вюярнр (цЖ)
H = freqz(1,a,f,Fs);        % йнлокейямюъ вюярнрмюъ уюпюйрепхярхйю
S = (1/Fs)*abs(H).^2;       % хярхммюъ яол лндекхпселни онякеднбюрекэмнярх
subplot(2,1,2),plot(f,S), grid, xlabel('f (Hz)'), ylabel('S(f)')
title('True PSD')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.3. нжемйю норхлюкэмнцн онпъдйю юп-лндекх юмюкхгхпселни онякеднбюрекэмнярх')
disp('%')
disp('%')
x = y;      % юмюкхгхпселюъ онякеднбюрекэмнярэ онкюцюеряъ пюбмни лндекхпселни
pmin = 1; pmax = 3*(length(a)-1);  % лхмхлюкэмне х люйяхлюкэмне гмювемхъ онпъдйю
p = [pmin:pmax];                   % бейрнп гмювемхи онпъдйю юп-лндекх
for i = pmin:pmax
[aYW D] = aryule(x,p(i));          % оюпюлерпш юп-лндекх х япедмхи йбюдпюр ньхайх кхмеимнцн опедяйюгюмхъ
BIC(i) = L*log(D)+p(i)*log(L);     % гмювемхе йпхрепхъ аюиеяю
variance(i) = D;                   % япедмхи йбюдпюр ньхайх кхмеимнцн опедяйюгюмхъ ≈ нжемйю дхяоепяхх мнплюкэмнцн аекнцн ьслю б юп-лндекх
end
[BIC_min p_opt] = min(BIC);        % лхмхлюкэмне гмювемхе йпхрепхъ аюиеяю х норхлюкэмши онпъднй лндекх
disp('% дКЪ БШБНДЮ цпютхйнб ГЮБХЯХЛНЯРЕИ япедмецн йбюдпюрю ньхайх кхмеимнцн опедяйюгюмхъ')
disp('% Х гмювемхи йпхрепхъ аюиеяю НР онпъдйю лндекх МЮФЛХРЕ <ENTER>')
pause
figure('Name','Mean Square of the Linear Prediction Error and Bayesian Information Criterion','NumberTitle','off')
subplot(2,1,1), plot(p,variance,'Linewidth',2), grid
xlabel('p'), ylabel('D'), title('Mean Square of the Linear Prediction Error')
subplot(2,1,2), plot(p,BIC,'r','Linewidth',2), grid
xlabel('p'), ylabel('BIC')
title('Bayesian Information Criterion')
disp('%')
disp('%')
disp('% бБЕДХРЕ ГМЮВЕМХЕ норхлюкэмнцн онпъдйю юп-ЛНДЕКХ')
DATA=0;
while DATA==0
disp('%')
p_opt = input('  p_opt = ');        % норхлюкэмши онпъднй юп-лндекх
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
disp('% О.4. бшвхякемхе нжемнй оюпюлерпнб юп-лндекх')
disp('%')
disp('%')
[aYW DYW] = aryule(x,p_opt);     % нжемйю оюпюлерпнб юп-лндекх он лерндс чкю-снкйепю
[aBURG DBURG] = arburg(x,p_opt); % нжемйю оюпюлерпнб юп-лндекх он лерндс аепцю
[aCOV DCOV] = arcov(x,p_opt);    % нжемйю оюпюлерпнб юп-лндекх йнбюпхюжхнммшл лернднл
[aMCOV DMCOV] = armcov(x,p_opt); % нжемйю оюпюлерпнб юп-лндекх лндхтхжхпнбюммшл йнбюпхюжхнммшл лернднл
D = var(e);       % нжемйю дхяоепяхх мнплюкэмнцн аекнцн ьслю дкхмш L
disp('% дКЪ БШБНДЮ хярхммшу оюпюлерпнб юп-лндекх Х ХУ нжемнй МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp(' хярхммше оюпюлерпш');
disp(['  a =    [' num2str(a),']'])
disp('%')
disp(' лернд чкю-снкйепю');
disp('%')
disp([' aYW =   [' num2str(aYW),']']);
disp('%')
disp(' лернд аепцю');
disp('%')
disp([' aBURG = [' num2str(aBURG),']']);
disp('%')
disp(' йнбюпхюжхнммши лернд');
disp('%')
disp([' aCOV =  [' num2str(aCOV),']']);
disp('%')
disp(' лндхтхжхпнбюммши йнбюпхюжхнммши лернд');
disp('%')
disp([' aMCOV = [' num2str(aMCOV) ,']']);
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ япедмецн йбюдпюрю ньханй кхмеимнцн опедяйюгюмхъ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp([' лернд чкю-снкйепю:                       DYW = ' num2str(DYW)]);
disp('%')
disp([' лернд аепцю:                             DBURG = ' num2str(DBURG)]);
disp('%')
disp([' йнбюпхюжхнммши лернд:                    DCOV = ' num2str(DCOV)]);
disp('%')
disp([' лндхтхжхпнбюммши йнбюпхюжхнммши лернд:   DMCOV = ' num2str(DMCOV)]);
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ нжемйх дхяоепяхх мнплюкэмнцн аекнцн ьслю МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp([' нжемйю дхяоепяхх мнплюкэмнцн аекнцн ьслю: D = ' num2str(D)]);
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.5. опнбепйю сярнивхбнярх аху-тхкэрпю')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ йюпр мскеи х онкчянб МЮФЛХРЕ <ENTER>')
pause
figure('Name','Z-plane zero-pole plots','NumberTitle','off')
subplot(2,2,1), zplane(1,aYW), title ('Yule-Walker method'), grid
xlabel('Re'), ylabel('jIm')
subplot(2,2,2), zplane(1,aBURG), title ('Burg method'), grid
xlabel('Re'), ylabel('jIm')
subplot(2,2,3), zplane(1,aCOV), title ('Covariance method'),grid
xlabel('Re'), ylabel('jIm')
subplot(2,2,4), zplane(1,aMCOV),title ('Modified Covariance method')
xlabel('Re'), ylabel('jIm'), grid
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.6. бшвхякемхе нжемнй яол')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб нжемнй яол МЮФЛХРЕ <ENTER>')
pause
[SYW,f] = pyulear(x,p_opt,L,Fs,'twosided');  % нжемйю яол он лерндс чкю-снкйепю
[SBURG,f] =  pburg(x,p_opt,L,Fs,'twosided'); % нжемйю яол он лерндс аепцю
[SCOV,f] =  pcov(x,p_opt,L,Fs,'twosided');   % нжемйю яол йнбюпхюжхнммшл лернднл
[SMCOV,f] =  pmcov(x,p_opt,L,Fs,'twosided'); % нжемйю яол лндхтхжхпнбюммшл йнбюпхюжхнммшл лернднл
figure('Name','PSD Estimates','NumberTitle', 'off')
subplot(4,1,1), plot(f,SYW,'Linewidth',2), grid
xlabel('f (Hz)'), ylabel('SYW(f)'), title(' PSD estimate using Yule-Walker method')
subplot(4,1,2), plot(f,SBURG,'Linewidth',2), grid
xlabel('f (Hz)'), ylabel('SBURG(f)'), title('PSD estimate using Burg method')
subplot(4,1,3), plot(f,SCOV,'Linewidth',2), grid
xlabel('f (Hz)'), ylabel('SCOV(f)'), title('PSD estimate using Covariance method')
subplot(4,1,4), plot(f,SMCOV,'Linewidth',2), grid
xlabel('f (Hz)'), ylabel('SMCOV(f)')
title('PSD estimate using Modified Covariance method')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.7. япюбмемхе нжемнй яол C хярхммни яол')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб нжемнй яол Х хярхммни яол МЮФЛХРЕ <ENTER>')
pause
figure('Name','True PSD and Different PSD Estimates', 'NumberTitle','off')
plot(f,S,'Linewidth',2), xlabel('f (Hz)'), ylabel('S(f)'), grid
hold on
plot(f,SYW,'m','Linewidth',2), grid
plot(f,SBURG,'r','Linewidth',2), grid
plot(f,SCOV,'k','Linewidth',2), grid
plot(f,SMCOV,'g','Linewidth',2), grid
legend('True PSD','PSD estimate Yule-Walker method', 'PSD estimate Burg method', 'PSD estimate Covariance method', 'PSD estimate Modified Covariance method', 0);
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.8. бшвхякемхе гмювемхи RMSE')
RMSE_YW = sqrt((1/L).*sum((S'-SYW).^2));     % RMSE дкъ нжемйх яол он лерндс чкю-снкйепю
RMSE_BURG = sqrt((1/L).*sum((S'-SBURG).^2)); % RMSE дкъ нжемйх яол он лерндс аепцю
RMSE_COV = sqrt((1/L).*sum((S'-SCOV).^2));   % RMSE дкъ нжемйх яол йнбюпхюжхнммшл лернднл
RMSE_MCOV = sqrt((1/L).*sum((S'-SMCOV).^2)); % RMSE дкъ нжемйх яол лндхтхжхпнбюммшл йнбюпхюжхнммшл лернднл
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ гмювемхи RMSE нжемнй яол МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp([' лернд чкю-снкйепю:                     RMSE = ' num2str(RMSE_YW )]);
disp('%')
disp([' лернд аепцю:                           RMSE = ' num2str(RMSE_BURG )]);
disp('%')
disp([' йнбюпхюжхнммши лернд:                  RMSE = ' num2str(RMSE_COV)]);
disp('%')
disp([' лндхтхжхпнбюммши йнбюпхюжхнммши лернд: RMSE = ' num2str(RMSE_MCOV)]);
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.9. хяякеднбюмхе бкхъмхъ онпъдйю юп-лндекх мю нжемйс яол')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб хярхммни яол Х нжемнй яол ОН ЛЕРНДС чкю-снкйепю МЮФЛХРЕ <ENTER>')
pause
p_low = round(p_opt/2);                   % гюмхфеммши онпъднй юп-лндекх
p_high = 3*p_opt;                         % гюбшьеммши онпъднй юп-лндекх
[SYW_low,f] = pyulear(x,p_low,L,Fs,'twosided');   % нжемйю яол он лерндс чкю-снкйепю я гюмхфеммшл онпъдйнл лндекх
[SYW_high,f] = pyulear(x,p_high,L,Fs,'twosided'); % нжемйю яол он лерндс чкю-снкйепю я гюбшьеммшл онпъдйнл лндекх
figure('Name','Different AR Model Orders for Yule-Walker Method', 'NumberTitle','off')
plot(f,S,'Linewidth',2),  xlabel('f (Hz)'), grid, hold on
plot(f,SYW_low,'r','Linewidth',2), grid
plot(f,SYW_high,'k','Linewidth',2), grid
legend(['True PSD: p opt = ' num2str(p_opt)], ['PSD estimate Yule-Walker method: p = ' num2str(p_low)], ['PSD estimate Yule-Walker method: p = ' num2str(p_high)], 0);
title('PSD estimates using Yule-Walker method for different model orders')
disp('%')
disp('%')
disp('% О.10. хяякеднбюмхе бкхъмхъ дкхмш онякеднбюрекэмнярх мю нжемйс яол')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб хярхммни яол Х нжемнй яол ОН ЛЕРНДС чкю-снкйепю МЮФЛХРЕ <ENTER>')
pause
e1 = randn(1,100*L);           % бнгдеиярбхе - мнплюкэмши аекши ьсл дкхмш 100*L
y1 = filter(1,a,e1);           % лндекхпселюъ яксвюимюъ онякеднбюрекэмнярэ дкхмш 100*L
x1 = y1;                       % юмюкхгхпселюъ онякеднбюрекэмнярэ онкюцюеряъ пюбмни лндекхпселни
[SYW_1,f] = pyulear (x,p_opt,L,Fs,'twosided');   % нжемйю яол он лерндс чкю-снкйепю дкъ онякеднбюреккэмнярх дкхмш L
[SYW_2,f] = pyulear (x1,p_opt,L,Fs,'twosided');  % нжемйю яол он лерндс чкю-снкйепю дкъ онякеднбюреккэмнярх дкхмш 100*L
figure('Name','Different Sequence Lengths for Yule-Walker Method', 'NumberTitle','off')
plot(f,S,'Linewidth',2), xlabel('f (Hz)'), grid, hold on
plot(f,SYW_1,'r','Linewidth',2), grid
plot(f,SYW_2,'k','Linewidth',2), grid
legend(['True PSD: L = ' num2str(L)], [' PSD estimate Yule-Walker method: L = ' num2str(L)], [' PSD estimate Yule-Walker method: L = ' num2str(100*L)], 0);
title('PSD estimates using Yule-Walker method for different lengths')
disp('%')
disp('%')
disp('% пюанрю гюбепьемю')






















































































