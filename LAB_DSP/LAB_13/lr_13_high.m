script
clc
clear
disp('% кп ╧13. яхмрег аху-тхкэрпю тбв лернднл ахкхмеимнцн Z-опенапюгнбюмхъ')
disp('%')
disp('%')
disp('% О.1. ббнд рпеанбюмхи й юву (Да) тбв')
disp('%')
disp('%')
disp('% бБЕДХРЕ мнлеп апхцюдш Х рпеанбюмхъ Й юву (Да)')
DATA=0;
while DATA==0;
Nb = input('Nb = ');         % мнлеп апхцюдш
Fs = input('Fs = ');         % вюярнрю дхяйперхгюжхх (цЖ)
fk = input('fk = ');         % цпюмхвмюъ вюярнрю ог (цЖ)
ft = input('ft = ');         % цпюмхвмюъ вюярнрю оо (цЖ)
rs = input('rs = ');         % лхмхлюкэмн дносярхлне гюрсуюмхе б ог
rp = input('rp = ');         % люйяхлюкэмн дносярхлне гюрсуюмхе б оо
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
disp('% О.2. яхмрег аху-тхкэрпнб аюррепбнпрю, веашьебю I х II пндю х гнкнрюпебю-йюсщпю')
disp('%')
disp('%')
disp('% дКЪ ЯХМРЕГЮ аху-ТХКЭРПНБ тбв МЮФЛХРЕ <ENTER>')
pause
WDp = ft/(Fs/2); WDs = fk/(Fs/2);    % цпюмхвмше мнплхпнбюммше вюярнрш оо Х ог
[R1,WDn1] = buttord(WDp,WDs,rp,rs);  % онпъднй х вюярнрю япегю аху-тхкэрпю тбв аюррепбнпрю
[R2,WDn2] = cheb1ord(WDp,WDs,rp,rs); % онпъднй х вюярнрю япегю аху-тхкэрпю тбв веашьебю I пндю
[R3,WDn3] = cheb2ord(WDp,WDs,rp,rs); % онпъднй х вюярнрю япегю аху-тхкэрпю тбв веашьебю веашьебю II пндю
[R4,WDn4] = ellipord(WDp,WDs,rp,rs); % онпъднй х вюярнрю япегю аху-тхкэрпю тбв гнкнрюпебю-йюсщпю
[b1,a1] = butter(R1,WDn1,'high');    % йнщттхжхемрш аху-тхкэрпю тбв аюррепбнпрю
[b2,a2] = cheby1(R2,rp,WDn2,'high'); % йнщттхжхемрш аху-тхкэрпю тбв веашьебю I пндю
[b3,a3] = cheby2(R3,rs,WDn3,'high'); % йнщттхжхемрш аху-тхкэрпю тбв веашьебю II пндю
[b4,a4] = cheby2(R4,rs,WDn4,'high'); % йнщттхжхемрш аху-тхкэрпю тбв веашьебю гнкнрюпебю-йюсщпю
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ОНПЪДЙНБ аху-ТХКЭРПНБ тбв МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp(['   R1 = ' num2str(R1),'   R2 = ' num2str(R2),'   R3 = ' num2str(R3),'   R4 = ' num2str(R4)])
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.3. юмюкхг уюпюйрепхярхй аху-тхкэрпнб тбв')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ уюпюйрепхярхй аху-тхкэрпнб тбв (вершпе цпютхвеяйху нймю) МЮФЛХРЕ <ENTER>')
pause
figure('Name','Highpass IIR Filter Butterworth','NumberTitle', 'off')
plot_iir(b1,a1,Fs)     % уюпюйрепхярхйх аху-тхкэрпю тбв аюррепбнпрю
figure('Name','Highpass IIR Filter Chebyshov I','NumberTitle', 'off')
plot_iir(b2,a2,Fs)     % уюпюйрепхярхйх аху-тхкэрпю тбв веашьебю II пндю
figure('Name','HighpassIIR Filter Chebyshov II','NumberTitle', 'off')
plot_iir(b3,a3,Fs)     % уюпюйрепхярхйх аху-тхкэрпю тбв веашьебю II пндю
figure('Name','Highpass IIR Filter Elliptic','NumberTitle', 'off')
plot_iir(b4,a4,Fs)     % уюпюйрепхярхйх аху-тхкэрпю тбв гнкнрюпебю-йюсщпю
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.4. яхмрег юто аюррепбнпрю, веашьебю I х II  пндю х гнкнрюпебю-йюсщпю')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ЦПЮМХВМШУ ВЮЯРНР юто тбв ог (Fk) Х оо (Ft) МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
Ft = (Fs/pi)*tan(pi*ft/Fs); Fk = (Fs/pi)*tan(pi*fk/Fs); % цпюмхвмше вюярнрш оо Х ог юто
disp(['   Fk = ' num2str(Fk),'   Ft = ' num2str(Ft)])
disp('%')
disp('%')
disp('% дКЪ ЯХМРЕГЮ юто тбв МЮФЛХРЕ <ENTER>')
pause
Wp = 2.*pi.*Ft; Ws = 2.*pi.*Fk;   % цпюмхвмше йпсцнбше вюярнрш оо Х ог юто
[Ra1,Wn1] = buttord(Wp,Ws,rp,rs,'s');   % онпъднй х вюярнрю япегю юто тбв аюррепбнпрю
[Ra2,Wn2] = cheb1ord(Wp,Ws,rp,rs,'s');  % онпъднй х вюярнрю япегю юто тбв веашьебю I пндю
[Ra3,Wn3] = cheb2ord(Wp,Ws,rp,rs,'s');  % онпъднй х вюярнрю япегю юто тбв веашьебю II пндю
[Ra4,Wn4] = ellipord(Wp,Ws,rp,rs,'s');  % онпъднй х вюярнрю япегю юто тбв гнкнрюпебю-йюсщпю
[bs1,as1] = butter(Ra1,Wn1,'high','s');      % йнщттхжхемрш юто тбв аюррепбнпрю
[bs2,as2] = cheby1(Ra2,rp,Wn2,'high','s');   % йнщттхжхемрш юто тбв веашьебю I пндю
[bs3,as3] = cheby2(Ra3,rs,Wn3,'high','s');   % йнщттхжхемрш юто тбв веашьебю II пндю
[bs4,as4] = ellip(Ra4,rp,rs,Wn4,'high','s'); % йнщттхжхемрш юто тбв веашьебю гнкнрюпебю-йюсщпю
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ОНПЪДЙНБ юто тбв МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp(['   Ra1 = ' num2str(Ra1),'   Ra2 = ' num2str(Ra2),'   Ra3 = ' num2str(Ra3),'   Ra4 = ' num2str(Ra4)])
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.5. бшбнд цпютхйнб юву юто аюррепбнпрю, веашьебю I х II пндю х гнкнрюпебю-йюсщпю')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юву юто МЮФЛХРЕ <ENTER>')
pause
f = 0:((Fs/2)/1000):Fs/2;           % яерйю вюярнр дкъ цпютхйю юву
W = 2.*pi.*f;
Ha1 = freqs(bs1,as1,W);             % ву юто аюррепбнпрю
Ha2 = freqs(bs2,as2,W);             % ву юто веашьебю I пндю
Ha3 = freqs(bs3,as3,W);             % ву юто веашьебю II пндю
Ha4 = freqs(bs4,as4,W);             % ву юто гнкнрюпебю-йюсщпю
figure('Name','Highpass Analog Filter - Magnitude','NumberTitle', 'off')
subplot(2,2,1),plot(f,abs(Ha1)),xlabel('f(Hz)'),grid,...
ylabel('MAGNITUDE'),title('Analog Filter Butterworth'),ylim([0 1.2])
subplot(2,2,2),plot(f,abs(Ha2)),xlabel('f(Hz)'),grid,...
ylabel('MAGNITUDE'),title('Analog Filter Chebyshov I'),ylim([0 1.2])
subplot(2,2,3),plot(f,abs(Ha3)),xlabel('f(Hz)'),grid,...
ylabel('MAGNITUDE'),title('Analog Filter Chebyshov II'),ylim([0 1.2])
subplot(2,2,4),plot(f,abs(Ha4)),xlabel('f(Hz)'),grid,...
ylabel('MAGNITUDE'),title('Analog Filter Elliptic'),ylim([0 1.2])
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.6. нохяюмхе рпеанбюмхи й юву аху-ТХКЭРПЮ б бхде назейрю fdesign')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ябниярб назейрю fdesign МЮФЛХРЕ <ENTER>')
pause
MAG_highpass = fdesign.highpass('Fst,Fp,Ast,Ap',fk,ft,rs,rp,[Fs]) % назейр fdesign дкъ тбв
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.7. яхмрег аху-тхкэрпю тбв гнкнрюпебю-йюсщпю б бхде назейрю dfilt')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ябниярб назейрю dfilt МЮФЛХРЕ <ENTER>')
pause
F_highpass = design(MAG_highpass,'ellip','MatchExactly', 'both','FilterStructure','df2sos')     % тбв б бхде назейрю dfilt
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.8. люяьрюахпнбюмхе б йюяйюдмни ярпсйрспе аху-тхкэрпю тбв гнкнрюпебю-йюсщпю')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ябниярб назейрю dfilt ОНЯКЕ ЛЮЯЬРЮАХПНБЮМХЪ МЮФЛХРЕ <ENTER>')
pause
F_highpass_scale = scale(F_highpass,'L2')   % пегскэрюр люяьрюахпнбюмхъ
disp('%')
disp('%')
disp('% яхмрег аху-тхкэрпю тбв гюбепьем')


































