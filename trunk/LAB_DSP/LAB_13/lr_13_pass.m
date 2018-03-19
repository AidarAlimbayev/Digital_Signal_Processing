script
clc
clear
disp('% кп ╧13. яхмрег аху-тхкэрпю от лернднл ахкхмеимнцн Z-опенапюгнбюмхъ')
disp('%')
disp('%')
disp('% О.1. ббнд рпеанбюмхи й юву (Да) от')
disp('%')
disp('%')
disp('% бБЕДХРЕ мнлеп апхцюдш Х рпеанбюмхъ Й юву (Да)')
DATA=0;
while DATA==0;
Nb = input('Nb = ');         % мнлеп апхцюдш
Fs = input('Fs = ');         % вюярнрю дхяйперхгюжхх (цЖ)
fk1 = input('fk1 = ');       % цпюмхвмюъ вюярнрю ог1 (цЖ)
ft1 = input('ft1 = ');       % цпюмхвмюъ вюярнрю оо1 (цЖ)
ft2 = input('ft2 = ');       % цпюмхвмюъ вюярнрю ог2 (цЖ)
fk2 = input('fk2 = ');       % цпюмхвмюъ вюярнрю ог2 (цЖ)
rp = input('rp = ');         % люйяхлюкэмн дносярхлне гюрсуюмхе б оо
rs = input('rs = ');         % лхмхлюкэмн дносярхлне гюрсуюмхе б ог
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
disp('% дКЪ ЯХМРЕГЮ аху-ТХКЭРПНБ от МЮФЛХРЕ <ENTER>')
pause
ft = [ft1 ft2]; fk = [fk1 fk2];       % бейрнпш цпюмхвмшу вюярнр оо Х ог
WDp = ft/(Fs/2); WDs = fk/(Fs/2);     % бейрнпш цпюмхвмшу мнплхпнбюммшу вюярнр оо Х ог
[R1,WDn1] = buttord(WDp,WDs,rp,rs);   % онпъднй х вюярнрш япегю аху-тхкэрпю от аюррепбнпрю
[R2,WDn2] = cheb1ord(WDp,WDs,rp,rs);  % онпъднй х вюярнрш япегю аху-тхкэрпю от веашьебю I пндю
[R3,WDn3] = cheb2ord(WDp,WDs,rp,rs);  % онпъднй х вюярнрш япегю аху-тхкэрпю от веашьебю II пндю
[R4,WDn4] = ellipord(WDp,WDs,rp,rs);  % онпъднй х вюярнрш япегю аху-тхкэрпю от гнкнрюпебю-йюсщпю
[b1,a1] = butter(R1,WDn1);            % йнщттхжхемрш аху-тхкэрпю от аюррепбнпрю
[b2,a2] = cheby1(R2,rp,WDn2);         % йнщттхжхемрш аху-тхкэрпю от веашьебю I пндю
[b3,a3] = cheby2(R3,rs,WDn3);         % йнщттхжхемрш аху-тхкэрпю от веашьебю II пндю
[b4,a4] = ellip(R4,rp,rs,WDn4);       % йнщттхжхемрш аху-тхкэрпю от гнкнрюпебю-йюсщпю
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ОНПЪДЙНБ аху-ТХКЭРПНБ от МЮФЛХРЕ <ENTER>')
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
disp('% О.3. юмюкхг уюпюйрепхярхй аху-тхкэрпнб от')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ уюпюйрепхярхй аху-тхкэрпнб от (вершпе цпютхвеяйху нймю) МЮФЛХРЕ <ENTER>')
pause
figure('Name','Bandpass IIR Filter Butterworth','NumberTitle', 'off')
plot_iir(b1,a1,Fs)       % уюпюйрепхярхйх аху-тхкэрпю от аюррепбнпрю
figure('Name','Bandpass IIR Filter Chebyshov I','NumberTitle', 'off')
plot_iir(b2,a2,Fs)       % уюпюйрепхярхйх аху-тхкэрпю от веашьебю II пндю
figure('Name','Bandpass IIR Filter Chebyshov II','NumberTitle', 'off')
plot_iir(b3,a3,Fs)       % уюпюйрепхярхйх аху-тхкэрпю от веашьебю II пндю
figure('Name','Bandpass IIR Filter Elliptic','NumberTitle', 'off')
plot_iir(b4,a4,Fs)       % уюпюйрепхярхйх аху-тхкэрпю от гнкнрюпебю-йюсщпю
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.4. яхмрег юто аюррепбнпрю, веашьебю I х II пндю х гнкнрюпебю-йюсщпю')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ЦПЮМХВМШУ ВЮЯРНР юто от ог1 (Fk1), оо1 (Ft1), оо2 (Ft2) Х ог2 (Fk2) МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
ft = [ft1 ft2]; fk = [fk1 fk2];       % бейрнпш цпюмхвмшу вюярнр оо Х ог аху-тхкэрпю
Ft = (Fs/pi)*tan(pi*ft/Fs); Fk = (Fs/pi)*tan(pi*fk/Fs); % бейрнпш цпюмхвмшу вюярнр оо Х ог юто
disp(['   Fk1 = ' num2str(Fk(1)),'   Ft1 = ' num2str(Ft(1)),'   Ft2 = ' num2str(Ft(2)),'   Fk2 = ' num2str(Fk(2))])
disp('%')
disp('%')
disp('% дКЪ ЯХМРЕГЮ юто от МЮФЛХРЕ <ENTER>')
pause
Wp = 2.*pi.*Ft; Ws = 2.*pi.*Fk;       % бейрнпш цпюмхвмшу йпсцнбшу вюярнр оо Х ог юто
[Ra1,Wn1] = buttord(Wp,Ws,rp,rs,'s'); % онпъднй х вюярнрш япегю юто от аюррепбнпрю
[Ra2,Wn2] = cheb1ord(Wp,Ws,rp,rs,'s');% онпъднй х вюярнрш япегю юто от веашьебю I пндю
[Ra3,Wn3] = cheb2ord(Wp,Ws,rp,rs,'s');% онпъднй х вюярнрш япегю юто от веашьебю II пндю
[Ra4,Wn4] = ellipord(Wp,Ws,rp,rs,'s');% онпъднй х вюярнрш япегю юто от гнкнрюпебю-йюсщпю
[bs1,as1] = butter(Ra1,Wn1,'s');      % йнщттхжхемрш юто от аюррепбнпрю
[bs2,as2] = cheby1(Ra2,rp,Wn2,'s');   % йнщттхжхемрш юто от веашьебю I пндю
[bs3,as3] = cheby2(Ra3,rs,Wn3,'s');   % йнщттхжхемрш юто от веашьебю II пндю
[bs4,as4] = ellip(Ra4,rp,rs,Wn4,'s'); % йнщттхжхемрш юто от гнкнрюпебю-йюсщпю
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ОНПЪДЙНБ юто от МЮФЛХРЕ <ENTER>')
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
figure('Name','Bandpass Analog Filter - Magnitude','NumberTitle', 'off')
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
MAG_bandpass = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',fk1,ft1,ft2,fk2,rs,rp,rs,[Fs]) % назейр fdesign дкъ от
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.7. яхмрег аху-тхкэрпю от гнкнрюпебю-йюсщпю б бхде назейрю dfilt')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ябниярб назейрю dfilt МЮФЛХРЕ <ENTER>')
pause
F_bandpass = design(MAG_bandpass,'ellip','MatchExactly', 'both','FilterStructure','df2sos')     % от б бхде назейрю dfilt
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.8. люяьрюахпнбюмхе б йюяйюдмни ярпсйрспе аху-тхкэрпю от гнкнрюпебю-йюсщпю')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ябниярб назейрю dfilt ОНЯКЕ ЛЮЯЬРЮАХПНБЮМХЪ МЮФЛХРЕ <ENTER>')
pause
F_bandpass_scale = scale(F_bandpass)   % пегскэрюр люяьрюахпнбюмхъ
disp('%')
disp('%')
disp('% яхмрег аху-тхкэрпю от гюбепьем')
































