script
clc
clear
disp('% кп ╧13. яхмрег аху-тхкэрпю пт лернднл ахкхмеимнцн Z-опенапюгнбюмхъ')
disp('%')
disp('%')
disp('% О.1. ббнд рпеанбюмхи й юву (Да) пт')
disp('%')
disp('%')
disp('% бБЕДХРЕ мнлеп апхцюдш Х рпеанбюмхъ Й юву (Да)')
DATA=0;
while DATA==0;
Nb = input('Nb = ');         % мнлеп апхцюдш
Fs = input('Fs = ');         % вюярнрю дхяйперхгюжхх Б цЖ
ft1 = input('ft1 = ');       % цпюмхвмюъ вюярнрю оо1 Б цЖ
fk1 = input('fk1 = ');       % цпюмхвмюъ вюярнрю ог1 Б цЖ
fk2 = input('fk2 = ');       % цпюмхвмюъ вюярнрю ог2 Б цЖ
ft2 = input('ft2 = ');       % цпюмхвмюъ вюярнрю ог2 Б цЖ
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
disp('% дКЪ ЯХМРЕГЮ аху-ТХКЭРПНБ пт МЮФЛХРЕ <ENTER>')
pause
ft = [ft1 ft2]; fk = [fk1 fk2];        % бейрнпш цпюмхвмшу вюярнр оо Х ог
WDp = ft/(Fs/2); WDs = fk/(Fs/2);      % бейрнпш цпюмхвмшу мнплхпнбюммшу вюярнр оо Х ог
[R1,WDn1] = buttord(WDp,WDs,rp,rs);    % онпъднй х вюярнрш япегю аху-тхкэрпю пт аюррепбнпрю
[R2,WDn2] = cheb1ord(WDp,WDs,rp,rs);   % онпъднй х вюярнрш япегю аху-тхкэрпю пт веашьебю I пндю
[R3,WDn3] = cheb2ord(WDp,WDs,rp,rs);   % онпъднй х вюярнрш япегю аху-тхкэрпю пт веашьебю II пндю
[R4,WDn4] = ellipord(WDp,WDs,rp,rs);   % онпъднй х вюярнрш япегю аху-тхкэрпю пт гнкнрюпебю-йюсщпю
[b1,a1] = butter(R1,WDn1,'stop');      % йнщттхжхемрш аху-тхкэрпю пт аюррепбнпрю
[b2,a2] = cheby1(R2,rp,WDn2,'stop');   % йнщттхжхемрш аху-тхкэрпю пт веашьебю I пндю
[b3,a3] = cheby2(R3,rs,WDn3,'stop');   % йнщттхжхемрш аху-тхкэрпю пт веашьебю II пндю
[b4,a4] = ellip(R4,rp,rs,WDn4,'stop'); % йнщттхжхемрш аху-тхкэрпю пт гнкнрюпебю-йюсщпю
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ОНПЪДЙНБ аху-ТХКЭРПНБ пт МЮФЛХРЕ <ENTER>')
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
disp('% О.3. юмюкхг уюпюйрепхярхй аху-тхкэрпнб пт')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ уюпюйрепхярхй аху-тхкэрпнб пт (вершпе цпютхвеяйху нймю) МЮФЛХРЕ <ENTER>')
pause
figure('Name','Bandstop IIR Filter Butterworth','NumberTitle', 'off')
plot_iir(b1,a1,Fs)     % уюпюйрепхярхйх аху-тхкэрпю пт аюррепбнпрю
figure('Name','Bandstop IIR Filter Chebyshov I','NumberTitle', 'off')
plot_iir(b2,a2,Fs)     % уюпюйрепхярхйх аху-тхкэрпю пт веашьебю II пндю
figure('Name','Bandstop IIR Filter Chebyshov II','NumberTitle', 'off')
plot_iir(b3,a3,Fs)     % уюпюйрепхярхйх аху-тхкэрпю пт веашьебю II пндю
figure('Name','Bandstop IIR Filter Elliptic','NumberTitle', 'off')
plot_iir(b4,a4,Fs)     % уюпюйрепхярхйх аху-тхкэрпю пт гнкнрюпебю-йюсщпю
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.4. яхмрег юто аюррепбнпрю, веашьебю I х II пндю х гнкнрюпебю-йюсщпю')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ЦПЮМХВМШУ ВЮЯРНР юто пт оо1 (Ft1), ог1 (Fk1), ог2 (Fk2) Х оо2 (Ft2) МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
ft = [ft1 ft2]; fk = [fk1 fk2];       % бейрнпш цпюмхвмшу вюярнр оо Х ог аху-тхкэрпю
Ft = (Fs/pi)*tan(pi*ft/Fs); Fk = (Fs/pi)*tan(pi*fk/Fs); % бейрнпш цпюмхвмшу вюярнр оо Х ог юто
disp(['   Ft1 = ' num2str(Ft(1)),'   Fk1 = ' num2str(Fk(1)),'   Fk2 = ' num2str(Fk(2)),'   Ft2 = ' num2str(Ft(2))])
disp('%')
disp('%')
disp('% дКЪ ЯХМРЕГЮ юто пт МЮФЛХРЕ <ENTER>')
pause
Wp = 2.*pi.*Ft; Ws = 2.*pi.*Fk;            % бейрнпш цпюмхвмшу йпсцнбшу вюярнр оо Х ог юто
[Ra1,Wn1] = buttord(Wp,Ws,rp,rs,'s');      % онпъднй х вюярнрш япегю юто пт аюррепбнпрю
[Ra2,Wn2] = cheb1ord(Wp,Ws,rp,rs,'s');     % онпъднй х вюярнрш япегю юто пт веашьебю I пндю
[Ra3,Wn3] = cheb2ord(Wp,Ws,rp,rs,'s');     % онпъднй х вюярнрш япегю юто пт веашьебю II пндю
[Ra4,Wn4] = ellipord(Wp,Ws,rp,rs,'s');     % онпъднй х вюярнрш япегю юто пт гнкнрюпебю-йюсщпю
[bs1,as1] = butter(Ra1,Wn1,'stop','s');      % йнщттхжхемрш юто пт аюррепбнпрю
[bs2,as2] = cheby1(Ra2,rp,Wn2,'stop','s');   % йнщттхжхемрш юто пт веашьебю I пндю
[bs3,as3] = cheby2(Ra3,rs,Wn3,'stop','s');   % йнщттхжхемрш юто пт веашьебю II пндю
[bs4,as4] = ellip(Ra4,rp,rs,Wn4,'stop','s'); % йнщттхжхемрш юто пт гнкнрюпебю-йюсщпю
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ОНПЪДЙНБ юто пт МЮФЛХРЕ <ENTER>')
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
f =0 :((Fs/2)/1000):Fs/2;            % яерйю вюярнр дкъ цпютхйю юву
W = 2.*pi.*f;
Ha1 = freqs(bs1,as1,W);              % ву юто аюррепбнпрю
Ha2 = freqs(bs2,as2,W);              % ву юто веашьебю I пндю
Ha3 = freqs(bs3,as3,W);              % ву юто веашьебю II пндю
Ha4 = freqs(bs4,as4,W);              % ву юто гнкнрюпебю-йюсщпю
figure('Name','Bandstop Analog Filter - Magnitude','NumberTitle', 'off')
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
MAG_bandstop = fdesign.bandstop('Fp1,Fst1,Fst2,Fp2,Ap1,Ast,Ap2',ft1,fk1,fk2,ft2,rp,rs,rp,[Fs]) % назейр fdesign дкъ пт
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.7. яхмрег аху-тхкэрпю пт гнкнрюпебю-йюсщпю б бхде назейрю dfilt')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ябниярб назейрю dfilt МЮФЛХРЕ <ENTER>')
pause
F_bandstop = design(MAG_bandstop,'ellip','MatchExactly', 'both','FilterStructure','df2sos')     % пт б бхде назейрю dfilt
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause
disp('%')
disp('%')
disp('% О.8. люяьрюахпнбюмхе б йюяйюдмни ярпсйрспе аху-тхкэрпю пт гнкнрюпебю-йюсщпю')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ябниярб назейрю dfilt ОНЯКЕ ЛЮЯЬРЮАХПНБЮМХЪ МЮФЛХРЕ <ENTER>')
pause
F_bandstop_scale = scale(F_bandstop,'L2') % пегскэрюр люяьрюахпнбюмхъ
disp('%')
disp('%')
disp('% яхмрег аху-тхкэрпю пт гюбепьем')




































