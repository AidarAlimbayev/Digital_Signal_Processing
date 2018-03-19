script
clc
clear
disp('% кп ╧8. кхмеимше дхяйпермше яхярелш')
disp('%')
disp('%')
disp('% бБЕДХРЕ хяундмше дюммше');
DATA=0;
while DATA==0
Nb = input('Nb = ');   % мнлеп апхцюдш
b = input('b = ');     % бейрнп йнщттхжхемрнб вхякхрекъ оепедюрнвмни тсмйжхх
a = input('a = ');     % бейрнп йнщттхжхемрнб гмюлемюрекъ оепедюрнвмни тсмйжхх
N1 = input('N1 = ');   % дкхмю хлоскэямни уюпюйрепхярхйх
N2 = input('N2 = ');   % дкхмю бнгдеиярбхъ
Fs = input('Fs = ');   % вюярнрю дхяйперхгюжхх
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
disp('% О.1. бшвхякемхе хлоскэямни уюпюйрепхярхйх - ТСМЙЖХЪ impz')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйю хлоскэямни уюпюйрепхярхйх МЮФЛХРЕ <ENTER>')
pause 
h1 = impz(b,a,N1);          % хлоскэямюъ уюпюйрепхярхйю
n = 0:(N1-1);               % дхяйпермне мнплхпнбюммне бпелъ дкъ ху
figure('Name','Impulse Response','NumberTitle', 'off')
subplot(2,1,1), stem(n,h1,'fill','MarkerSize',3), grid
xlabel('n'), ylabel('h(n)')
title('Impulse Response h(n) - impz')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause 
disp('%')
disp('%')
disp('% О.2. бшвхякемхе хлоскэямни уюпюйрепхярхйх - ТСМЙЖХЪ filter')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйю хлоскэямни уюпюйрепхярхйх МЮФЛХРЕ <ENTER>')
pause 
u0 = [1 zeros(1,(N1-1))];    % жхтпнбни едхмхвмши хлоскэя
h2 = filter(b,a,u0);         % хлоскэямюъ уюпюйрепхярхйю
subplot(2,1,2), stem(n,h2,'fill','MarkerSize',3), grid
xlabel('n'), ylabel('h(n)'), title('Impulse Response h(n) - filter')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause 
disp('%')
disp('%')
disp('% О.3. бшвхякемхе пеюйжхх он тнплске ябепрйх')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб бнгдеиярбхъ х пеюйжхх, БШВХЯКЕММНИ ОН тнплске ябепрйх,  МЮФЛХРЕ <ENTER>')
pause
x = input_1(N2);    % бнгдеиярбхе (дхяйпермши опълнсцнкэмши хлоскэя)  
y1 = conv(x,h1);    % пеюйжхъ дкхмш, пюбмни дкхме ябепрйх
L = N1+N2-1;        % дкхмю ябепрйх
n = 0:(N2-1);       % дхяйпермне мнплхпнбюммне бпелъ дкъ бнгдеиярбхъ
n1 = 0:(L-1);       % дхяйпермне мнплхпнбюммне бпелъ дкъ ябепрйх
figure('Name','Input and Output Signals','NumberTitle', 'off')
subplot(4,1,1), stem(n,x,'fill','MarkerSize',3), grid, xlabel('n')
ylabel('x(n)'), title('Input Signal - Discrete Rectangular Impulse x(n)')
subplot(4,1,2), stem(n1,y1,'fill','MarkerSize',3), grid
ylabel('y(n)'), title('Output Signal y1(n) √ conv (length = L)')
subplot(4,1,3), stem(n,y1(1:N2),'fill','MarkerSize',3), grid
xlabel('n'), ylabel('y1(n)')
title('Output Signal y1(n) √ conv (length = N2)')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause 
disp('%')
disp('%')
disp('% О.4. бшвхякемхе пеюйжхх он пюгмнярмнлс спюбмемхч')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйю пеюйжхх, БШВХЯКЕММНИ ОН пюгмнярмнлс спюбмемхч, МЮФЛХРЕ <ENTER>')
pause 
y2 = filter(b,a,x);  % пеюйжхъ кдя 
subplot(4,1,4), stem(n,y2,'fill','MarkerSize',3), grid
xlabel('n'), ylabel('y(n)')
title('Output Signal y2(n) √ filter (length = N2)')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause 
disp('%')
disp('%')
disp('% О.5. бшвхякемхе оюпюлерпнб оепедюрнвмни тсмйжхх б бхде опнхгбедемхъ опняреиьху лмнфхрекеи')
disp('%')
disp('%')
disp('%дКЪ бшбндю МСКЕИ (q) Х ОНКЧЯНБ (p) б юкцеапюхвеяйни тнпле Х ЙНЩТТХЖХЕМРЮ СЯХКЕМХЪ (K) МЮФЛХРЕ <ENTER>')
pause
[q,p,K] = tf2zpk(b,a)  % мскх (q) х онкчяш (p) б юкцеапюхвеяйни тнпле х йнщттхжхемр сяхкемхъ (K)
disp('%')
disp('%дКЪ бшбндю МСКЕИ (q) Б онйюгюрекэмни тнпле МЮФЛХРЕ <ENTER>')
pause
disp('% rq - пюдхсяш, wq - юпцслемрш МСКЕИ')
rq = abs(q)           % пюдхсяш йнлокейямн янопъфеммшу мскеи 
wq = angle(q)         % юпцслемрш йнлокейямн янопъфеммшу мскеи
disp('%дКЪ бшбндю ОНКЧЯНБ (p) Б онйюгюрекэмни тнпле МЮФЛХРЕ <ENTER>')
pause
disp('% rp - пюдхсяш, wp - юпцслемрш ОНКЧЯНБ')
rp = abs(p)           % пюдхсяш йнлокейямн янопъфеммшу онкчянб
wp = angle(p)         % юпцслемрш йнлокейямн янопъфеммшу онкчянб
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause 
disp('%')
disp('%')
disp('% О.6. бшвхякемхе оюпюлерпнб оепедюрнвмни тсмйжхх б бхде опнхгбедемхъ лмнфхрекеи брнпнцн онпъдйю') 
disp('%')
disp('%')
disp('% дКЪ бшбндю ЛЮРПХЖШ ЙНЩТТХЖХЕМРНБ (s) Х ЙНЩТТХЖХЕМРЮ СЯХКЕМХЪ (G) МЮФЛХРЕ <ENTER>')
pause
[s,G] = tf2sos(b,a) % йнщттхжхемрш (s) х йнщттхжхемр сяхкемхъ (G)
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause 
disp('%')
disp('%')
disp('% О.7. бшвхякемхе оюпюлерпнб оепедюрнвмни тсмйжхх б бхде ясллш опняршу дпнаеи')
disp('%')
disp('%')
disp('% дКЪ бшбндю ЙНЩТТХЖХЕМРНБ ПЮГКНФЕМХЪ (r), ОНКЧЯНБ (p) Х ЖЕКНИ ВЮЯРХ (c) МЮФЛХРЕ <ENTER>')
pause
[r,p,c] = residuez(b,a)  % йнщттхжхмрш пюгкнфемхъ (r) Х онкчяш (p) б юкцеапюхвеяйни тнпле х жекюъ вюярэ (c)
disp('%')
disp('%дКЪ бшбндю йнщттхжхемрнб пюгкнфемхъ (r) Б онйюгюрекэмни тнпле МЮФЛХРЕ <ENTER>')
pause
rr = abs(r)           % пюдхсяш йнлокейямн янопъфеммшу йнщттхжхемрнб пюгкнфемхъ (r)
wr = angle(r)         % юпцслемрш йнлокейямн янопъфеммшу йнщттхжхемрнб пюгкнфемхъ (r)
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause 
disp('%')
disp('%')
disp('% О.8. бшбнд йюпрш мскеи х онкчянб')
disp('%')
disp('%')
disp('% дКЪ бшбндю йюпрш мскеи х онкчянб МЮФЛХРЕ <ENTER>')
pause
figure('Name',' Z-plane zero-pole plot','NumberTitle', 'off')
zplane(b,a), title('Z-plane zero-pole plot'), grid
xlabel('Re'), ylabel('jIm')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause 
disp('%')
disp('%')
disp('% О.9.бшвхякемхе юву Х тву б ьйюке мнплхпнбюммшу вюярнр')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юву Х тву Я ЬЙЮКЕ мнплхпнбюммшу ВЮЯРНР МЮФЛХРЕ <ENTER>')
pause
w = 0:pi/100:pi;       % бейрнп мнплхпнбюммшу вюярнр (пюд)
H_w = freqz(b,a,w);    % йнлокейямюъ вюярнрмюъ уюпюйрепхярхйю
MAG_w = abs(H_w);      % юву
PHASE_w = angle(H_w);  % тву
figure('Name','Magnitude and Phase Responses','NumberTitle', 'off')
subplot(2,2,1), plot(w,MAG_w), grid, xlabel('w (rad)'), title('MAGNITUDE - |м(w)|')
subplot(2,2,3), plot(w,PHASE_w), grid, xlabel('w (rad)'), title('PHASE √ arg [H(w)]  (rad)')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause 
disp('%')
disp('%')
disp('% О.10. бшвхякемхе юву Х тву б ьйюке юаянкчрмшу вюярнр')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ цпютхйнб юву Х тву Б ЬЙЮКЕ юаянкчрмшу ВЮЯРНР МЮФЛХРЕ <ENTER>')
pause 
f = 0:Fs/100:Fs/2;     % бейрнп юаянкчрмшу вюярнр (цЖ)
H = freqz(b,a,f,Fs);   % йнлокейямюъ вюярнрмюъ уюпюйрепхярхйю
MAG = abs(H);          % юву
PHASE = angle(H);      % тву
subplot(2,2,2), plot(f,MAG), grid, xlabel('f (Hz)'), title('MAGNITUDE - |м(f)|')
subplot(2,2,4), plot(f,PHASE), grid, xlabel('f (Hz)'), title('PHASE √ arg [H(f)] (rad)')
disp('%')
disp('%')
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause 
disp('%')
disp('%')
disp('% О.11. нохяюмхе ярпсйрспш пейспяхбмнцн гбемю')
disp('%')
disp('%')
disp('% дКЪ БШБНДЮ ябниярб назейрнб dfilt МЮФЛХРЕ <ENTER>')
pause
Hd1 = dfilt.df1(b,a)   % опълюъ ярпсйрспю (Direct-Form I)
Hd2 = dfilt.df2(b,a)   % опълюъ йюмнмхвеяйюъ ярпсйрспю (Direct-Form II) 
Hd3 = dfilt.df1t(b,a)  % опълюъ рпюмяонмхпнбюммюъ ярпсйрспю (Direct-Form I Transposed)
Hd4 = dfilt.df2t(b,a)  % опълюъ йюмнмхвеяйюъ рпюмяонмхпнбюммюъ ярпсйрспю (Direct-Form I Transposed) 
disp('% дКЪ ОПНДНКФЕМХЪ МЮФЛХРЕ <ENTER>')
pause 
disp('%')
disp('%')
disp('% О.12. юмюкхг бкхъмхъ мскеи х онкчянб мю бхд юву')
disp('%')
disp('%')
b(1,:) = [1 0 0];      % йнщттхжхемрш вхякхрекъ √ 1-Ъ ярпнйю люрпхжш
b(2,:) = [1 0 0];      % йнщттхжхемрш вхякхрекъ √ 2-Ъ ярпнйю люрпхжш
b(3,:) = [1 0 0];      % йнщттхжхемрш вхякхрекъ √ 3-Ъ ярпнйю люрпхжш
b(4,:) = [1 1 0];      % йнщттхжхемрш вхякхрекъ √ 4-Ъ ярпнйю люрпхжш
a(1,:) = a;                   % йнщттхжхемрш гмюлемюрекъ √ 1-Ъ ярпнйю люрпхжш 
a(2,:)=[1 -a(1,2) a(1,3)];    % йнщттхжхемрш гмюлемюрекъ √ 2-Ъ ярпнйю люрпхжш
a(3,:)=[1 a(1,2) 1.2*a(1,3)]; % йнщттхжхемрш гмюлемюрекъ √ 3-Ъ ярпнйю люрпхжш
a(4,:)=[1 a(1,2) a(1,3)];     % йнщттхжхемрш гмюлемюрекъ √ 4-Ъ ярпнйю люрпхжш
w = 0:pi/100:pi;              % бейрнп мнплхпнбюммшу вюярнр (пюд)
for i=1:4 
H3(:,i) = freqz(b(i,:),a(i,:),w);   % вюярнрмюъ уюпюйрепхярхйю √ i-И ярнкаеж люрпхжш  
MAG3(:,i) = abs(H3(:,i)); MAX(:,i) = max(MAG3(:,i));   % юву  √ i-И ярнкаеж люрпхжш √ х люйяхлсл юву 
MAGN(:,i) = MAG3(:,i)/MAX(:,i);     % мнплхпнбюммюъ юву √ i-И ярнкаеж люрпхжш
end
disp('%  дКЪ БШБНДЮ йюпрш мскеи х онкчянб Х мнплхпнбюммни юву МЮФЛХРЕ <ENTER>')
pause
figure('Name','Z-plane zero-pole plots and Normalized Magnitudes','NumberTitle', 'off')
for i = 1:4
subplot(4,2,2*i-1), zplane(b(i,:),a(i,:)), title('Z-plane zero-pole plot'), grid
xlabel('Re'), ylabel('jIm')
subplot(4,2,2*i), plot(w,MAGN(:,i)), grid
xlabel('w (rad)'), title('Normalized Magnitude A(w)')
end
disp('%')
disp('%')
disp('% пюанрю гюбепьемю')

































