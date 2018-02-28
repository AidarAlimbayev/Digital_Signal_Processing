script
clc
clear
disp('% �� �11. ������ ���-������� �� ������� ����')
disp('%')
disp('%')
disp('% �.1. ���� ���������� � ��� ��')
disp('%')
disp('%')
disp('% ������� ����� ������� � ���������� � ���')
DATA=0;
while DATA==0;
Nb = input('Nb = ');       % ����� �������
Fs = input('Fs = ');       % ������� ������������� (��)
ft1 = input('ft1 = ');     % ��������� ������� ��1 (��)
fk1 = input('fk1 = ');     % ��������� ������� ��1 (��)
fk2 = input('fk2 = ');     % ��������� ������� ��2 (��)
ft2 = input('ft2 = ');     % ��������� ������� ��2 (��)
d11 = input('d11 = ');     % ����������� ���������� ���������� � ��1
d2 = input('d2 = ');       % ����������� ���������� ���������� � ��
d12 = input('d12 = ');     % ����������� ���������� ���������� � ��2 
disp('% ��������� ������������ ����� �������� ������')
disp('% ��� ���������� �������� ������ ������� 1')
disp('% ��� ������������ �������� ������ ������� 0 � ��������� ����')
DATA = input('--> ');
end
disp('%')
disp('%')
disp('% ��� ����������� ������� <ENTER>')
pause 
disp('%')
disp('%')
disp('% �.2. ���������� ���������� ������� kaiserord')
disp('%')
disp('%')
disp('% ��� ������ ���������� ������� kaiserord ������� <ENTER>')
pause 
m = [1 0 1];             % ������ �������� ��������� ��� 
f = [ft1 fk1 fk2 ft2];   % ������ ��������� ������
ripple = [d11 d2 d12];   % ������ ����������� ���������� ����������  
[R,wc,beta,ftype] = kaiserord(f,m,ripple,Fs);    % ���������� ���������� ���� �������
disp(['R = ' num2str(R)])               % ������ ������� ���-�������
disp(['wc(1) = ' num2str(wc(1)) '      wc(2) = ' num2str(wc(2))]) % ������ ������������� ������ �������
disp(['beta = ' num2str(beta)])         % �������� ���� �������
disp(['ftype = ' char(ftype)])          % ��� ���-�������
disp('%')
disp('%')
disp('% ��� ����������� ������� <ENTER>')
pause
disp('%')
disp('%')
disp('% �.3. ������ ���-������� ��')
ORDER = 0; % ������� ������������� ������� ���-�������: 0 � �������������; 1 - �����������
while ORDER==0;
disp('%')
disp('%')
disp('% ��� ������� ���-������� �� ������� <ENTER>')
pause
b4 = fir1(R,wc,ftype,kaiser(R+1,beta),'noscale');     % ������������ ���-������� ��
disp('%')
disp('%')
disp(['  ������������ ���-������ �� ������� R = ' num2str(R)])
disp('%')
disp('%')
disp('% ��� ������ ����������� ������������ ���������� ���')
disp('% � ��1 (dp1), �� (ds) � ��2 (dp2) � �������� ���������� d11, d2 � d12 ������� <ENTER>')
pause
[dp1,ds,dp2] = check_stop(b4,ft1,fk1,fk2,ft2,Fs);     % ���������� ����������� ������������ �� ������ ���������� � ��1, �� � ��2
disp('%')
disp(['dp1=' num2str(dp1) '      ds = ' num2str(ds) '      dp2 = ' num2str(dp2)])
disp(['d11 = ' num2str(d11) '      d2 = ' num2str(d2) '      d12 = ' num2str(d12)])
disp('%')
disp('%')
disp('% �������� ����������� ���������� � ���������')
disp('%')
disp('% ���� ������� ������������� ������������, ������� 1')
disp('% ���� �� �������������, ������� 0 � ����� ������� R')
ORDER = input('--> ');
if ORDER==0 
R = input('R = ');                            % ������� ���-������� 
while rem(R,2)~=0 
disp('% ������� ������� ������ �����������')
R = input('R = ');                            % ������� ���-������� 
end
end
end
disp('%')
disp(['  ������������ �� ������������ ������� R = ' num2str(R)])
disp('%')
disp('%')
disp('% ��� ����������� ������� <ENTER>')
pause
disp('%')
disp('%')
disp('% �.4. ������ ������������� ���-������� ��')
disp('%')
disp('%')
disp('% ��� ������ �������� ��, ��� � ��� ������� <ENTER>')
pause
figure('Name','Bandstop FIR Filter - Impulse Response, Magnitude, Phase','NumberTitle', 'off')
plot_fir(R,b4,Fs)             % ���������� �������� ��, ��� � ���  
disp('%')
disp('%')
disp('% ��� ����������� ������� <ENTER>')
pause 
disp('%')
disp('%')
disp('% �.5. �������� ��������� ���-������� �� � ���� ������� dfilt')
disp('%')
disp('%')
disp('% ��� ������ ������� ������� dfilt ������� <ENTER>')
pause
F_bandstop = dfilt.dfsymfir(b4)     % ������ dfilt � ���-������ ��
disp('%')
disp('%')
disp('% ������ ���-������� �� ��������')

































