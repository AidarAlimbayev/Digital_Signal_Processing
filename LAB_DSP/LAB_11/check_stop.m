function [dp1,ds,dp2] = check_stop(b,ft1,fk1,fk2,ft2,Fs)
% �������� ���������� ���������� � ��� �� 
%  
% b � ������ ������������� ���-������� �� 
% ft1,fk1,fk2,ft2 � ��������� ������� ��1, �� � ��2
% Fs � ������� ������������� (��)
%  
% dp1,ds,dp2 � ������������ ���������� ��� � �� � ��
% fp1,fs,fp2 � ������� ������ (��) ��� ��1, �� � ��2 (������ �����)
% H � ��������� ��������������
% a = [1] � ����������� ����������� ������������ ������� 
%  
a = [1];
fp1 = 0:ft1/1000:ft1;
H = freqz(b,a,fp1,Fs);
dp1 = max([max(abs(H))-1 1-min(abs(H))]);
fs = fk1:(fk2-fk1)/1000:fk2;
H = freqz(b,a,fs,Fs);
ds = max(abs(H));
fp2 = ft2:(Fs/2-ft2)/1000:Fs/2;
H = freqz(b,a,fp2,Fs);
dp2 = max([max(abs(H))-1 1-min(abs(H))]);


