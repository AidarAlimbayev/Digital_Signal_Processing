% % Example 1:
% %   Design a lowpass FIR filter with normalized cut-off frequency at
% %   0.3 and determine its impulse response.
%  
% b=fircls1(54,0.3,0.02,0.008);
% impz(b)b

% Example 2:
%   Design a 5th order lowpass elliptic IIR filter and determine its
%   impulse response.
 
[b,a] = ellip(5,0.5,20,0.4);
impz(b,a)

% % Example 3:
%     %   Design a Butterworth highpass IIR filter, represent its coefficients
%     %   using second order sections, and display its impulse response.
%  
%     [z,p,k] = butter(6,0.7,'high');
%     SOS = zp2sos(z,p,k);
%     impz(SOS)

%   % Example 4:
%     %   Use the designfilt function to design a highpass IIR digital filter 
%     %   with order 8, passband frequency of 75 KHz, and a passband ripple 
%     %   of 0.2 dB. Sample rate is 200 KHz. Visualize the impulse response 
%     %   using 256 samples.
%    
%     D = designfilt('highpassiir', 'FilterOrder', 8, ...
%              'PassbandFrequency', 75e3, 'PassbandRipple', 0.2,...
%              'SampleRate', 200e3);
%  
%     impz(D,256)
 