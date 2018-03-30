 fs = 100;                               % Sampling frequency
 t = 0:1/fs:1;                           % Time vector
 x = sin(2*pi*t*3)+.25*sin(2*pi*t*40);   % Input Signal
 b = ones(1,10)/10;  % 10 point averaging filter
 y = fftfilt(b,x);   % FIR filtering using overlap-add method
 plot(t,x,t,y,'--');
 legend('Original Signal','Filtered Signal')