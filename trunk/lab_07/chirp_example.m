 t=0:0.001:2;                    % 2 secs @ 1kHz sample rate
 y=chirp(t,0,1,150);             % Start @ DC, cross 150Hz at t=1sec
 spectrogram(y,256,250,256,1E3); % Display the spectrogram
 