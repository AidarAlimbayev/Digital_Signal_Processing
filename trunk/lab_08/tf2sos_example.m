[b,a] = butter(3,.4);
      [z,p,k] = tf2zpk(b,a)

% Example:
    %   Create a 5th order butterworth filter and convert its transfer 
    %   function to second order-sections form.
 
    fc = 1000;                      % Cut-off frequency (Hz)
    fs = 9000;                      % Sampling rate (Hz)
    [b,a] = butter(5,2*fc/fs);      % Normalized butterworth filter
    [sos,g] = tf2sos(b,a);          % Second Order Section conversion
    fvtool(sos,g)                   % Visualize the filter
    
    