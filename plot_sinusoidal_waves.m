   %% Time specifications:
   Fs = 12500;                   % samples per second
   dt = 1/Fs;                   % seconds per sample
   StopTime = 0.05;             % seconds
   t = (0:dt:StopTime-dt)';     % seconds
   ph = 0;                     % phase in rad
   %% Sine wave:
   Fc = 10;                     % hertz
   x = sin(2*pi*Fc*t+ph);
   % Plot the signal versus time:
   figure;
   plot(t,x);
   xlabel('time (in seconds)');
   title('Signal versus Time');
   zoom xon;