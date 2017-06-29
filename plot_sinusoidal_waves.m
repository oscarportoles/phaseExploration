   %% Time specifications:
   Fs = 1000;                   % samples per second
   dt = 1/Fs;                   % seconds per sample
   StopTime = 0.5;             % seconds
   t = (0:dt:StopTime-dt)';     % seconds
   ph = 0;                     % phase in rad
   %% Sine wave:
   Fc = 100;                     % hertz
   x = sin(2*pi*Fc*t+ph);
   % Plot the signal versus time:
   figure(1);
   plot(t,x,'b');
   xlabel('time (in seconds)');
   title('Signal versus Time');
   zoom xon;