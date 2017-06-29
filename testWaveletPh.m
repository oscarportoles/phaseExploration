% load data from Cohen's book
load sampleEEGdata
EEGo = EEG;

 % Time specifications:
Fs = 250;                   % samples per second
dt = 1/Fs;                   % seconds per sample
StopTime = 1;             % seconds
t = (0:dt:StopTime-dt)';     % seconds
ph = 0;                     % phase in rad
step = [0, 10, -10];

for r = 1:3
     % Sine wave:
    Fc = 5;                     % hertz
%     x = sin(2*pi*Fc*t+ph) + step(r);
    EEG.data = 0.1*squeeze(EEGo.data(47,:,1));% + step(r) * ones(1,EEG.pnts);
%     EEG.srate = Fs;
%     EEG.data = x';
%     EEG.pnts = length(EEG.data);
%     EEG.times = t;

    % create wavelet
    frequency = 5; % in Hz, as usual
    time = -1:1/EEG.srate:1;
    s    = (4/(2*pi*frequency))^2; % note that s is squared here rather than in the next line...
    wavelet = exp(2*1i*pi*frequency.*time) .* exp(-time.^2./(2*s)/frequency);

    % FFT parameters
    n_wavelet            = length(wavelet);
    n_data               = EEG.pnts;
    n_convolution        = n_wavelet+n_data-1;
    half_of_wavelet_size = (length(wavelet)-1)/2;

    % FFT of wavelet and EEG data
    fft_wavelet = fft(wavelet,n_convolution);
    fft_data    = fft(squeeze(EEG.data),n_convolution); % FCz, trial 1

    convolution_result_fft = ifft(fft_wavelet.*fft_data,n_convolution) * sqrt(s);

    % cut off edges
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);

    % % plot for comparison
    % figure
    % subplot(211)
    % % plot(EEG.times,real(convolution_result_fft))
    % plot(EEG.times,EEG.data)
    % hold on
    % plot(EEG.times,zeros(size(EEG.data)),'k')
    % xlabel('Time (ms)'), ylabel('Voltage (\muV)')
    % title([ 'Projection onto real axis is filtered signal at ' num2str(frequency) ' Hz.' ])
    % 
    % % subplot(312)
    % % plot(EEG.times,abs(convolution_result_fft).^2)
    % % xlabel('Time (ms)'), ylabel('Power (\muV^2)')
    % % title([ 'Magnitude of projection vector squared is power at ' num2str(frequency) ' Hz.' ])
    % 
    % subplot(212)
    % plot(EEG.times,angle(convolution_result_fft))
    % hold on
    % plot(EEG.times,(pi/2).*ones(size(EEG.data)),'k')
    % plot(EEG.times,(-pi/2).*ones(size(EEG.data)),'k')
    % xlabel('Time (ms)'), ylabel('Phase angle (rad.)')
    % title([ 'Angle of vector is phase angle time series at ' num2str(frequency) ' Hz.' ])

    figure(r)
%     subplot(3,1,r)
    plot(EEG.times,real(convolution_result_fft))
%     plot(EEG.times,EEG.data,'b')
    hold on
    plot(EEG.times,angle(convolution_result_fft),'.g')
%     plot(EEG.times,angle(hilbert(EEG.data)),'.r')
    plot(EEG.times,angle(hilbert(real(convolution_result_fft))),'.r')
    legend({'amplitude','wavelet','hilbert'})
    plot(EEG.times,(pi/2).*ones(size(EEG.data)),'k')
    plot(EEG.times,(-pi/2).*ones(size(EEG.data)),'k')
    plot(EEG.times,zeros(size(EEG.data)),'k')
    xlabel('Time (ms)'), ylabel('Voltage (\muV) ; angle (rad)')
    title(['EEG bandpass fileterd with a wavelet kernel at ' num2str(frequency) ' Hz'])
%     title(['sin(x) + ' num2str(step(r))])
end