% Example of fitting an AR model

offset = 1*3e5;										% time offset
plot_ind = 50;										% plot series
binwidth = 30;										% in samples (30 samples = 1ms)

x = sum(reshape(sample1((1:3e4)+offset), 30, []));	% bin data
figure(plot_ind + 2); plot(x)						% plot binned

[xc,lags] = xcorr(x,100,'coeff');					% Get AC sequence
figure(plot_ind + 3);								% ... and plot
stem(lags(51:end),xc(51:end),'filled')
xlabel('Lag')
ylabel('ACF')
title('Sample Autocorrelation Sequence')
grid; hold off

p = 200;											% Fit AR(p) model
[arcoefs,E,K] = aryule(x,p);
pacf = -K;											% PAC sequence is negative reflection coeff's
figure(plot_ind + 4);								% ... and plot
stem(pacf,'filled')
xlabel('Lag'); ylabel('Partial ACF')
title('Partial Autocorrelation Sequence')
xlim([1 p])
uconf = 1.96/sqrt(length(x));						% ... with confidence bounds
lconf = -uconf;
hold on
plot([1 p],[1 1]'*[lconf uconf],'r')
grid; hold off

figure(plot_ind + 5)
stem(arcoefs);
title('AR Coefficients')
xlim([1 p])
grid on