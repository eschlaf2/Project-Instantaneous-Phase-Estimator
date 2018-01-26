% Phase tracker 

%% Load and filter the sample LFP data

PLOT = true; % debugging plots

load('/Volumes/NO NAME/DATA/LFPsamples/sample1.mat'); % LFP data from Alik/Ethan
k = 10;						 % downsampling factor
sample1 = sample1(1:k:end);  % subsample
sample1 = medfilt1(sample1); % median filter to reduce spurious spikes

fs = 3e4 / k;					% sample rate
low = 6;						% lower frequency bound
high = 10;						% upper frequency bound

low2 = 30;
high2 = 40;						% second frequency band

% Filter the signal
d = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',low,'HalfPowerFrequency2',high, ...
    'SampleRate',fs);
d2 = designfilt('bandpassiir', 'FilterOrder', 4, ...
	'HalfPowerFrequency1', low2, 'HalfPowerFrequency2', high2, ...
	'SampleRate', fs);

ts_filt = filtfilt(d, sample1); % Get the filtered time series
ts_filt2 = filtfilt(d2, sample1); % Get the high frequency filtered series

% Get the analytic signal, phase and amplitude
ts_analytic = hilbert(ts_filt);
ts_analytic2 = hilbert(ts_filt2);

%% Load Uri's data
RUN = false;
if RUN
load('/Volumes/NO NAME/GITHUB/Project-Instantaneous-Phase-Estimator/sample_data.mat') % hippocampal data from uri
ts_filt = EEGfilt;                        % ... and get the filtered data.
ts_analytic = EEGanalytic;
fs = 3e3;
end

%% Clip

% Select a chunk to analyze
t = (0:length(ts_filt)-1)/fs;
inds = 10001:20000;                     % choose an interval of time,
ts_filt = ts_filt(inds);				% clip the filtered signal
ts_analytic = ts_analytic(inds);		% clip the analytic signal

ts_filt2 = ts_filt2(inds);
ts_analytic2 = ts_analytic2(inds);

t = t(inds);							% time

%% Get truth values
phi = angle(ts_analytic);	% calculate the phases
M = abs(ts_analytic);		% ... and amplitudes

phi2 = angle(ts_analytic2);
M2 = abs(ts_analytic2);

if PLOT % Look at filtered and analytic signals
	figure(1); fullwidth(); ax1 = subplot(211); plot(t, ts_filt, t, M); 
	ax2 = subplot(212); plot(t, phi); linkaxes([ax1, ax2], 'x')
end

dphi_FUN = @(phi, w) diff(unwrap(phi) + w * k^2 * (rand(size(phi)) - .5)); % Compute the empirical changes in phase from the data.
dM_FUN = @(M, w) diff(M + w * k * (rand(size(M)) - .5));	% Compute the empirical changes in amplitude envelope from the data.
dphi_FUN = @(phi, w) diff(unwrap(phi) + w/2 * k^2 * (randn(size(phi)))); % Compute the empirical changes in phase from the data.
dM_FUN = @(M, w) diff(M + w/2 * k * (randn(size(M))));	% Compute the empirical changes in amplitude envelope from the data.
% dphi_FUN = @(phi, dphi, w) diff(unwrap(phi) + 1 * dphi(unidrnd(length(phi)-1, length(phi), 1))); % Compute the empirical changes in phase from the data.

dphi = dphi_FUN(phi, 5e-4);
dM = dM_FUN(M, .4);
dphi2 = dphi_FUN(phi2, 5e-4);
dM2 = dM_FUN(M2, .2);

if PLOT
	figure(11); fullwidth(); subplot(221); histogram(dM); hold on; histogram(diff(M)); hold off; title('M')
	subplot(222); histogram(dphi); hold on; histogram(diff(unwrap(phi))); hold off; title('phi');
	subplot(223); histogram(dM2); hold on; histogram(diff(M2)); hold off; title('M2')
	subplot(224); histogram(dphi2); hold on; histogram(diff(unwrap(phi2))); hold off; title('phi2');
end

%% Preprocess (IP)
RUN = false; % Show/hide section	
if RUN 
	inds = 102:1101;
	data = EEGfilt;
	phi = mod(angle(EEGanalytic(inds)),2*pi); % Compute the phase.;

	fs = 3e2;
	t = (0:length(data) - 1) / fs;

	LOWER_BOUND = [-Inf, -Inf, 0, 1, -pi/2, -Inf];
	UPPER_BOUND = [Inf, Inf, Inf, Inf, pi/2, Inf];
	A0 = [rand(1, 2) - .5, 100, 8, 0, 0];
	options = optimoptions('lsqcurvefit', 'Display', 'off');

	ahat = zeros(length(A0), length(inds));
	resnorm = zeros(1, length(inds));
	theta = zeros(1, length(inds));

	for i = 1:length(inds)
		curve_inds = ((inds(i) - 44) : (inds(i)));
		y1 = data(curve_inds - 1); 
		y2 = data(curve_inds - 2); 
		ts_filt = t(curve_inds)';
		y = data(curve_inds);
		fun = @(a, x) a(1) * y1 + a(2) * y2 + a(3) * cos(2 * pi * x * a(4) + a(5)) + a(6);
		[ahat(:, i), resnorm(i)] = lsqcurvefit(fun, A0, ts_filt, y, LOWER_BOUND, UPPER_BOUND, options);
		theta(i) = mod(2 * pi * ts_filt(end) * ahat(4, i) + ahat(5, i), 2*pi);

	end

	figure(1);
	ahat_5 = ahat(5,:);
	ahat(5,:) = ahat_5;
	subplot(131); plot(ahat', 'linewidth', 2); legend({'A1', 'A2', 'M', 'f', 'phi', 'Int'});
	subplot(132); plot(ts_filt, y, 'ko', ts_filt, fun(ahat(:, end), ts_filt), 'b-')
	subplot(133); plot(resnorm);
	set(gcf, 'units', 'normalized', 'position', [0 .5 1 .5]);
end

%% Run sequential monte carlo (SMC) a.k.a., particle filter.

nsteps = length(inds);                                  % # time steps of data.
nparts = 1e4;											% # of particles
noisestd = 1;
samples = [M(1); phi(1); M2(1); phi2(1)]*ones(1,nparts);% Create initial set of particles, they're all the same, at the first observed amp & phase.
sig = ts_filt+ts_filt2+normrnd(0,noisestd,size(ts_filt));        % ... add noise to observed signal, to vary difficulty of tracking.

Mest = M;
phiest = phi;
M2est = M2;
phi2est = phi2;

if ~exist('pf', 'var'), pf = figure(); end
for i = 2:nsteps										% For each time point of data,
    
	samples = SMC_update(...							% Update particles
		samples, sig(i), dM, dphi, noisestd, dM2, dphi2);	% Note: was 10*dM, 1.4*dphi
	Mest(i) = mean(samples(1,:));                       %   Estimate amplitude from all particles.
    phiest(i) = angle(sum(exp(1i*samples(2,:))));       %   Estimate phase from all particles.
	M2est(i) = mean(samples(3, :));
	phi2est(i) = angle(sum(exp(1i*samples(4,:))));
	
    if 0                                                %   Plot everything.
    plot(samples(1,:),samples(2,:),'.',M(i),phi(i),'r*'); axis([0 max(M) 0 2*pi]); title(num2str(i));
    xlabel('Amplitude'); ylabel('Phase')
    hold on;
    figure(pf); plot(Mest(i), phiest(i), '*g')
    hold off
    drawnow; pause(1e-2);
	end
end

figure(2); fullwidth(); plot(t,sig,t,Mest.*cos(phiest) + M2est.*cos(phi2est));
figure(3); fullwidth(); ax1 = subplot(211); plot(t, M, t, Mest); ax2 = subplot(212); plot(t, mod(phi, 2*pi), t, mod(phiest, 2*pi));
figure(4); fullwidth(); ax1 = subplot(211); plot(t, M2, t, M2est); ax2 = subplot(212); plot(t, mod(phi2, 2*pi), t, mod(phi2est, 2*pi));

%% Variable assignment

RUN = false; 
if RUN
data = s1_300_filt;
data_analytic = s1_300_analytic;

inds = 10001:20000;                       % choose an interval of time,
ts_filt = data(inds);                        % ... and get the filtered data.
M = abs(data_analytic(inds));               % Compute the amplitude envelope.
phi = mod(angle(data_analytic(inds)),2*pi); % Compute the phase.;

dphi = diff(unwrap(phi));                 % Compute the empirical changes in phase from the data.
dM = diff(M);   
end

%% Animation

RUN = false; 
if RUN
clear 1; figure(1);
ax1 = subplot(121);
ax2 = subplot(122);
L = 200;
c = winter(L + 1);
for i = 1:2:length(dM)
m = max(1, i - L);
plot(ax1, M(1:end-1), dM, 'color', .75 * [1 1 1]); hold(ax1, 'on'); scatter(ax1, M(m:i), dM(m:i), 10, c(1:i - m + 1,:), 'filled'); set(ax1, 'xlim', [quantile(M, .01), quantile(M, .99)], 'ylim', [quantile(dM, .01), quantile(dM, .99)]); hold(ax1, 'off')
plot(ax2, (-499:500), data(inds(1:1000) - 500 + i), (-499:500), abs(data_analytic(inds(1:1000) - 500 + i))); xticks(0); ylim([min(ts_filt), max(ts_filt)]); grid on;
drawnow limitrate nocallbacks; 
end
end

%% Supplementary functions

function samples = SMC_update(samples, observation, dM, dphi, noisestd, dM2, dphi2)
	
	nsteps = length(dM);
	nparts = length(samples);
	eps = @() unidrnd(nsteps-1, 1, nparts);                 %   Get a random set of indices.
    
                                                        %   State transition function
    samples(1,:)=abs(samples(1,:)+dM(eps())');         %   ... perturb each amplitude by a random amount derived from the data.
    samples(2,:)=samples(2,:)+dphi(eps())';     %   ... perturb each phase by a random amount derived from the data.
	
	if exist('dM2', 'var')
		samples(3,:) = abs(samples(3,:) + dM2(eps())');
		samples(4,:) = samples(4,:) + dphi2(eps())';
	end

                                                        %   Determine the likelihood of each particle.
                                                        %   ... the estimated signal is amp*cos(phase). compute this for each particle and compare to signal.
                                                        %   ... assign likelihood using a Gaussian, with mean 0 for difference of particle and data, and fixed std.
                                                        %   ... add a fixed small probability to prevent particles from disappearing.
	est = samples(1,:).*cos(samples(2,:));
	if size(samples,1) == 4
		est = est + samples(3,:).*cos(samples(4,:));
	end
    p = normpdf(observation-est, 0 ,noisestd)+1e-6;
    p = p/sum(p);                                       %   Normalize the probability to sum to 1.

                                                        %   Resample using inverse CDF method. Make the CDF from p, choose a random CDF value [0 1],
                                                        %   ... then choose the particle with this CDF value. This is the posterior distribution.
    samples = samples(:,floor(interp1(cumsum(p),1:nparts,unifrnd(0,1,1,nparts),'linear',0))+1);
end
