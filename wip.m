% Phase tracker 

%% Load and filter the sample LFP data

PLOT = true; % debugging plots
DATA = 'LFP'; % 'LFP' for Alik/Ethan's sample data; 'HC' for Uri's

switch DATA
	% Each band should be in a separate row. For the LFP data, for example,
	% to look at signals from multiple bands, extract each signal
	% separately and store them in rows of ts_*.
	
	case 'LFP' % LFP data from Alik/Ethan
		load('/Volumes/NO NAME/DATA/LFPsamples/sample1.mat'); 
		k = 10;						% downsampling factor
		fs = 3e4 / k;				% sample rate
		low = 6; % [6 30];				% lower frequency bound(s)
		high = 10; % [10 40];				% upper frequency bound(s)

		sample1 = sample1(1:k:end);  % subsample
		sample1 = medfilt1(sample1); % median filter to reduce spurious spikes

		% Filter the signal
		ts_filt = zeros(length(low), length(sample1));
		ts_analytic = zeros(length(low), length(sample1));
		for i = 1:length(low)
			d = designfilt('bandpassiir','FilterOrder',4, ...
				'HalfPowerFrequency1',low(i),'HalfPowerFrequency2',high(i), ...
				'SampleRate',fs);						% Fourth order Butterworth
			ts_filt(i, :) = filtfilt(d, sample1);		% Get the filtered time series
			ts_analytic(i, :) = hilbert(ts_filt(i, :));	% Get the analytic signal
		end
		
	case 'HC' % hippocampal data from Uri
		load('/Volumes/NO NAME/GITHUB/Project-Instantaneous-Phase-Estimator/sample_data.mat') 
		ts_filt = EEGfilt(:)';              % ... and get the filtered data.
		ts_analytic = EEGanalytic(:)';
		fs = 3e3;
		k = 1;								% No downsampling
end		
t = (0:length(ts_filt)-1)/fs;				% time (for plotting)
nbands = size(ts_filt, 1);					% Number of frequency bands

%% Clip

% Select a chunk to analyze
inds = 10001:20000;                     % choose an interval of time,
ts_filt = ts_filt(:, inds);				% clip the filtered signal
ts_analytic = ts_analytic(:, inds);		% ... and the analytic signal
t = t(inds);							% ... and time

%% Get truth values
phi = unwrap(angle(ts_analytic), [], 2);	% calculate the phases
M = abs(ts_analytic);		% ... and amplitudes

if PLOT % Look at filtered and analytic signals
	figure(1); fullwidth(); ax1 = subplot(211); plot(t, ts_filt); hold on
	ax1.ColorOrderIndex = 1; plot(t, M); hold off;
	axis('tight'); title('Magnitude'); legend(cellstr(num2str((1:nbands)')));
	ax2 = subplot(212); plot(t, mod(phi, 2*pi)); linkaxes([ax1, ax2], 'x')
	title('Phase'); legend(cellstr(num2str((1:nbands)')));
	linkaxes([ax1 ax2], 'x')
end

%% Distribution of differences
% There is some tuning that needs to be done here. I'm not sure what makes
% the most sense. Originally, Uri scaled dM by a factor of 10 and left dphi
% as is. I thought it might also make sense to add some noise to the
% magnitude and phase traces and then take the differences... I dunno.
% Either way, the noise parameter for both dphi and dM needs to be
% optimized. Different sample rates and frequency bands change how well the
% PF does using a given set of noise parameters.

% dphi_FUN = @(phi, w) diff(unwrap(phi) + w * k^2 * (rand(size(phi)) - .5)); % Compute the empirical changes in phase from the data.
% dM_FUN = @(M, w) diff(M + w * k * (rand(size(M)) - .5));	% Compute the empirical changes in amplitude envelope from the data.
% dphi_FUN = @(phi, w) diff(unwrap(phi) + w/2 * k^2 * (randn(size(phi)))); % Compute the empirical changes in phase from the data.
% dM_FUN = @(M, w) diff(M + w/2 * k * (randn(size(M))));	% Compute the empirical changes in amplitude envelope from the data.
% dphi_FUN = @(phi, dphi, w) diff(unwrap(phi) + 1 * dphi(unidrnd(length(phi)-1, length(phi), 1))); % Compute the empirical changes in phase from the data.

% Add noise (chi2 and gaussian) to phases and magnitudes and calculate
% differences
dphi_FUN = @(phi, w) diff(phi + w(:) * ones(1, length(phi)) .* random('chi2', 1, size(phi)), [], 2);
dM_FUN = @(M, w) diff(M + w(:) * ones(1, length(phi)) .* randn(size(phi)), [], 2);

dphi = dphi_FUN(phi, 1e-3*sqrt(std(phi, [], 2)));
dM = dM_FUN(M, .1*sqrt(std(M, [], 2)));

if PLOT
	figure(11); fullwidth(); 
	for i = 1:nbands
		subplot(nbands, 2, 2*i - 1); 
		histogram(dM(i,:)); hold on; histogram(diff(M(i,:)')); hold off; 
		title(['M', num2str(i)])
		subplot(nbands, 2, 2*i); 
		histogram(dphi(i,:)); hold on; histogram(diff(phi(i,:)')); hold off; 
		title(['\phi', num2str(i)])
	end
end
drawnow;

% Run sequential monte carlo (SMC) a.k.a., particle filter (PF).

nsteps = length(inds);                                  % # time steps of data.
nparts = 1e4;											% # of particles
noisestd = 1;
samples = [median(dM, 2); median(dphi, 2)] * ones(1, nparts);		% Create initial set of particles, 
samples = [M(:, 1); phi(:, 2)] * ones(1, nparts);		% Create initial set of particles, 
														% they're all the same, at the 
														% first observed amp & phase.
sig = sum(ts_filt, 1) + noisestd * randn(size(ts_filt)); % ... add noise to observed signal, to vary difficulty of tracking.

Mest = M;
phiest = phi;


for i = 2:nsteps										% For each time point of data,
    
	samples = PF_update(...								% Update particles
		samples, sig(i), [dM; dphi], noisestd);			% Note: was 10*dM, 1.4*dphi
	Mest(:, i) = mean(samples(1:nbands,:), 2);          % Estimate amplitude from all particles.
    phiest(:, i) = angle(sum(exp(1i*samples((nbands+1):end,:)), 2)); % Estimate phase from all particles.
	
    if 0                                     %   Plot everything.
		if ~exist('pf', 'var'), pf = figure(99); fullwidth(); end
		figure(pf);
		for ii = 1:nbands
			subplot(1, nbands, ii);
			plot(samples(ii,:),samples(ii + nbands,:),'.',M(ii, i),phi(ii, i),'r*'); 
			axis([0 max(M(i, :)) 0 2*pi]); title(num2str(i));
			xlabel('Amplitude'); ylabel('Phase')
			hold on;
			plot(Mest(ii, i), phiest(ii, i), '*g')
			hold off
		end
		drawnow; pause(1e-2);
	end
end

figure(2); fullwidth(); plot(t,sig,t,sum(Mest.*cos(phiest), 1)); 
title('Signal'); legend('True', 'Estimated')

figure(3); fullwidth(); 
title('Components');
for i = 1:nbands
	plt = 2 * nbands * 100 + 10 + 2 * i - 1;
	ax1 = subplot(plt); plot(t, M(i, :), t, Mest(i, :)); 
	title(sprintf('M%d', i)); legend('True', 'Estimated')
	ax2 = subplot(plt + 1); plot(t, mod(phi(i, :), 2*pi), t, mod(phiest(i,:), 2*pi));
	title(sprintf('phi%d', i)); legend('True', 'Estimated')
	linkaxes([ax1, ax2], 'x');
end

%% Animation
% Don't remember what I was looking at but it was fun to look at so I'm not
% ready to throw it away yet. Probably doesn't workk right now, but i might
% still try to look at it again.
RUN = false; 
if RUN
clear 999; figure(999);
ax1 = subplot(121);
ax2 = subplot(122);
L = 200;
c = winter(L + 1);
for i = 1:2:length(dM)
m = max(1, i - L);
inds = max(i - 500, 1):(min(i + 500 - 1, length(dM)));
plot(ax1, M(1,1:end-1), dM, 'color', .75 * [1 1 1]); hold(ax1, 'on'); scatter(ax1, M(1,m:i), dM(1, m:i), 10, c(1:i - m + 1,:), 'filled'); set(ax1, 'xlim', [quantile(M(1,:), .01), quantile(M(1,:), .99)], 'ylim', [quantile(dM(1,:), .01), quantile(dM(1,:), .99)]); hold(ax1, 'off')
plot(ax2, ts_filt(1,inds)); hold on;
plot(ax2, abs(ts_analytic(1,inds))); 
xlim([-499 500]);xticks(0); ylim([min(ts_filt), max(ts_filt)]); grid on;
drawnow limitrate nocallbacks; 
end
end

%% Supplementary functions

function samples = PF_update(samples, observation, transfer, noisestd)
	
	nsteps = length(transfer);
	nparts = length(samples);
	nbands = size(transfer, 1) / 2;
								
	%   Get a random set of transfers
	deltas = cell2mat(arrayfun(@(row) transfer(row, unidrnd(nsteps-1, 1, nparts)), 1:size(samples,1), 'uniformoutput', false)');	
    
                                                    %   State transition function
	samples = samples + deltas;			%   ... perturb each particle
	samples(1:nbands,:) = abs(samples(1:nbands,:));	

													%   Determine the likelihood of each particle.
													%   ... the estimated signal is amp*cos(phase). compute this for each particle and compare to signal.
													%   ... assign likelihood using a Gaussian, with mean 0 for difference of particle and data, and fixed std.
													%   ... add a fixed small probability to prevent particles from disappearing.
	est = sum(samples(1:nbands,:).*cos(samples(nbands + 1:end,:)), 1);
	
	
	
    p = normpdf(observation-est, 0 ,noisestd)+1e-6;
    p = p/sum(p);                                       %   Normalize the probability to sum to 1.

                                                        %   Resample using inverse CDF method. Make the CDF from p, choose a random CDF value [0 1],
                                                        %   ... then choose the particle with this CDF value. This is the posterior distribution.
    samples = samples(:,floor(interp1(cumsum(p),1:nparts,unifrnd(0,1,1,nparts),'linear',0))+1);
end
