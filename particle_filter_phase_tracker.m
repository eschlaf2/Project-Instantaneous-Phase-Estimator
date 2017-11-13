%% Load the data and select a chunk to analyze.

clear
load('sample_data.mat')                                                         % Load hippocampal data to get a representative slow rhythm.

inds = 10001:20000;                                                             % choose an interval of time,
x = EEGfilt(inds);                                                              % ... and get the filtered data.
M = abs(EEGanalytic(inds));                                                     % Compute the amplitude envelope.
phi = mod(angle(EEGanalytic(inds)),2*pi);                                       % Compute the phase.;

dphi = diff(unwrap(phi));                                                       % Compute the empirical changes in phase from the data.
dM = diff(M);                                                                   % Compute the empirical changes in amplitude envelope from the data.

%% Run sequencial monte carlo (SMC) a.k.a., particle filter.

nparts = 10000;                                         % # particles.
nsteps = length(inds);                                  % # time steps of data.
samples = [M(1); phi(1)]*ones(1,nparts);                % Create initial set of particles, they're all the same, at the first observed amp & phase.
noisestd = 10;                                          % Noise level to add to our observed signal,
sig = x+normrnd(0,noisestd,size(x));                    % ... add noise to observed signal, to vary difficulty of tracking.
for i = 1:length(inds),                                 % For each time point of data,
    eps = unidrnd(nsteps-1, 1, nparts);                 %   Get a random set of indices.
    
                                                        %   State transition function
    samples(1,:)=abs(samples(1,:)+10*dM(eps)');         %   ... pertrub each amplitude by a random amount derived from the data.  [? Don't need "abs" ?]
    samples(2,:)=mod(samples(2,:)+dphi(eps)',2*pi);     %   ... pertrub each phase by a random amount derived from the data.

                                                        %   Determine the likelihood of each particle.
                                                        %   ... the estimated signal is amp*cos(phase). compute this for each particle and compare to signal.
                                                        %   ... assign likelihood using a Gaussian, with mean 0 for difference of particle and data, and fixed std.
                                                        %   ... add a fixed small probability to prevent particles from disappearing.
    p = normpdf(sig(i)-samples(1,:).*cos(samples(2,:)), 0 ,noisestd)+1e-6;
    p = p/sum(p);                                       %   Normalize the probability to sum to 1.

                                                        %   Resample using inverse CDF method. Make the CDF from p, choose a random CDF value [0 1],
                                                        %   ... then choose the particle with this CDF value. This is the posterior distribution.
    samples = samples(:,floor(interp1(cumsum(p),1:nparts,unifrnd(0,1,1,nparts),'linear',0))+1);
                                                        %   Plot everything.
    plot(samples(1,:),samples(2,:),'.',M(i),phi(i),'ro'); axis([0 300 0 2*pi]); title(num2str(i));
    xlabel('Amplitude'); ylabel('Phase')
    Mest(i) = mean(samples(1,:));                       %   Estimate amplitude from all particles.
    phiest(i) = angle(sum(exp(1i*samples(2,:))));       %   Estimate phase from all particles.
    hold on;
    plot(Mest(i), phiest(i), 'xk')
    hold off
    drawnow; pause(.001);
end;
plot(inds,sig,inds,Mest.*cos(phiest));

