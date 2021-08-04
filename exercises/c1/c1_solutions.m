% 1.1, 1.2, 1.4, 1.8, 1.9

%% 1.1 
% Generate spikes for 10 s (or longer if you want better statistics)
% using a Poisson spike generator with a constant rate of 100 Hz, and
% record their times of occurrence. Compute the coefficient of variation of
% the interspike intervals, and the Fano factor for spike counts obtained
% over counting intervals ranging from 1 to 100 ms. Plot the interspike
% interval histogram.

t = 60;  % total time to generate spikes in seconds
dt = 0.001;  % bin size for generating spikes in seconds
rt = 100;  % firing rate in hz

% Generate spike process.
poisson_spk_process = poissrnd((rt * dt), 1, (t / dt));
% Can't have more than one spike in a bin.
spk_mask = poisson_spk_process > 0;
% Get spike times.
spk_times = find(spk_mask) / length(spk_mask) * t;
% Get ISIs and plot em.
isis = diff(spk_times);
figure, histogram(isis); xlabel('Time (s)'); ylabel('Counts'); 
title(sprintf('ISI Counts (total spikes = %i)', length(spk_times)));
% Get ISI coefficient of variation.
isi_cv = std(isis) / mean(isis);
% Get spike count fano factors.
spk_ct_bins = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1];  % in seconds
spk_ct_ffs = zeros(1, length(spk_ct_bins));
for j = 1 : length(spk_ct_bins)
    spk_ct_win = spk_ct_bins(j) / dt;  % spike count window
    spk_cts = movsum(spk_mask, spk_ct_win);  % spike counts
    spk_ct_ffs(j) = var(spk_cts) / mean(spk_cts);
end

% *Note*: Could have also generated the ISIs first (from an exponential
% distribution) and then found the spike times from these ISIs.

%% 1.2
% Add a refractory period to the Poisson spike generator by allowing the
% firing rate to depend on time (inhomogenous poisson process). Initially,
% set the firing rate to a constant value, r(t) = r0. After every spike,
% set r(t) to 0, and then allow it to recover exponentially back to r0 with
% a time constant T_ref that controls the refractory recovery rate. In
% other words, have r(t) obey the equation T_ref * dr/dt = r0 - r except
% immediately after a spike, when it is set to 0. Plot the coefficient of
% variation as a function of T_ref over the range 1 ms < T_ref < 20 ms, and
% plot interspike interval histograms for a few different values of T_ref
% in this range. Compute the Fano factor for spike counts obtained over
% counting intervals ranging from 1 to 100 ms for the case T_ref = 10 ms.

t = 60;  % total time to generate spikes, in seconds
dt = 0.001;  % bin size for generating spikes, in seconds
r0 = 100;  % firing rate in hz
spk_mask = zeros(1, (t / dt));  % boolean spike array, for every time bin.
% Time constants for refractory recovery rate, in seconds.
t_refs = linspace(0.001, 0.02, 5);
% Initialize cells for ISIs, coefficient of variations of ISIs, and spike
% count fann factors for each `t_ref`.
isis = cell(length(t_refs), 1);
isis_hist_fig = figure;
isis_cv = cell(length(t_refs), 1);
spk_ct_ffs = cell(length(t_refs), 1);

for i = 1:length(t_refs)
    % Get current `t_ref`.
    t_ref = t_refs(i);
    % The difference from the last bin in which there was a spike to the 
    % current bin (initialize high).
    last_spk_bin = length(spk_mask);
    % Initialize firing rate.
    r = r0;
    % For each possible time bin, see if a spike occurred.
    for j = 2 : (t/dt)
        % Refractory penalty: recently occuring spikes reduce firing rate.
        ref_pen = r0 / exp(last_spk_bin * dt / t_ref);
        % Compute current firing rate.
        r = r0 - ref_pen;
        spk_mask(j) = poissrnd((r * dt));
        if spk_mask(j) > 0
            last_spk_bin = 0;
        else
            last_spk_bin = last_spk_bin + 1;
        end
    end
    % Get spike times.
    spk_times = find(spk_mask) * dt;
    isis{i} = diff(spk_times);
    isis_hist_fig; hold on; histogram(isis{i});
    isis_cv{i} = std(isis{i}) / mean(isis{i});
    % Get spike count fano factors.
    spk_ct_bins = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1];  % in secs
    spk_ct_ffs_cur = zeros(1, length(spk_ct_bins));
    for k = 1 : length(spk_ct_bins)
        spk_ct_win = spk_ct_bins(k) / dt;  % window to sum spike counts
        spk_cts = movsum(spk_mask, spk_ct_win);  % spike counts
        spk_ct_ffs_cur(k) = var(spk_cts) / mean(spk_cts);
    end
    spk_ct_ffs{j} = spk_ct_ffs_cur;
end
isis_hist_fig; 
legend(num2str(t_refs(1)), num2str(t_refs(2)), num2str(t_refs(3)), ...
       num2str(t_refs(4)), num2str(t_refs(5)));
xlabel('Time (s)'); ylabel('counts');
title(sprintf(['ISI Histograms for Inhomogenous Poisson Processes \n' ...
               'with Different Time Constants for Refractory Recovery ' ...
               'Rate']));

%% 1.4
% Generate a Poisson spike train with a time-dependent firing rate
% r(t) = 100(1 + cos(2*pi*t/300 ms)) Hz. Approximate the firing rate from
% this spike train using a variable r_approx that satisfies: 
% T_approx * drapprox/dt = -r_approx, except that 
% rapprox -> rapprox + 1/T_approx every time a spike occurs. Make plots of 
% the true rate, the spike sequence generated, and the estimated rate. 
% Experiment with a few different values of T_approx in the range of 1 to 
% 100 ms. Determine the best value of T_approx by computing the average
% squared error of the estimate, the integral of: 
% dt(r(t) - rapprox(t))^2 , for different values of T_approx, and finding
% the value of T_approx that minimizes this error.

% r(t) = 100 * (1 + cos(2*pi*t/.3))
% dr/dt = 100 * ((6.67 * pi) * (-sin(6.67 * pi * t)))
%       ~= 2095 * -sin(20.95 * t)

t = 60;  % total time to generate spikes, in seconds
dt = 0.001;  % bin size for generating spikes, in seconds
% The actual firing rate.
r = 100 * (1 + cos(2 * pi * ([1:(t / dt)] * dt)  / .3));
r_spk_mask = poissrnd((r * dt), 1, (t / dt));
r_spk_times = find(r_spk_mask) * dt;

% Possible `t_approxs` to use to approximate true firing rate.
t_approxs = linspace(0.001, 0.1, 10);
spk_mask = zeros(1, (t / dt));  % boolean spike array, for every time bin.

% Initialize the firing rate estimates and mean-squared-errors from true
% firing rate for each `t_approx`
r_ests = cell(1, length(t_approxs));
r_errs = zeros(1, length(t_approxs));

% Compute spike processes for each `t_approx`.
for i = 1:length(t_approxs)
    % Get current `t_approx`.
    t_approx = t_approxs(i);
    % Current estimated firing rate computed with `t_approx`.
    r_est = -(t_approx * (2095 * -sin(20.95 * ([1:(t / dt)] * dt)))) ... 
            + (t_approx * length(r_spk_times));
    r_ests{i} = r_est;
%     % For each possible time bin, see if a spike occurred.
%     for j = 2 : (t / dt)
%         % Compute current firing rate.
%         spk_mask(j) = poissrnd((r_est * dt));
%         % Check whether to update `r_est` if there was a spike.
%         if spk_mask(j)
%             r_est = r_est + t_approx;
%         end
%         % Don't allow for refractory violations.
%         if spk_mask(j) && spk_mask(j - 1) && (dt < 0.002)
%             spk_mask(j) = 0;
%         end
%     end
%     % Get spike times.
%     spk_times = find(spk_mask) * dt;
    % Plot spike times, `r`, and `r_est`.
    figure; hold on;
    scatter(r_spk_times, (mean(r) * ones(1, length(r_spk_times))));
    plot((dt * [1 : (t / dt)]), r);
    plot((dt * [1 : (t / dt)]), r_est);
    xlim([0, 2]); xlabel('Time (s)'); ylabel('Firing rate (spks/s)');
    title(sprintf('t_approx = %d', t_approx));
    % Compute the mean squared error between `r` and `r_est`.
    r_err = sum((r - r_est).^2);
    r_errs(i) = r_err;
end

% Find the best `t_approx`: the `t_approx` that gives lowest `r_err`.
[~, best_t_approx_idx] = (min(r_errs));
best_t_approx = t_approxs(best_t_approx_idx);


%% 1.8
% File c1p8.mat contains data collected and provided by Rob de Ruyter van
% Steveninck from a fly H1 neuron responding to an approximate white-noise
% visual motion stimulus. Data were collected for 20 minutes at a sampling
% rate of 500 Hz. In the file, rho is a vector that gives the sequence of
% spiking events or nonevents at the sampled times (every 2 ms). When an
% element of rho is one, this indicates the presence of a spike at the
% corresponding time, whereas a zero value indicates no spike. The variable
% stim gives the sequence of stimulus values at the sampled times (ratio of
% black-to-white pixels). Calculate and plot the spike-triggered average
% from these data over the range from 0 to 300 ms (150 time steps).

load('c1p8.mat');

dt = 0.002;  % sampling time in s
s_t_a_t = 0.3;  % STA time range in s
% Build stim design matrix. 
% number of params in `x` design matrix (bins in past)
n_p_x = s_t_a_t / dt;
% Pad early bins of stimulus with zero.
padded_stim = [zeros(n_p_x - 1, 1); stim];
x = hankel(padded_stim(1 : (end - n_p_x + 1)), ...
           padded_stim((end - n_p_x + 1) : end));

% Let's visualize a small part of the design matrix just to see it.
n_obs_disp = 100;  % number of observations to display
% Display an `n_r_x` by `n_obs_disp` portion of `x`.
figure;
imagesc(-n_p_x + 1 : 0, 1 : n_obs_disp, x(1 : n_obs_disp, :));
h_cb = colorbar;
colormap(gray)
h_cb.Label.String = 'stim intensity values';
xlabel('lags before spike time bin');
ylabel('time bin of response');
title('design matrix');

% Compute and plot STA.
n_spks = length(find(rho));
s_t_a = (x' * rho) / n_spks;
s_t_a_t_bins = [(-n_p_x + 1) : 1 : 0] * dt;
figure; plot(s_t_a_t_bins, s_t_a);
title('STA');
xlabel('time before spike (s)');
ylabel('stim intensity');


%% 1.9 
% Using the data of 1.8, calculate and plot stimulus averages triggered on
% events consisting of a pair of spikes (which need not necessarily be
% adjacent) separated by a given interval (as in figure 1.10). Plot these
% two-spike-triggered average stimuli for various separation intervals
% ranging from 2 to 100 ms. (Hint: use convolution for pattern matching:
% e.g. `find(conv(rho,[1 0 1]) == 2)` will contain the indices of all the
% events with two spikes separated by 4 ms.) Plot, as a function of the
% separation between the two spikes, the magnitude of the difference
% between the two-spike-triggered average and the sum of two
% single-spike-triggered averages (obtained in 1.8) separated by the same
% time interval. At what temporal separations is this difference
% negligible? (i.e. at what separations does the two-spike STA *not* tell 
% us anything more than a single-spike STA?)

% separation intervals in which to look for two spikes
sep_intrvls = linspace(0.002, 0.1, 10);
% array to hold two-spike-triggered-averages for each separation interval
two_spk_s_t_as = zeros(length(sep_intrvls), n_p_x);

stas_fig = figure; hold on;
stas_diff_fig = figure; hold on;

% Compute two-spike STA for each separation interval
for i = 1 : length(sep_intrvls)
    sep = sep_intrvls(i);
    % Get filter for current separation interval
    two_spk_filt = ones(1, (1 + ceil(sep / dt)));
    rho2 = conv(rho, two_spk_filt, 'same') == 2;
    n_evnts = length(find(rho2));
    % Compute two-spike STA.
    s_t_a2 = (x' * rho2) / n_evnts;
    two_spk_s_t_as(i, :) = s_t_a2;
    % Plot two-spike STA.
    figure(stas_fig); plot(s_t_a_t_bins, s_t_a2);
    % Compute magnitude of diff between two-spike STA and sum of 2
    % one-spike STA.
    s_t_a_diff = sum(((s_t_a * 2) - s_t_a2).^2);
    figure(stas_diff_fig); scatter(sep, s_t_a_diff, 'b');
end

figure(stas_fig);
legend(num2str(sep_intrvls(1)), num2str(sep_intrvls(2)), ...
       num2str(sep_intrvls(3)), num2str(sep_intrvls(4)), ...
       num2str(sep_intrvls(5)), num2str(sep_intrvls(6)), ...
       num2str(sep_intrvls(7)), num2str(sep_intrvls(8)), ...
       num2str(sep_intrvls(9)), num2str(sep_intrvls(10)));
title('Two-spike STAs for different separation intervals');
xlabel('time before spike (s)');
ylabel('stim intensity');

figure(stas_diff_fig);
xlabel('two-spike separation interval');
ylabel('Sum of Squared Diffs: 2 vs. 1 spike STA');
title('Sum of Squared Diffs: 2 vs. 1 spike STA');
