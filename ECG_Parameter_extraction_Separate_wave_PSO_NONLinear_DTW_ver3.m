clear                     % Remove all variables from workspace
close all                 % Close all open figure windows

% Load ECG data from MAT file
[file, path] = uigetfile('*.mat','Select ECG mat file');
data = load(file);

L_num_of_Gaussian_kernels_P_wave = 3;
L_num_of_Gaussian_kernels_QRS_wave = 5;
L_num_of_Gaussian_kernels_T_wave = 3;


fs = data.fs;                % Sampling frequency of ECG signal (Hz)

x = data.x;               % Extract stored signal matrix
ecg = x(1,1:end);             % Use first channel as ECG signal
length_sig = length(ecg); % Total number of ECG samples

ecg_bins = 200;           % Number of phase bins for ECG mean calculation
Length_of_Corr_window_for_whole_ecg = round(fs*0.1);
Corr_window_for_whole_ecg  = triang(Length_of_Corr_window_for_whole_ecg);
Length_of_Corr_window_for_ecg_mean = round(fs*0.1*(ecg_bins/fs));
Corr_window_for_ecg_mean  = triang(Length_of_Corr_window_for_ecg_mean);

%% 'Wavelet from Triangular Scaling Function'
% t = linspace(-1,2,Length_of_Corr_window_for_whole_ecg);
% 
% Corr_window_for_whole_ecg = psi_triangular(t);
% 
% t = linspace(-1,2,Length_of_Corr_window_for_ecg_mean);
% 
% Corr_window_for_ecg_mean = psi_triangular(t);



%% -------- R-peak detection using Pan–Tompkins algorithm
[qrs_positions] = pantompkins_qrs(abs(ecg),fs);

figure(1),plot(ecg,'b'),hold on,plot(qrs_positions,ecg(qrs_positions),'*r'),hold off
legend({'ECG Signal','R Peaks'})
title(file)
%% -------- Phase calculation
% Linear phase based on RR intervals
[Linearphase,~] = calculate_linear_phase_ver2(qrs_positions,length_sig,fs);

% Alternative options:
% NonlinearPhase = calculate_dtw_phase(ecg, qrs_positions);   % DTW-based nonlinear phase

NonlinearPhase = Linearphase; 

ecg_beat_cell = {};
counter = 0;
for k=1:length(qrs_positions)-1
    if (qrs_positions(k+1)-qrs_positions(k))>4
        counter = counter +1;
        
        ecg_beat_cell{counter,1} = [ecg(:,qrs_positions(k):qrs_positions(k+1));Linearphase(:,qrs_positions(k):qrs_positions(k+1))];

        sample_index = qrs_positions(k):qrs_positions(k+1);
        ecg_beat_cell{counter,2} = sample_index;
    end
end
%===============nonlinear phase assignment
a= ecg_beat_cell{3,1}; % reference beat

for i=1:size(ecg_beat_cell,1)
    b= ecg_beat_cell{i,1}; % test beat
 

    phase_a = a(2,:)';
    phase_b = b(2,:)';
 

% matlab 2015==============

% [DTW_matrix,Dist,optimal_path]=dtw_ver1((a(1,:))',(b(1,:))');
% phase_b(optimal_path(:,2)) = phase_a(optimal_path(:,1));
% matlab 2015==============

% matlab 2017==============


% [d,i1,i2] = dtw((a(1,:)),(b(1,:)),'squared',round(fs/4));
[d,i1,i2] = dtw((a(1,:)),(b(1,:)));

phase_b(i2) = phase_a(i1);
%=======================
 

 sample_index = ecg_beat_cell{i,2};

 NonlinearPhase(1,sample_index) = phase_b;
%    str = ['beats' '   ' num2str(i) '/' num2str(size(ecg_beat_cell,1))];
%    disp(str)
end

% NonlinearPhase = Linearphase;




%% -------- ECG mean and standard deviation in phase domain
[~,ECGmean,~] = ecgsd_extractor_ver1(ecg,Linearphase,ecg_bins);
[ECGsd,ECGmean_nonlinear_phase,meanphase] = ecgsd_extractor_ver1(ecg,NonlinearPhase,ecg_bins);
% ECGmean_nonlinear_phase =  wdenoise(ECGmean_nonlinear_phase,7,Wavelet="bior4.4",DenoisingMethod="SURE",NoiseEstimate="LevelDependent");
figure(2),plot(ECGmean,'r'),hold on,plot(ECGmean_nonlinear_phase,'b'),hold off
title(file)
legend({'ECG mean Linear','ECG mean nonlinear'})

%% Finding best Gaussian Kernels for P wave in ECG mean

P_wave_ECG_mean = (meanphase<-pi/6 & meanphase>-pi/6-pi/2).*ECGmean_nonlinear_phase;
% further smoothing of ECG mean using wavelet
% P_wave_ECG_mean = wdenoise(P_wave_ECG_mean,7,Wavelet="bior4.4",DenoisingMethod="SURE",NoiseEstimate="LevelDependent");

P_wave_ECG_mean_smooth = sgolayfilt(P_wave_ECG_mean,3,11);

%% 'Wavelet from Triangular Scaling Function'

% t = linspace(-1,2,Length_of_Corr_window_for_ecg_mean);
% 
% Corr_window_for_ecg_mean = psi_triangular(t);



P_mask = (meanphase<-pi/6 & meanphase>-pi/6-pi/2);
% Instead of findpeaks(corr_P_win), use abs after removing baseline

corr_P_win = conv(P_wave_ECG_mean_smooth,Corr_window_for_ecg_mean,"same");
corr_P_abs = abs(corr_P_win - medfilt1(corr_P_win,round(fs/10)));

[~,idx_P_peak_ecg_mean] = max(corr_P_abs);
figure(3),plot(P_wave_ECG_mean),hold on, plot(idx_P_peak_ecg_mean,P_wave_ECG_mean(idx_P_peak_ecg_mean),'*'),plot(ECGmean_nonlinear_phase,'k'),hold off
title('P peak of ECG mean')
phase_P = meanphase(idx_P_peak_ecg_mean);
phase_P = min(max(phase_P,-pi/2),0);
%% -------- P wave parameter extraction using Gaussian mixture model

% ========================building new myfun based on L Gaussians
ecg_mean_temp = 0;
ai_total = [];
bi_total = [];
tetai_total  = [];
ai_Pwave = [];
bi_Pwave = [];
tetai_Pwave  = [];
for i=1:L_num_of_Gaussian_kernels_P_wave
% disp(num2str(i))
ecg_mean_temp1 = P_wave_ECG_mean - ecg_mean_temp;
 
lb = [-1.5*max(abs(ecg)).*ones(1,1)   0.000001*ones(1,1)   max(phase_P-0.8,-pi)*ones(1,1)  ];

ub = [1.5*max(abs(ecg)).*ones(1,1)  5*ones(1,1)  min(phase_P+0.8,0)  ];
myfun1 = @(params)  norm(ecg_mean_temp1'-sum((repmat(params(1:1),ecg_bins,1).*exp(-(rem(repmat(meanphase,1,1)'-repmat(params(3),ecg_bins,1)+pi,2*pi)-pi) .^2 ./ (2*(repmat(params(2),ecg_bins,1)) .^ 2))),2));


% options = optimoptions('particleswarm','SwarmSize',30,'HybridFcn',@fmincon,'MaxIter',1000);
options = optimoptions('particleswarm','SwarmSize',2000,'MaxIter',200,'Display','off');

OptimumParams = particleswarm(myfun1,3*1,lb,ub,options);

% L = (length(OptimumParams)/3);
if isempty(OptimumParams)
    break
end
ai_1 = OptimumParams(1);
bi_1 = OptimumParams(2);
tetai_1 = OptimumParams(3);

ai_Pwave = [ai_Pwave ai_1];
bi_Pwave = [bi_Pwave bi_1];
tetai_Pwave  = [tetai_Pwave tetai_1];

ai_total = [ai_total ai_1];
bi_total = [bi_total bi_1];
tetai_total  = [tetai_total tetai_1];
dtetai_1 = rem(meanphase - tetai_1 + pi,2*pi)-pi;
ecg_mean_temp = ecg_mean_temp + ai_1 .* exp(-dtetai_1 .^2 ./ (2*bi_1 .^ 2));
figure(39),plot(ecg_mean_temp,'b'),hold on,plot(P_wave_ECG_mean,'r')
legend({'Synthetic P wave','P wave Mean'}),hold off
title(['Synthetic P wave   ' num2str(i) 'th' '  Gaussian'])
% pause(3)
end
%% -------- T wave parameter extraction using Gaussian mixture model

T_wave_ECG_mean = (meanphase>pi/6 & meanphase<pi/6+3*pi/4).*ECGmean_nonlinear_phase;



T_mask = (meanphase>pi/6 & meanphase<pi/6+2*pi/3);
% Instead of findpeaks(corr_T_win), use abs after removing baseline

corr_T_win = conv(T_wave_ECG_mean,Corr_window_for_ecg_mean,"same");
corr_T_abs = abs(corr_T_win - median(corr_T_win(T_mask)));

[~,idx_T_peak_ecg_mean] = max(corr_T_abs);
figure,plot(T_wave_ECG_mean),hold on, plot(idx_T_peak_ecg_mean,T_wave_ECG_mean(idx_T_peak_ecg_mean),'*'),hold off
phase_T = meanphase(idx_T_peak_ecg_mean);
title('T peak of ECG mean')

phase_T = min(max(phase_T,pi/6),2*pi/3);


% ========================building new myfun based on L Gaussians

ecg_mean_temp = 0;
ai_Twave = [];
bi_Twave = [];
tetai_Twave  = [];
for i=1:L_num_of_Gaussian_kernels_T_wave
% disp(num2str(i))
ecg_mean_temp1 = T_wave_ECG_mean - ecg_mean_temp;
lb = [-1.5*max(abs(ecg)).*ones(1,1)   0.000001*ones(1,1)   max(phase_T-0.5,pi/6)*ones(1,1)  ];
ub = [(1.5*max(abs(ecg))).*ones(1,1)  5*ones(1,1)  min(phase_T+0.5,3*pi/4)*ones(1,1)  ];
myfun1 = @(params)  norm(ecg_mean_temp1'-sum((repmat(params(1:1),ecg_bins,1).*exp(-(rem(repmat(meanphase,1,1)'-repmat(params(3),ecg_bins,1)+pi,2*pi)-pi) .^2 ./ (2*(repmat(params(2),ecg_bins,1)) .^ 2))),2));


% options = optimoptions('particleswarm','SwarmSize',30,'HybridFcn',@fmincon,'MaxIter',1000);
options = optimoptions('particleswarm','SwarmSize',1000,'MaxIter',200,'Display','off');

OptimumParams = particleswarm(myfun1,3*1,lb,ub,options);

% L = (length(OptimumParams)/3);
if isempty(OptimumParams)
    break
end
ai_1 = OptimumParams(1);
bi_1 = OptimumParams(2);
tetai_1 = OptimumParams(3);

ai_Twave = [ai_Twave ai_1];
bi_Twave = [bi_Twave bi_1];
tetai_Twave  = [tetai_Twave tetai_1];

ai_total = [ai_total ai_1];
bi_total = [bi_total bi_1];
tetai_total  = [tetai_total tetai_1];
dtetai_1 = rem(meanphase - tetai_1 + pi,2*pi)-pi;
ecg_mean_temp = ecg_mean_temp + ai_1 .* exp(-dtetai_1 .^2 ./ (2*bi_1 .^ 2));
figure(41),plot(ecg_mean_temp,'b'),hold on,plot(T_wave_ECG_mean,'r')
legend({'Synthetic T wave','T wave Mean'}),hold off
title(['Synthetic T wave  ' num2str(i) 'th' '  Gaussian'])
% pause(3)
end



%% -------- QRS wave parameter extraction using Gaussian mixture model

QRS_wave_ECG_mean = (meanphase>-pi/6 & meanphase<pi/6).*ECGmean_nonlinear_phase;
 
QRS_mask = (meanphase>-pi/6 & meanphase<pi/6);
% Instead of findpeaks(corr_T_win), use abs after removing baseline

corr_QRS_win = conv(QRS_wave_ECG_mean,Corr_window_for_ecg_mean,"same");
corr_QRS_abs = abs(corr_QRS_win - median(corr_QRS_win(QRS_mask)));

[~,idx_QRS_peak_ecg_mean] = max(corr_QRS_abs);
figure,plot(QRS_wave_ECG_mean),hold on, plot(idx_QRS_peak_ecg_mean,QRS_wave_ECG_mean(idx_QRS_peak_ecg_mean),'*'),hold off
phase_QRS = meanphase(idx_QRS_peak_ecg_mean);
title('QRS peak of ECG mean')

phase_QRS = max(min(phase_QRS,pi/6),-pi/6);


% ========================building new myfun based on L Gaussians
ecg_mean_temp = 0;

ai_QRSwave = [ ];
bi_QRSwave = [ ];
tetai_QRSwave  = [ ];


for i=1:L_num_of_Gaussian_kernels_QRS_wave
% disp(num2str(i))
ecg_mean_temp1 = QRS_wave_ECG_mean - ecg_mean_temp;
lb = [-1.5*max(abs(ecg)).*ones(1,1)   0.000001*ones(1,1)   max(phase_QRS-0.5,-pi/3)*ones(1,1)  ];
ub = [1.5*max(abs(ecg)).*ones(1,1)  5*ones(1,1)  min(phase_QRS+0.5,pi/3)*ones(1,1)  ];
myfun1 = @(params)  norm(ecg_mean_temp1'-sum((repmat(params(1:1),ecg_bins,1).*exp(-(rem(repmat(meanphase,1,1)'-repmat(params(3),ecg_bins,1)+pi,2*pi)-pi) .^2 ./ (2*(repmat(params(2),ecg_bins,1)) .^ 2))),2));


% options = optimoptions('particleswarm','SwarmSize',30,'HybridFcn',@fmincon,'MaxIter',1000);
options = optimoptions('particleswarm','SwarmSize',500,'MaxIter',200,'Display','off');

OptimumParams = particleswarm(myfun1,3*1,lb,ub,options);

% L = (length(OptimumParams)/3);
if isempty(OptimumParams)
    break
end
ai_1 = OptimumParams(1);
bi_1 = OptimumParams(2);
tetai_1 = OptimumParams(3);

ai_QRSwave = [ai_QRSwave ai_1];
bi_QRSwave = [bi_QRSwave bi_1];
tetai_QRSwave  = [tetai_QRSwave tetai_1];

ai_total = [ai_total ai_1];
bi_total = [bi_total bi_1];
tetai_total  = [tetai_total tetai_1];

dtetai_1 = rem(meanphase - tetai_1 + pi,2*pi)-pi;
ecg_mean_temp = ecg_mean_temp + ai_1 .* exp(-dtetai_1 .^2 ./ (2*bi_1 .^ 2));
figure(42),plot(ecg_mean_temp,'b'),hold on,plot(QRS_wave_ECG_mean,'r')
legend({'Synthetic QRS wave','QRS wave Mean'}),hold off
title(['Synthetic QRS wave   ' num2str(i) 'th' '  Gaussian'])
% pause(3)
end

Synthetic_P = 0;
for j=1:length(ai_Pwave)
    dtetai = rem(NonlinearPhase - tetai_Pwave(j) + pi,2*pi)-pi;

    Synthetic_P = Synthetic_P+ai_Pwave(j) .* exp(-dtetai .^2 ./ (2*bi_Pwave(j) .^ 2));
end

figure(48);
plot(Synthetic_P,'b','LineWidth',2), hold on
plot(ecg,'--r')
iptsetpref('ImshowBorder','tight')
legend({'Synthetic Nonlinear P wave','Original ECG'})
hold off
title(file)

Synthetic_T = 0;
for j=1:length(ai_Twave)
    dtetai = rem(NonlinearPhase - tetai_Twave(j) + pi,2*pi)-pi;

    Synthetic_T = Synthetic_T+ai_Twave(j) .* exp(-dtetai .^2 ./ (2*bi_Twave(j) .^ 2));
end

figure(49);
plot(Synthetic_T,'b','LineWidth',2), hold on
plot(ecg,'--r')
iptsetpref('ImshowBorder','tight')
legend({'Synthetic Nonlinear T wave','Original ECG'})
hold off
title(file)


%% %% ُsynthetic QRS wave
Synthetic_QRS = 0;
for j=1:length(ai_QRSwave)
    dtetai = rem(NonlinearPhase - tetai_QRSwave(j) + pi,2*pi)-pi;

    Synthetic_QRS = Synthetic_QRS+ai_QRSwave(j) .* exp(-dtetai .^2 ./ (2*bi_QRSwave(j) .^ 2));
end

figure(51);
plot(Synthetic_QRS,'b','LineWidth',2), hold on
plot(ecg,'--r')
iptsetpref('ImshowBorder','tight')
legend({'Synthetic QRS wave','Original ECG'})
hold off
title(file)








Alpha_i = ai_total;   % Gaussian amplitudes
Beta_i  = bi_total;   % Gaussian widths
Theta_i = tetai_total;   % Gaussian centers
[Theta_i,idx] = sort(Theta_i,'ascend');
Alpha_i = Alpha_i(idx);
Beta_i = Beta_i(idx);
OptimumParams = [Alpha_i Beta_i Theta_i];

Synthetic_ECG_mean = 0;
for j=1:length(ai_total)
    dtetai = rem(meanphase - tetai_total(j) + pi,2*pi)-pi;

    Synthetic_ECG_mean = Synthetic_ECG_mean+ai_total(j) .* exp(-dtetai .^2 ./ (2*bi_total(j) .^ 2));
end

%% -------- Plot ECG mean vs synthetic ECG mean
figure(46);
plot(Synthetic_ECG_mean,'b','LineWidth',2), hold on
plot(ECGmean_nonlinear_phase,'--r')
iptsetpref('ImshowBorder','tight')
legend({'Synthetic ECG mean','Original ECG mean Nonlinear'})
hold off
title(file)

%% ploting synthetic ECG

Synthetic_ECG = 0;
for j=1:length(ai_total)
    dtetai = rem(NonlinearPhase - tetai_total(j) + pi,2*pi)-pi;

    Synthetic_ECG = Synthetic_ECG+ai_total(j) .* exp(-dtetai .^2 ./ (2*bi_total(j) .^ 2));
end

figure(47);
plot(Synthetic_ECG,'b','LineWidth',2), hold on
plot(ecg,'--r')
iptsetpref('ImshowBorder','tight')
legend({'Synthetic Nonlinear ECG','Original ECG'})
hold off
title(file)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions

function [ecgsd,ecg_mean,phase_mean] = ecgsd_extractor_ver1(ecg,phase,bins)

x1 = ecg;                        % ECG signal
meanPhase = zeros(1,bins);       % Mean phase per bin
ECGmean = zeros(1,bins);         % Mean ECG per bin
ECGsd = zeros(1,bins);           % ECG standard deviation per bin

% Handle wrap-around phase bin near -pi / +pi
I = find( phase >= (pi-pi/bins) | phase < (-pi+pi/bins) );

if(~isempty(I))
    meanPhase(1) = -pi;
    ECGmean(1) = mean(x1(I));
    ECGsd(1) = std(x1(I));
else
    ECGsd(1) = -1;               % Mark empty bins
end

% Loop over phase bins
for i = 1 : bins-1
    I = find( phase >= 2*pi*(i-0.5)/bins - pi & ...
        phase <  2*pi*(i+0.5)/bins - pi );

    if(~isempty(I))
        meanPhase(i+1) = mean(phase(I));
        ECGmean(i+1) = mean(x1(I));
        ECGsd(i+1) = std(x1(I));
    else
        ECGsd(i+1) = -1;
    end
end

% Interpolate missing bins
K = find(ECGsd==-1);

for i = 1:length(K)
    switch K(i)
        case 1
            meanPhase(1) = -pi;
            ECGmean(1) = ECGmean(2);
            ECGsd(1) = ECGsd(2);
        case bins
            meanPhase(bins) = pi;
            ECGmean(bins) = ECGmean(bins-1);
            ECGsd(bins) = ECGsd(bins-1);
        otherwise
            meanPhase(K(i)) = mean(meanPhase([K(i)-1 K(i)+1]));
            ECGmean(K(i))   = mean(ECGmean([K(i)-1 K(i)+1]));
            ECGsd(K(i))     = mean(ECGsd([K(i)-1 K(i)+1]));
    end
end

phase_mean = meanPhase;
ecg_mean   = ECGmean;
ecgsd      = ECGsd;

end



function [Phase,Omega] = calculate_linear_phase_ver2(locs,length_sig,fs)

% locs       : qrs_positionss of detected R-peaks
% length_sig : total number of ECG samples
% fs         : sampling frequency

ind = locs(:)';                % Convert R‑peak qrs_positionss to row vector

Phase = zeros(1,length_sig);   % Phase of each ECG sample
Omega = zeros(1,length_sig);   % Instantaneous angular frequency

RR = mean(diff(ind));          % Mean RR interval (samples)

%% -------- Phase before the first R‑peak

stepTheta = 2*pi/RR;           % Average phase increment per sample
omega_val = fs*stepTheta;      % Instantaneous angular frequency

theta = 0;                     % Initialize phase

for j = ind(1)-1:-1:1          % Move backward from first R‑peak
    theta = theta - stepTheta; % Decrease phase
    theta = mod(theta+pi,2*pi)-pi; % Wrap phase into [-pi , pi]

    Phase(j) = theta;          % Store phase
    Omega(j) = omega_val;      % Store frequency
end

%% -------- Phase between consecutive R‑peaks

for k = 1:length(ind)-1

    bins = ind(k+1)-ind(k);    % Number of samples between R-peaks

    stepTheta = 2*pi/bins;     % Phase increment so phase spans one cycle
    omega_val = fs*stepTheta;  % Corresponding angular frequency

    theta = 0;
    Phase(ind(k)) = 0;         % Define phase at R‑peak as zero

    for j = ind(k)+1 : ind(k+1)-1
        theta = theta + stepTheta; % Linear phase progression
        if theta>pi
            theta = -pi;
        end
        Phase(j) = theta;
        Omega(j) = omega_val;
    end

    Phase(ind(k+1)) = 0;       % Next R‑peak also set to zero phase
end

%% -------- Phase after the last R‑peak

stepTheta = 2*pi/RR;           % Use mean RR again
omega_val = fs*stepTheta;

theta = 0;

for j = ind(end)+1:length_sig
    theta = theta + stepTheta; % Continue phase linearly
    theta = mod(theta+pi,2*pi)-pi;

    Phase(j) = theta;
    Omega(j) = omega_val;
end

end



function [Phase, Omega] = calculate_ECG_phase_RTP_custom(P_peaks, R_peaks, T_peaks, length_sig, fs)

% -------------------------------------------------------------------------
% Custom ECG phase mapping:
%
%   R peak = 0
%   T peak = pi/2
%   P peak = -pi/2
%
% Since P occurs before the next R, in an increasing phase representation:
%
%   R(k)     = 2*pi*n
%   T(k)     = pi/2 + 2*pi*n
%   P(k+1)   = 3*pi/2 + 2*pi*n
%   R(k+1)   = 2*pi*(n+1)
%
% Therefore the chronological phase sequence is:
%
%   R -> T -> P -> R
%
% This produces a stable cardiac phase for ECG cycle analysis.
% -------------------------------------------------------------------------

Phase = zeros(1, length_sig);
Omega = zeros(1, length_sig);

% ----------------------------
% 1) Combine and sort peaks
% ----------------------------
all_peaks = [P_peaks(:); R_peaks(:); T_peaks(:)];

% labels:
%   1 = P
%   2 = R
%   3 = T
labels = [ones(numel(P_peaks),1); ...
          2*ones(numel(R_peaks),1); ...
          3*ones(numel(T_peaks),1)];

[pk, idx] = sort(all_peaks);
lb = labels(idx);

% Remove peaks outside valid range
valid_idx = pk >= 1 & pk <= length_sig;
pk = pk(valid_idx);
lb = lb(valid_idx);

% If no valid peaks exist
if isempty(pk)
    warning('No valid ECG peaks were provided.');
    return;
end

% ----------------------------
% 2) Ensure starting with R
% ----------------------------
% Since phase reference is R = 0,
% remove all peaks before the first R.
first_R_idx = find(lb == 2, 1, 'first');

if isempty(first_R_idx)
    warning('No R peaks found. Cannot construct ECG phase.');
    return;
end

pk = pk(first_R_idx:end);
lb = lb(first_R_idx:end);

% ----------------------------
% 3) Enforce exact chronological pattern:
%
%       R -> T -> P -> R -> T -> P ...
%
% labels:
%       R = 2
%       T = 3
%       P = 1
% ----------------------------

expected_pattern = [2 3 1];  % R, T, P

clean_pk = [];
clean_lb = [];

pattern_idx = 1;

i = 1;
while i <= numel(lb)

    expected_label = expected_pattern(pattern_idx);

    if lb(i) == expected_label
        clean_pk(end+1) = pk(i); %#ok<AGROW>
        clean_lb(end+1) = lb(i); %#ok<AGROW>

        pattern_idx = pattern_idx + 1;

        if pattern_idx > numel(expected_pattern)
            pattern_idx = 1;
        end
    end

    i = i + 1;
end

pk = clean_pk(:);
lb = clean_lb(:);

% Need at least two peaks for interpolation
if numel(pk) < 2
    warning('Not enough valid ordered peaks for phase interpolation.');
    return;
end

% ----------------------------
% 4) Assign custom target phases
%
% Chronological sequence:
%
%   R(k) = 2*pi*n
%   T(k) = pi/2 + 2*pi*n
%   P(k+1) = 11*pi/6 + 2*pi*n
%
% where n is the cardiac cycle index.
% ----------------------------

target_phase = zeros(size(pk));

cycle_idx = 0;

for k = 1:numel(pk)

    if lb(k) == 2
        % R peak
        target_phase(k) = 2*pi*cycle_idx;

    elseif lb(k) == 3
        % T peak
        target_phase(k) = pi/2 + 2*pi*cycle_idx;

    elseif lb(k) == 1
        % P peak
        % P is before the next R, equivalent to -pi/2
        % relative to next R, or 3*pi/2 in current cycle.
        target_phase(k) = 3*pi/2 + 2*pi*cycle_idx;

        % After P, the next R belongs to the next cycle
        cycle_idx = cycle_idx + 1;
    end
end

% ----------------------------
% 5) Linear interpolation between consecutive peaks
% ----------------------------

for k = 1:numel(pk)-1

    i1 = pk(k);
    i2 = pk(k+1);

    th1 = target_phase(k);
    th2 = target_phase(k+1);

    bins = i2 - i1;

    if bins <= 0
        continue;
    end

    dtheta = (th2 - th1) / bins;
    omega_val = fs * dtheta;

    theta = th1;

    Phase(i1) = wrapToPi(theta);
    Omega(i1) = omega_val;

    for j = i1+1:i2
        theta = theta + dtheta;
        Phase(j) = wrapToPi(theta);
        Omega(j) = omega_val;
    end
end

% ----------------------------
% 6) After last peak — continue with average cardiac cycle
% ----------------------------

R_locs = pk(lb == 2);

if numel(R_locs) > 1
    RR = mean(diff(R_locs));
else
    RR = fs;  % fallback: assume 1 second cardiac period
end

dtheta = 2*pi / RR;
omega_val = fs * dtheta;

theta = target_phase(end);

for j = pk(end)+1:length_sig
    theta = theta + dtheta;
    Phase(j) = wrapToPi(theta);
    Omega(j) = omega_val;
end

% ----------------------------
% 7) Before first peak
% ----------------------------

theta = target_phase(1);

for j = pk(1)-1:-1:1
    theta = theta - dtheta;
    Phase(j) = wrapToPi(theta);
    Omega(j) = omega_val;
end

end


function [Phase, Omega] = calculate_ECG_phase_RTP_custom_ver1(P_peaks, R_peaks, T_peaks, length_sig, fs)

% -------------------------------------------------------------------------
% Robust custom ECG phase mapping:
%
%   R peak = 0
%   T peak = pi/2
%   P peak = -pi/2
%
% In increasing unwrapped phase:
%
%   R(k)     = 2*pi*n
%   T(k)     = pi/2 + 2*pi*n        if available
%   P(k+1)   = 3*pi/2 + 2*pi*n      if available
%   R(k+1)   = 2*pi*(n+1)
%
% The R peaks are the main anchors.
% Missing T or P peaks do not break phase construction.
% -------------------------------------------------------------------------

Phase = zeros(1, length_sig);
Omega = zeros(1, length_sig);

% ----------------------------
% 1) Clean input peaks
% ----------------------------

P_peaks = P_peaks(:);
R_peaks = R_peaks(:);
T_peaks = T_peaks(:);

P_peaks = P_peaks(P_peaks >= 1 & P_peaks <= length_sig);
R_peaks = R_peaks(R_peaks >= 1 & R_peaks <= length_sig);
T_peaks = T_peaks(T_peaks >= 1 & T_peaks <= length_sig);

P_peaks = unique(sort(P_peaks));
R_peaks = unique(sort(R_peaks));
T_peaks = unique(sort(T_peaks));

if isempty(R_peaks)
    warning('No R peaks found. Cannot construct ECG phase.');
    return;
end

% ----------------------------
% 2) Build robust phase anchors
% ----------------------------

anchor_pk = [];
anchor_phase = [];

num_R = numel(R_peaks);

for r = 1:num_R-1

    R1 = R_peaks(r);
    R2 = R_peaks(r+1);
    RR = R2 - R1;

    if RR <= 0
        continue;
    end

    cycle_idx = r - 1;

    % Main R anchor
    anchor_pk(end+1) = R1; %#ok<AGROW>
    anchor_phase(end+1) = 2*pi*cycle_idx; %#ok<AGROW>

    % ----------------------------
    % T peak inside current R-R interval
    % Expected phase = pi/2
    % Expected time roughly R1 + RR/4
    % ----------------------------
    T_candidates = T_peaks(T_peaks > R1 & T_peaks < R2);

    if ~isempty(T_candidates)
        expected_T = R1 + RR/4;
        [~, idx_T] = min(abs(T_candidates - expected_T));
        T_sel = T_candidates(idx_T);

        anchor_pk(end+1) = T_sel; %#ok<AGROW>
        anchor_phase(end+1) = pi/2 + 2*pi*cycle_idx; %#ok<AGROW>
    end

    % ----------------------------
    % P peak inside current R-R interval
    % P before next R:
    % P phase = -pi/3 relative to next R
    %         = 5*pi/3 relative to current R
    %
    % Expected time roughly R1 + 5*RR/6
    % ----------------------------
    P_candidates = P_peaks(P_peaks > R1 & P_peaks < R2);

    if ~isempty(P_candidates)
        expected_P = R1 + 5*RR/6;
        [~, idx_P] = min(abs(P_candidates - expected_P));
        P_sel = P_candidates(idx_P);

        anchor_pk(end+1) = P_sel; %#ok<AGROW>
        anchor_phase(end+1) = 5*pi/3 + 2*pi*cycle_idx; %#ok<AGROW>
    end

end

% Add last R anchor
last_cycle_idx = num_R - 1;
anchor_pk(end+1) = R_peaks(end);
anchor_phase(end+1) = 2*pi*last_cycle_idx;

% ----------------------------
% 3) Sort anchors by time
% ----------------------------

[anchor_pk, sort_idx] = sort(anchor_pk);
anchor_phase = anchor_phase(sort_idx);

% Remove duplicate sample qrs_positionss if any
[anchor_pk, unique_idx] = unique(anchor_pk, 'stable');
anchor_phase = anchor_phase(unique_idx);

if numel(anchor_pk) < 2
    warning('Not enough anchors for phase interpolation.');
    return;
end

% ----------------------------
% 4) Linear interpolation between anchors
% ----------------------------

for k = 1:numel(anchor_pk)-1

    i1 = anchor_pk(k);
    i2 = anchor_pk(k+1);

    th1 = anchor_phase(k);
    th2 = anchor_phase(k+1);

    bins = i2 - i1;

    if bins <= 0
        continue;
    end

    dtheta = (th2 - th1) / bins;
    omega_val = fs * dtheta;

    theta = th1;

    Phase(i1) = wrapToPi(theta);
    Omega(i1) = omega_val;

    for j = i1+1:i2
        theta = theta + dtheta;
        Phase(j) = wrapToPi(theta);
        Omega(j) = omega_val;
    end
end

% ----------------------------
% 5) Estimate average RR for before/after edges
% ----------------------------

if numel(R_peaks) > 1
    RR_mean = mean(diff(R_peaks));
else
    RR_mean = fs;
end

dtheta_edge = 2*pi / RR_mean;
omega_edge = fs * dtheta_edge;

% ----------------------------
% 6) After last anchor
% ----------------------------

theta = anchor_phase(end);

for j = anchor_pk(end)+1:length_sig
    theta = theta + dtheta_edge;
    Phase(j) = wrapToPi(theta);
    Omega(j) = omega_edge;
end

% ----------------------------
% 7) Before first anchor
% ----------------------------

theta = anchor_phase(1);

for j = anchor_pk(1)-1:-1:1
    theta = theta - dtheta_edge;
    Phase(j) = wrapToPi(theta);
    Omega(j) = omega_edge;
end

end

function Phase = calculate_dtw_phase(ecg, r_peaks)

N = length(ecg);            % Total signal length
Phase = zeros(1,N);         % Phase array

nBeats = length(r_peaks);

%% -------- Extract beats centered on R-peaks
beats = {};                 % Cell array of beats
beat_index = {};            % qrs_positionss of each beat in original ECG
L = [];                     % Beat lengths

for i = 2:nBeats-1

    RR_prev = r_peaks(i) - r_peaks(i-1);   % Previous RR interval
    RR_next = r_peaks(i+1) - r_peaks(i);   % Next RR interval

    % Beat boundaries centered on R-peak
    start_idx = round(r_peaks(i) - RR_prev/2);
    end_idx   = round(r_peaks(i) + RR_next/2);

    start_idx = max(1,start_idx);
    end_idx   = min(N,end_idx);

    beat = ecg(start_idx:end_idx);

    beats{end+1} = beat;
    beat_index{end+1} = start_idx:end_idx;
    L(end+1) = length(beat);
end

%% -------- Build average template beat
% template = zeros(1,Lref);
template = beats{3};
Lref = length(template);          % Reference length

% for i=1:length(beats)
%     b = resample(beats{i},Lref,length(beats{i}));
%     template = template + b;
% end
% 
% template = template / length(beats);

%% -------- Define template phase (R at center)
phase_template = linspace(-pi,pi,Lref);

%% -------- DTW alignment and phase mapping
for i=1:length(beats)

    beat = beats{i};

    % DTW alignment qrs_positionss
    [~,ix,iy] = dtw(beat,template);

    phi = zeros(1,length(beat));

    % Assign phase using warping path
    for k=1:length(ix)
        phi(ix(k)) = phase_template(iy(k));
    end

    % Interpolate missing phase samples
    phi = interp1(1:length(phi),phi,1:length(phi),'linear','extrap');

    Phase(beat_index{i}) = phi;
end

%% -------- Phase extrapolation before first beat
first_idx = beat_index{1}(1);
RR_mean = mean(diff(r_peaks));
stepTheta = 2*pi/RR_mean;

theta = Phase(first_idx);

for j = first_idx-1:-1:1
    theta = theta - stepTheta;
    theta = mod(theta+pi,2*pi)-pi;
    Phase(j) = theta;
end

%% -------- Phase extrapolation after last beat
last_idx = beat_index{end}(end);
theta = Phase(last_idx);

for j = last_idx+1:N
    theta = theta + stepTheta;
    theta = mod(theta+pi,2*pi)-pi;
    Phase(j) = theta;
end
end







function [DTW_matrix,d,optimal_path]=dtw_ver2(sig1,sig2,w)
% Copyright (C) 2013 Quan Wang <wangq10@rpi.edu>,
% Signal Analysis and Machine Perception Laboratory,
% Department of Electrical, Computer, and Systems Engineering,
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

% dynamic time warping of two signals
% s: signal 1, size is ns*k, row for time, colume for channel 
% t: signal 2, size is nt*k, row for time, colume for channel 
% w: window parameter
%      if s(i) is matched with t(j) then |i-j|<=w
% d: resulting distance

if nargin<3
    w=Inf;
end

N=length(sig1);
M=length(sig2);
if size(sig1,2)~=size(sig2,2)
    error('Error in dtw(): the dimensions of the two input signals do not match.');
end
w=max(w, abs(N-M)); % adapt window size

%% initialization
DTW_matrix=zeros(N+1,M+1)+Inf; % cache matrix
DTW_matrix(1,1)=0;

%% begin dynamic programming
for i=1:N
    for j=max(i-w,1):min(i+w,M)
        cost=norm(sig1(i,:)-sig2(j,:));
        DTW_matrix(i+1,j+1)=cost+min( [DTW_matrix(i,j+1), DTW_matrix(i+1,j), DTW_matrix(i,j)] );
        
    end
end
d=DTW_matrix(N+1,M+1);

n=N;
m=M;
k=1;
optimal_path=[];
optimal_path(1,:)=[N,M];
while ((n+m)~=2)
    if (n-1)==0
        m=m-1;
    elseif (m-1)==0
        n=n-1;
    else 
      [values,number]=min([1*DTW_matrix(n-1,m),1*DTW_matrix(n,m-1),1*DTW_matrix(n-1,m-1)]);
      switch number
      case 1
        n=n-1;
      case 2
        m=m-1;
      case 3
        n=n-1;
        m=m-1;
      end
  end
    k=k+1;
    optimal_path=cat(1,optimal_path,[n,m]);
end
end

function psi = psi_triangular(t)

phi = @(x) max(1-abs(x),0);   % triangular function

psi = -0.5*phi(2*t) + phi(2*t-1) - 0.5*phi(2*t-2);

end