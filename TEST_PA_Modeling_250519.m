%% PA Modeling and DPD Implementation 
% The example show to implement the PA model and apply DPD methoology for correction.
% History
%% 
% * 2025-05-19, Draft
% Set Test Waveform

clear paInSig paOutSig paOutSig_ideal
load('WIFI_EHT320_MCS13_15.80dBm_NumPkts4_1280MHz_0527SIMCW_0dB.mat') % load waveform
wvParams;
paInSig = wvParams.PaInSignal;
try
    paOutSig = wvParams.PaOutSignal;
catch
    paOutSig = wvParams.PaOutSignal_ideal;
end
try
    paOutSig_ideal = wvParams.PaOutSignal_ideal;
catch
    paOutSig_ideal = paOutSig;
end
try
    PaModelFlatnessPerCBW = wvParams.PaModelFlatnessPerCBW;
catch
    PaModelFlatnessPerCBW = 0;
end
try
    isPaModelFlatCompensateLoss = wvParams.PaModelFlatCompensateLoss;
catch
    isPaModelFlatCompensateLoss = 0;
end
if PaModelFlatnessPerCBW~=0
    isPaFlatness = 1;
    b_coefs_flat = wvParams.PaModelFlatFirCoefs;
else
    isPaFlatness = 0;
    b_coefs_flat = 1;
end
aclrVal_wv = wvParams.PlotACLR; % plot
pltComm = wvParams.PlotComm; % plot
dlSlots = wvParams.usefulWvfm;
fsMHz = wvParams.SampleRateMHz;
signalFormat = wvParams.SignalFormat;
MCS_index = wvParams.MCS;

% Synchronize Waveform between Input and Output

delay_set = 5
x = paInSig;
y = circshift(paOutSig, delay_set);
if 1
    [corr_val, lag] = xcorr(x, y);
    [~, idx_max_corr] = max(abs(corr_val));
    delay_get = -lag(idx_max_corr);
else
    corr_val = ifft( conj(fft(x)) .* fft(y) );
    [~, idx_max_corr] = max(abs(corr_val))
    delay_get = idx_max_corr-1;
end
y_shift = circshift(y, -delay_get); % correction
% Plot Waveform

nSamps = numel(paOutSig);
abs_pa_output = abs(fftshift(fft(paOutSig)))/nSamps;
figure(051901), plot(20*log10(abs_pa_output*sqrt(1000)))
% Plot AMAM and AMPM

min_power_dBm = -70;
max_power_dBm = Inf;
% sample sort and select
pwr_pa_input_dBm = 10*log10(abs(paInSig).^2/1 * 1000);
idx_min = (pwr_pa_input_dBm>min_power_dBm);
idx_max = (pwr_pa_input_dBm<=max_power_dBm);
pa_input_sig = paInSig(idx_min & idx_max);
pa_output_sig = paOutSig(idx_min & idx_max);
% comput power
pwr_pa_input_dBm = 10*log10(abs(pa_input_sig).^2/1 * 1000);
pwr_pa_output_dBm = 10*log10(abs(pa_output_sig).^2/1 * 1000);

% plot amam
figure(051902), plot(pwr_pa_input_dBm, pwr_pa_output_dBm)

% plot ampm
angle_diff = angle(pa_output_sig ./ pa_input_sig);
phs_diff_deg = (angle_diff) / pi * 180;
figure(051903), scatter(pwr_pa_input_dBm, phs_diff_deg), hold on

% Plot AMAM and AMPM

polyOrder = 7
memoryDepth = 2
coefs = PAMemoryPoly.get_coefs_memory_poly(pa_input_sig, pa_output_sig, polyOrder, memoryDepth)
paOutSigFit = PAMemoryPoly.fit_pa_memory_poly(paInSig, coefs, (1:1:polyOrder), (0:1:memoryDepth));
% export
pa_output_sig_fit = paOutSigFit(1+memoryDepth:end,:);
pa_output_sig = paOutSig(1+memoryDepth:end,:);
pa_input_sig = paInSig(1+memoryDepth:end,:);
pwr_pa_input_dBm = 10*log10(abs(pa_input_sig).^2/1 * 1000);
pwr_pa_output_fit_dBm = 10*log10(abs(pa_output_sig_fit).^2/1 * 1000);
% plot amam
figure(051902), scatter(pwr_pa_input_dBm, pwr_pa_output_fit_dBm), xlim([-70, Inf])
% plot ampm
angle_diff_fit = angle(pa_output_sig_fit ./ pa_input_sig);
phs_diff_deg_fit = (angle_diff_fit) / pi * 180;
figure(051903), scatter(pwr_pa_input_dBm, phs_diff_deg_fit), xlim([-70, Inf])
% NMSE
x = pa_output_sig;
y = pa_output_sig_fit;
NMSE = sum(abs(x-y).^2) / sum(abs(x).^2);
NMSE_dB = 10*log10(NMSE);
EVM = sqrt(mean(abs(x-y).^2) / mean(abs(x).^2));
if 0 % debug
end