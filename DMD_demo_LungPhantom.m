addpath(genpath('Utils')); % Add necessary files to directory

%% Simulate phantom for m=480 and add noise
load('LungPhantom.mat'); % Load phantom (a single snapshot)

rng(160123); % random number generator seed for repeatability

% Parameters for modified Lujan formula from Bauman & Bieri (10.1002/mrm.26096). 
A_DC = 80;
A_R = 16;
A_C = 4;
phase_r = pi/4;
phase_c = pi/2;

dt = 0.25; 
T = 120; % Acquisition for 120 seconds
m = T/dt; % Number of acquisitions (with default parameters m=480)

t = dt:dt:T;
Fs = 1/dt; % Sampling rate

f_Resp_start = 0.2; % Starting respiration frequency (in Hz)
f_Resp_end = 0.18; % Ending respiration frequency (in Hz)

f_Card_start = 1.0; % Starting cardiac frequency (in Hz)
f_Card_end = 0.9; % Ending cardiac frequency (in Hz)

tt = 5/dt; % Transition time = 5 seconds
tstr = (m-tt)/2; % Transition start
tend = tstr + tt; % Transition end
% Angular frequencies for respiration and pulsation
omega_r = 2*pi*[f_Resp_start*ones(1,tstr) linspace(f_Resp_start,f_Resp_end,tt) f_Resp_end*ones(1,tend)];
omega_c = 2*pi*[f_Card_start*ones(1,tstr) linspace(f_Card_start,f_Card_end,tt) f_Card_end*ones(1,tend)];

backDens = 1; % Assume 1 signal density for background tissue
lungDens = 0.2; % Assume 0.2 signal density for lung 
vessDens = 0.9; % Assume 0.9 signal density for large blood vessels

% Create a combined segmentation map
segMap = backGround + 2*lungParencCor + 3*largeVesselsCor;

tissueIndx = segMap==2; % Indices of lung tissue

% Simulate noiseless phantom for m=480
noiselessPhantom = zeros(sx,sy,m);
for xIdx = 1:sx
    for yIdx = 1:sy
        if segMap(xIdx,yIdx)==1 % Background. Only DC
            for nIdx = 1:m
                noiselessPhantom(xIdx,yIdx,nIdx) = backDens*A_DC;
            end
        elseif segMap(xIdx,yIdx)==2 % Lung Parenchyma. DC + Vent + Perf
            for nIdx = 1:m
                noiselessPhantom(xIdx,yIdx,nIdx) = lungDens*(A_DC - A_R * cos(nIdx*(omega_r(nIdx)*dt/2)+phase_r ).^2 +  A_C * cos(nIdx*(omega_c(nIdx)*dt/2)+phase_c).^2);
            end
        elseif segMap(xIdx,yIdx)==3 % Large Vessels. DC + Perf
            for nIdx = 1:m
                noiselessPhantom(xIdx,yIdx,nIdx) = vessDens*(A_DC +  A_C * cos(nIdx*(omega_c(nIdx)*dt/2)+phase_c).^2);
            end
        end
    end
end

% Visualize noiseless phantom
% implay(noiselessPhantom/80)

% Create noise patterns and noisy phantom
snr = 50;
noiseLevel = max(abs(noiselessPhantom(95,95,:)))/snr./sqrt(2); % Pick a large vessel for SNR level

% Generate complex Gaussian noise
noisePattern = noiseLevel*randn([sx,sy,m])+sqrt(-1)*noiseLevel*randn([sx,sy,m]); 
% Add its magnitude to the phantom
phantom = noiselessPhantom + abs(noisePattern);

% Visualize phantom
% implay(phantom/80)

%% Dynamic Mode Decomposition (DMD)
X = reshape(phantom,[sx*sy,m]); % Vectorize the snapshots

stackNum = 5; % Data stacking number
DMDvarin = {'dt',dt,'nstacks',stackNum};
[Phi_DMD, omega_DMD, lambda_DMD, b_DMD, freq_DMD, Xdmd_DMD , rDMD] = DynamicModeDecomp(X, DMDvarin);

res_DMD = reshape(Phi_DMD(1:sx*sy,:),[sx,sy,rDMD]); % extract the current-state DMD modes from Phi_DMD by pulling out the first sx*sy rows.

%% Find ventilation and perfusion related signal changes via their frequencies
ventRange = [0.05 0.35]; % Look for ventilation around these frequencies
perfRange = [0.75 1.25]; % Look for perfusion around these frequencies
vent_DMD_idx = find(freq_DMD>ventRange(1) & freq_DMD<ventRange(2)); % Ventilation indices
perf_DMD_idx = find(freq_DMD>perfRange(1) & freq_DMD<perfRange(2) & abs(lambda_DMD)>0.80); % Perfusion indices
% Here with "abs(lambda_DMD)>0.80" decaying modes can be quickly omitted
% via magnitude of eigenvalues.

idxDC_DMD = find(abs(freq_DMD)<5e-4); % Find DC frequencies // DC indices
dc_DMD = reconstructFreqImage(b_DMD/2,res_DMD,idxDC_DMD); % zero-freq image (Division by 2 due to missing symmetry in zero frequency)
vent_DMD = reconstructFreqImage(b_DMD,res_DMD,vent_DMD_idx); % ventilation image
perf_DMD = reconstructFreqImage(b_DMD,res_DMD,perf_DMD_idx); % perfusion image

%% Generate functional maps
% Calculate fractional ventilation map
BGr = std(dc_DMD(1:30,1:30)); % A region from background
BG = std(BGr(:)); % noise level in background 
ventMap = abs((vent_DMD)./((vent_DMD/2)+dc_DMD-BG));

% Calculate normalized perfusion map (w.r.t large vessels)
perfp = prctile(perf_DMD(:),99); % Here I assume large vessels contribute to the highest signal levels in the area of interest. For in vivo data, looking at an ROI might be needed. 
perf_DMD(perf_DMD>perfp) = perfp; % Reduce outliers 
perfMap = perf_DMD/perfp;  % Normalize

%% Display results
ventColorMap = ocean; % Ventilation colormap (Perceptually uniform colormaps)
perfColorMap = blackbody; % Perfusion colormap (Perceptually uniform colormaps)

figure;
ax1 = subplot(2,2,1); imshow(phantom(:,:,1),[]); colorbar, title('Phantom [a.u.]')
ax2 = subplot(2,2,2); imshow(dc_DMD,[]); colorbar,         title('DC Component [a.u.]')
ax3 = subplot(2,2,3); imshow(ventMap,[0 0.2]); colorbar,   title('Fractional Ventilation [ml/ml]'), 
ax4 = subplot(2,2,4); imshow(perfMap,[0 0.4]); colorbar,   title('Perfusion [normalized]'), 
 
colormap(ax1,gray)
colormap(ax2,gray)
colormap(ax3,ventColorMap)
colormap(ax4,perfColorMap)