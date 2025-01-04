%% Task:
% This exercise starts with 1D GRE sequence (ex04_fid2gre1d). Could you try
% to extend the 1D GRE to 2D GRE with matrix size of 256*256 by
% implementing a set of phase-encoding (PE) gradients?

clear; close all; clc;
seq_name = "ex16_gre_demo_2d_label";
userID = 'zheng';

%% Set system and sequnece parameters
sys = mr.opts('MaxGrad', 22, ...
    'GradUnit', 'mT/m', ...
    'MaxSlew', 120, ...
    'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, ...
    'adcDeadTime', 20e-6 ...
);

seq = mr.Sequence(sys);
fov = 256e-3; Nx = 256; Ny = 1;
alpha = 10;
Nslices = 1;
sliceThickness = 3e-3;
sliceGap = 0;
TR = 20e-3;
TE = 6e-3;
roDuration = 5.12e-3;

%% Create events
% Create alpha-degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(alpha * pi / 180, 'Duration', 3e-3, ...
    'SliceThickness', sliceThickness, 'apodization', 0.42, 'timeBwProduct', 4, 'system', sys);

% Define other gradients and ADC events
deltak = 1 / fov;
gx = mr.makeTrapezoid('x', 'FlatArea', Nx * deltak, 'FlatTime', roDuration, 'system', sys);
adc = mr.makeAdc(Nx, 'Duration', gx.flatTime, 'Delay', gx.riseTime, 'system', sys);
gxPre = mr.makeTrapezoid('x', 'Area', -gx.amplitude * (adc.dwell * (adc.numSamples / 2 + 0.5) + 0.5 * gx.riseTime), ...
    'Duration', 1e-3, 'system', sys);
gzReph = mr.makeTrapezoid('z', 'Area', -gz.area / 2, 'Duration', 1e-3, 'system', sys);
phaseAreas = ((0:Ny - 1) - Ny / 2) * deltak;
gy = mr.makeTrapezoid('y', 'Area', max(abs(phaseAreas)), 'Duration', mr.calcDuration(gxPre), 'system', sys);
peScales = phaseAreas / gy.area;

for iY = 1:Ny
    gyPre(iY) = mr.scaleGrad(gy, peScales(iY));
end

%% Calculate timing
delayTE = ceil((TE - mr.calcDuration(gxPre) - gz.fallTime - gz.flatTime / 2 ...
    - mr.calcDuration(gx) / 2) / seq.gradRasterTime) * seq.gradRasterTime;
delayTR = ceil((TR - mr.calcDuration(gz) - mr.calcDuration(gxPre) ...
    - mr.calcDuration(gx) - delayTE) / seq.gradRasterTime) * seq.gradRasterTime;
assert(delayTE >= 0);
assert(delayTR >= 0);

%% Loop over phase encodes and define sequence blocks
tic;

for i = 1:Ny
    seq.addBlock(rf, gz); % slice-selective excitation
    seq.addBlock(gxPre, gzReph); % hint: add gyPre(i) here
    seq.addBlock(mr.makeDelay(delayTE));
    seq.addBlock(gx, adc);
    seq.addBlock(mr.makeDelay(delayTR));
end

toc;

%% check whether the timing of the sequence is correct
[ok, error_report] = seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% prepare sequence export
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'gre');

seq.write(['ex11_gre1d2gre2d_' userID '.seq']); % Write to pulseq file

%% plot sequence and k-space diagrams

seq.plot('timeRange', [0 5] * TR);

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(ktraj(1, :), ktraj(2, :), 'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold; plot(ktraj_adc(1, :), ktraj_adc(2, :), 'r.'); % plot the sampling points
