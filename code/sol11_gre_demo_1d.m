%% Task: 
% Make a 1D GRE sequence

clear; close all; clc;
root = "E:\UIHBranch\BRANCHES\uMRV10_3T_R3\UIH\appdata\MRSiteData\Share\Imaging\pulseq\";
seq_name = "ex11_gre_demo_1d";
userID = '_lz';
seqName = strcat(seq_name, userID, ".seq");

%% 1. Define system and sequnece parameters
sys = mr.opts( ...
    'MaxGrad', mr.convert(40,'mT/m'), ...
    'MaxSlew', mr.convert(109,'mT/m/ms'), ...
    'gradRasterTime', 10e-6,...
    'rfRingdownTime', 10e-6, ...
    'rfDeadTime', 400e-6, ...
    'adcDeadTime', 70e-6, ...
    'adcRasterTime', 200e-9, ...
    'blockDurationRaster', 10e-6,...
    'B0', 3 ...
    );

seq = mr.Sequence(sys);
fov = 256e-3; Nx = 64;
alpha = 90;
Nslices = 1;
sliceThickness = 3e-3;
TR = 20e-3;
TE = 6e-3;
roDuration = 5.12e-3;

%% 2. Define events
% Create alpha-degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(alpha * pi / 180, 'Duration', 3e-3, ...
    'SliceThickness', sliceThickness, 'apodization', 0.42, 'timeBwProduct', 4, 'system', sys);
gzRephase = mr.makeTrapezoid('z', 'Area', -gz.area/2, 'Duration', 1e-3, 'system', sys);

% Define readout gradients and ADC events
deltak = 1 / fov;
gx = mr.makeTrapezoid('x', 'FlatArea', Nx * deltak, 'FlatTime', roDuration, 'system', sys);
adc = mr.makeAdc(Nx, 'Duration', gx.flatTime, 'Delay', gx.riseTime, 'system', sys);
gxPre = mr.makeTrapezoid('x', 'Area', -gx.area/2, 'Duration', 1e-3, 'system', sys);

%% 3. Calculate timing
delayTE = ceil((TE - mr.calcDuration(gxPre) - gz.fallTime - gz.flatTime / 2 ...
    - mr.calcDuration(gx) / 2) / seq.gradRasterTime) * seq.gradRasterTime;
delayTR = ceil((TR - mr.calcDuration(gz) - mr.calcDuration(gx) - delayTE) / seq.gradRasterTime) * seq.gradRasterTime;
assert(delayTE >= 0);
assert(delayTR >= 0);

%% 4. Construct sequence blocks

seq.addBlock(rf, gz); 
seq.addBlock(gzRephase, gxPre); 
seq.addBlock(mr.makeDelay(delayTE));
seq.addBlock(gx, adc);
seq.addBlock(mr.makeDelay(delayTR));

%% 5. Check whether the timing of the sequence is correct
[ok, error_report] = seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% 6. Prepare sequence export
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'gre');

seq.write(strcat("seq\", seqName));
seq.write(strcat(root, seqName));

%% 7. Plot sequence and k-space diagrams
seq.plot('timeRange', [0 5] * TR);

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-space trajectories
figure; plot(ktraj(1, :), ktraj(2, :), 'b');
axis('equal'); 
hold; plot(ktraj_adc(1, :), ktraj_adc(2, :), 'r.');
