%DEMORATIOCALIBRATION  Demo: ratiometric measurements and roGFP redox calibration.
%
% Demonstrates ratioCompute and redoxCalibrate on synthetic two-channel
% fluorescence data representing grx1-roGFP2 imaging:
%   Section 1 – Basic ratio computation with spatial gradient
%   Section 2 – Autofluorescence correction effect
%   Section 3 – Henderson-Hasselbalch calibration for pH probes
%   Section 4 – redoxCalibrate: OxD normalisation curve vs iFactor
%   Section 5 – redoxCalibrate: Nernst calibration curve (OxD → Eh)
%   Section 6 – Apply full roGFP calibration to synthetic image, RGB output
%
% Run with:
%   demoRatioCalibration
%
% Path requirements:
%   Physiology_sandbox/src/ and Physiology_sandbox/utils/ must be on path.

close all; clc;
fprintf('=== demoRatioCalibration ===\n\n');

% =========================================================================
% Section 1 — Basic ratio with spatial gradient
% =========================================================================
fprintf('Section 1: Basic ratio computation...\n');

nY = 64; nX = 64; nT = 6;

% Simulate grx1-roGFP2: 405nm (oxidised, Ch1) and 488nm (reduced, Ch2)
% Create a spatial gradient: left = more reduced, right = more oxidised
[X, ~] = meshgrid(linspace(0,1,nX), linspace(0,1,nY));
oxdTrue = 0.1 + 0.8 * X;   % true OxD: 0.1 (left) to 0.9 (right)

% roGFP2 parameters
rMin = 0.1; rMax = 1.0; iFactor = 5.5;
% Invert OxD formula to get R: R = [iFactor*rMax*OxD + rMin*(1-OxD)] / [iFactor*OxD + (1-OxD)]
oxdTrue_safe = max(0.001, min(0.999, oxdTrue));
Rtrue = (iFactor * rMax .* oxdTrue_safe + rMin .* (1 - oxdTrue_safe)) ./ ...
        (iFactor .* oxdTrue_safe + (1 - oxdTrue_safe));

% Build two-channel image: num (405nm signal proportional to OxD),
%   den (488nm signal inversely proportional to OxD)
F405 = single(0.5 + 0.3 * oxdTrue);       % [0.5, 0.8] range
F488 = single(0.6 - 0.2 * oxdTrue);       % [0.6, 0.4] range
imIn = cat(3, repmat(F405,[1 1 1 1 nT]), repmat(F488,[1 1 1 1 nT]));
imIn = reshape(imIn, [nY nX 2 1 nT]);
imIn = imIn + 0.01 * randn(size(imIn), 'single');

p_ratio = struct('numeratorChannel',1, 'denominatorChannel',2);
[ratioSeries, ratioMean, ~, ~] = ratioCompute(imIn, p_ratio);

figure('Name','Section 1: Basic Ratio','Position',[50 500 750 260]);
subplot(1,3,1); imagesc(squeeze(imIn(:,:,1,1,1))); axis image off; colorbar;
title('Ch1 (405nm, numerator)'); colormap(gca,'hot');
subplot(1,3,2); imagesc(squeeze(imIn(:,:,2,1,1))); axis image off; colorbar;
title('Ch2 (488nm, denominator)'); colormap(gca,'hot');
subplot(1,3,3); imagesc(ratioMean(:,:,1,1)); axis image off; colorbar;
title('Ratio Ch1/Ch2'); colormap(gca,'parula');
fprintf('  → Ratio range: %.2f–%.2f (expected %.2f–%.2f)\n', ...
    min(ratioMean(:)), max(ratioMean(:)), min(Rtrue(:)), max(Rtrue(:)));

% =========================================================================
% Section 2 — Autofluorescence correction
% =========================================================================
fprintf('\nSection 2: Autofluorescence correction...\n');

% Add an AF channel (e.g., 488nm excitation with a different emission filter)
AFsig = single(0.1 * ones(nY, nX));   % uniform AF at 0.1
imIn3 = cat(3, imIn, reshape(repmat(AFsig,[1 1 1 1 nT]),[nY nX 1 1 nT]));

p_noAF = struct('numeratorChannel',1,'denominatorChannel',2);
p_AF   = struct('numeratorChannel',1,'denominatorChannel',2, ...
                'autofluorescence',struct('use',true,'channel',3, ...
                'scaleNumerator',0,'scaleDenominator',1.0));

[~, rMeanNoAF, ~, ~] = ratioCompute(imIn3, p_noAF);
[~, rMeanAF,   ~, ~] = ratioCompute(imIn3, p_AF);

figure('Name','Section 2: AF Correction','Position',[50 180 700 260]);
subplot(1,3,1); imagesc(rMeanNoAF(:,:,1,1)); axis image off; colorbar;
title('Ratio: no AF correction'); colormap(gca,'parula');
subplot(1,3,2); imagesc(rMeanAF(:,:,1,1)); axis image off; colorbar;
title('Ratio: AF subtracted from Ch2'); colormap(gca,'parula');
subplot(1,3,3); imagesc(rMeanAF(:,:,1,1) - rMeanNoAF(:,:,1,1));
axis image off; colorbar; title('\Delta ratio (AF correction)'); colormap(gca,'rdbu');
fprintf('  → Mean ratio change from AF correction: %.3f\n', ...
    mean(rMeanAF(:) - rMeanNoAF(:)));

% =========================================================================
% Section 3 — Henderson-Hasselbalch pH calibration
% =========================================================================
fprintf('\nSection 3: Henderson-Hasselbalch pH calibration...\n');

% Simulate a pH probe: ratio varies from Racid to Rbase across image
Racid = 0.2; Rbase = 0.9; pKa = 6.8;
R_sweep = single(linspace(Racid + 0.01, Rbase - 0.01, nX));
imPH = zeros(nY, nX, 2, 1, 1, 'single');
imPH(:,:,1,1,1) = repmat(R_sweep, [nY 1]);   % numerator = ratio directly
imPH(:,:,2,1,1) = ones(nY, nX, 'single');    % denominator = 1

p_ph = struct('numeratorChannel',1,'denominatorChannel',2, ...
              'calibration', struct('method','henderson-hasselbalch', ...
                                    'pKa',pKa,'Racid',Racid,'Rbase',Rbase));
[~, ~, ~, calibMean] = ratioCompute(imPH, p_ph);

Raxis = double(R_sweep);
pHfromHH = pKa + log10((Raxis - Racid) ./ (Rbase - Raxis));

figure('Name','Section 3: H-H pH Calibration','Position',[800 500 600 280]);
subplot(1,2,1);
plot(Raxis, pHfromHH, 'b-', 'LineWidth', 2);
xlabel('Ratio R'); ylabel('pH'); grid on;
xline(Racid,'k--','R_{acid}'); xline(Rbase,'k--','R_{base}');
yline(pKa, 'r--', sprintf('pKa = %.1f', pKa));
title('Henderson-Hasselbalch calibration curve');

subplot(1,2,2);
imagesc(calibMean(:,:,1,1)); axis image off; colorbar;
title('Calibrated pH map');
colormap(gca, jet(128));
fprintf('  → pH range in image: %.2f – %.2f (pKa = %.1f)\n', ...
    min(calibMean(isfinite(calibMean(:)))), ...
    max(calibMean(isfinite(calibMean(:)))), pKa);

% =========================================================================
% Section 4 — OxD normalisation curves vs iFactor
% =========================================================================
fprintf('\nSection 4: OxD normalisation curves for different iFactor values...\n');

R_vals  = linspace(rMin, rMax, 200);
iFactors = [1.0, 2.0, 5.5, 10.0];
colors   = lines(numel(iFactors));

figure('Name','Section 4: OxD vs Ratio','Position',[50 500 520 380]);
hold on;
for k = 1:numel(iFactors)
    delta   = iFactors(k);
    OxD_k   = (R_vals - rMin) ./ (delta .* (rMax - R_vals) + (R_vals - rMin));
    plot(R_vals, OxD_k, '-', 'Color', colors(k,:), 'LineWidth', 2);
end
xlabel('Measured ratio R');  ylabel('Oxidised fraction (OxD)');
title('OxD normalisation: sensitivity to instrument factor \delta');
legend(arrayfun(@(d) sprintf('\\delta = %.1f',d), iFactors, 'UniformOutput',false), ...
    'Location','nw');
xline(rMin,'k:'); xline(rMax,'k:');
yline(0.5, 'k--', 'OxD = 0.5 (E_{hh} = E_0'')');
grid on; box on; ylim([0 1]);
fprintf('  → Lower iFactor compresses OxD range; higher iFactor stretches low-OxD region.\n');

% =========================================================================
% Section 5 — Nernst calibration curve (OxD → Eh)
% =========================================================================
fprintf('\nSection 5: Nernst calibration curves at different temperatures...\n');

OxD_vals = linspace(0.01, 0.99, 200);
temps     = [277, 293, 310];   % 4°C, 20°C, 37°C
midpoint  = -280;              % E°' for roGFP2 (mV)
n = 2; R_gas = 8.314; F_fara = 96485;
colors2 = [0 0.4 1; 0 0.7 0; 1 0.2 0];

figure('Name','Section 5: Nernst Curves','Position',[600 500 560 380]);
subplot(1,2,1);
hold on;
for k = 1:numel(temps)
    Eh_k = midpoint - (1000*R_gas*temps(k)/(n*F_fara)) .* log((1-OxD_vals)./OxD_vals);
    plot(OxD_vals, Eh_k, '-', 'Color', colors2(k,:), 'LineWidth', 2);
end
xlabel('OxD'); ylabel('E_{hh} (mV vs SHE)');
title('Nernst equation: temperature sensitivity');
legend(arrayfun(@(T) sprintf('T = %d K (%d°C)',T,T-273), temps,'UniformOutput',false), ...
    'Location','ne');
xline(0.5,'k--','OxD = 0.5'); yline(midpoint,'k:','E_0'''); grid on; box on;

subplot(1,2,2);
hold on;
pHvals    = [6.5, 7.0, 7.4, 8.0];
T_ref = 293;
for k = 1:numel(pHvals)
    E0prime_k = midpoint - (1000*(2.303*R_gas*T_ref)/F_fara) * (pHvals(k) - 7.0);
    Eh_k = E0prime_k - (1000*R_gas*T_ref/(n*F_fara)) .* log((1-OxD_vals)./OxD_vals);
    plot(OxD_vals, Eh_k, '-', 'LineWidth', 2);
end
xlabel('OxD'); ylabel('E_{hh} (mV vs SHE)');
title('Nernst equation: pH dependence');
legend(arrayfun(@(pH) sprintf('pH = %.1f',pH), pHvals,'UniformOutput',false), ...
    'Location','ne');
grid on; box on;
fprintf('  → pH shift of 59 mV per pH unit visible as vertical offset.\n');

% =========================================================================
% Section 6 — Full roGFP calibration: ratio → OxD → Eh with RGB output
% =========================================================================
fprintf('\nSection 6: Full roGFP calibration pipeline with RGB output...\n');

p_calib = struct('rMin',rMin, 'rMax',rMax, 'iFactor',iFactor, ...
                 'midpoint',-280, 'temp',293, 'pH',7.0, 'pH0',7.0, 'n',2);

[oxdSeries, ehSeries] = redoxCalibrate(ratioMean, p_calib);

% ER mask from blob structure
erMask = logical(repmat((blobA + blobB) > 0.05, [1 1 1 1 1]));

p_rgb = struct('ratioMin',rMin,'ratioMax',rMax, ...
               'ehMin',-350,'ehMax',-200,'whiteUse',true);

[rRgb, oxdRgb, ehRgb] = ratioRgb(imIn(:,:,:,:,1), erMask, ratioMean, ...
                                   oxdSeries, ehSeries, p_rgb);

figure('Name','Section 6: roGFP RGB Outputs','Position',[50 50 1100 320]);
subplot(1,3,1); imshow(squeeze(rRgb(:,:,1,1,:)));
title('Ratio (R2 colormap)');
subplot(1,3,2); imshow(squeeze(oxdRgb(:,:,1,1,:)));
title('OxD 0→1 (R2 colormap)');
subplot(1,3,3); imshow(squeeze(ehRgb(:,:,1,1,:)));
title('E_{hh} mV (D7 diverging)');

fprintf('  → OxD range: %.3f – %.3f\n', min(oxdSeries(:)), max(oxdSeries(:)));
fprintf('  → Eh range:  %.1f – %.1f mV\n', min(ehSeries(:)), max(ehSeries(:)));
fprintf('\n=== demoRatioCalibration complete ===\n');
