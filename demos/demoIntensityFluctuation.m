%DEMOINTENSITYFLUCTUATION  Demo: temporal intensity statistics for calcium imaging.
%
% Demonstrates intensityFluctuationCompute on a synthetic calcium transient:
%   Section 1 – Generate synthetic single-channel image with a step transient
%   Section 2 – Compare variance, CV, and dF/F methods side-by-side
%   Section 3 – Effect of baseline window choice on dF/F
%   Section 4 – Visualise dF/F map with physiologyConvertToRgb
%
% Run with:
%   demoIntensityFluctuation
%
% Path requirements:
%   Physiology_sandbox/src/ and Physiology_sandbox/utils/ must be on path.

close all; clc;
fprintf('=== demoIntensityFluctuation ===\n\n');

% -------------------------------------------------------------------------
% Section 1 — Synthetic calcium transient
% -------------------------------------------------------------------------
% Image: 64×64 px, 1 channel, 1 Z, 20 time frames.
% Spatial structure: two Gaussian "ER" blobs with different response amplitudes.
% Temporal structure: flat baseline (frames 1-5), step increase at frame 6
%   (+80% in blob A, +40% in blob B), gradual recovery with time constant 4 frames.

nY = 64; nX = 64; nT = 20;
fprintf('Section 1: Generating synthetic calcium transient (%d frames)...\n', nT);

% Spatial basis: two blobs
[X, Y] = meshgrid(1:nX, 1:nY);
blobA = exp(-((X-20).^2 + (Y-32).^2) / (2*8^2));   % left blob
blobB = exp(-((X-44).^2 + (Y-32).^2) / (2*8^2));   % right blob
baseline = 0.3 * blobA + 0.3 * blobB + 0.02 * rand(nY, nX);  % resting intensity

imIn = zeros(nY, nX, 1, 1, nT, 'single');
tau = 4;   % recovery time constant (frames)
for iT = 1:nT
    if iT < 6
        frame = baseline;
    else
        t = iT - 6;
        dA = 0.8 * exp(-t / tau);   % 80% step, decaying in blob A
        dB = 0.4 * exp(-t / tau);   % 40% step, decaying in blob B
        frame = baseline + dA * blobA + dB * blobB;
    end
    imIn(:,:,1,1,iT) = single(frame) + 0.01 * rand(nY, nX, 'single');
end
fprintf('  → Done. Baseline mean: %.3f; peak frame mean: %.3f\n', ...
    mean(imIn(:,:,1,1,1), 'all'), max(mean(reshape(imIn,[],nT))));

% -------------------------------------------------------------------------
% Section 2 — Compare variance, CV, dF/F methods
% -------------------------------------------------------------------------
fprintf('\nSection 2: Computing variance, CV, and dF/F...\n');

p_var = struct('channelIdx',1, 'method','variance', ...
               'baseline', struct('frames',[1 5]));
p_cv  = struct('channelIdx',1, 'method','cv', ...
               'baseline', struct('frames',[1 5]));
p_dff = struct('channelIdx',1, 'method','dff', ...
               'baseline', struct('frames',[1 5]));

[~, varMean, ~]  = intensityFluctuationCompute(imIn, p_var);
[~, cvMean,  ~]  = intensityFluctuationCompute(imIn, p_cv);
[dffSeries, dffMean, F0] = intensityFluctuationCompute(imIn, p_dff);

figure('Name','Section 2: Method Comparison','Position',[50 400 900 280]);
colormap('hot');

subplot(1,3,1);
imagesc(varMean(:,:,1,1));  axis image off; colorbar;
title('Variance'); clim([0 max(varMean(:))]);
xlabel('Higher = more variable');

subplot(1,3,2);
imagesc(cvMean(:,:,1,1));   axis image off; colorbar;
title('CV (σ/μ)'); clim([0 max(cvMean(:))]);
xlabel('Normalised for intensity');

subplot(1,3,3);
imagesc(dffMean(:,:,1,1));  axis image off; colorbar;
title('Mean ΔF/F₀');
xlabel('Baseline-normalised response');

fprintf('  → Variance peak: %.4f | CV peak: %.4f | Mean dF/F peak: %.2f\n', ...
    max(varMean(:)), max(cvMean(:)), max(dffMean(:)));

% -------------------------------------------------------------------------
% Section 3 — Effect of baseline window choice on dF/F
% -------------------------------------------------------------------------
fprintf('\nSection 3: Comparing baseline window choices...\n');

% Correct window (pre-stimulus): frames 1-5
p_correct = struct('channelIdx',1,'method','dff','baseline',struct('frames',[1 5]));
[~, dffCorrect, ~] = intensityFluctuationCompute(imIn, p_correct);

% Wrong window (includes stimulus): frames 1-10
p_wrong = struct('channelIdx',1,'method','dff','baseline',struct('frames',[1 10]));
[~, dffWrong, ~] = intensityFluctuationCompute(imIn, p_wrong);

% Time course at peak pixel of blob A
[~, pkIdx] = max(dffCorrect(:));
[pyA, pxA, ~, ~] = ind2sub([nY nX 1 1], pkIdx);

tAxis = 1:nT;
dffCurveCorrect = squeeze(dffSeries(pyA, pxA, 1, 1, :));
dffSeries2 = intensityFluctuationCompute(imIn, p_wrong);
dffCurveWrong = squeeze(dffSeries2(pyA, pxA, 1, 1, :));

figure('Name','Section 3: Baseline Window Effect','Position',[50 60 700 280]);

subplot(1,2,1);
plot(tAxis, dffCurveCorrect, 'b-o', 'LineWidth',1.5, 'MarkerSize',4); hold on;
plot(tAxis, dffCurveWrong,   'r--s','LineWidth',1.5, 'MarkerSize',4);
xline(5.5, 'k--'); xline(10.5, 'r:');
xlabel('Frame'); ylabel('ΔF/F₀');
legend('Correct window (1-5)','Inflated window (1-10)','Stimulus onset','Location','ne');
title('Peak pixel ΔF/F₀ time course');

subplot(1,2,2);
imagesc(dffCorrect(:,:,1,1) - dffWrong(:,:,1,1)); axis image off; colorbar;
title('Correct − Inflated mean ΔF/F₀');
xlabel('Positive = baseline window underestimates response');

fprintf('  → Peak dF/F (correct baseline): %.2f; (inflated baseline): %.2f\n', ...
    max(dffCurveCorrect), max(dffCurveWrong));

% -------------------------------------------------------------------------
% Section 4 — Visualise dF/F with physiologyConvertToRgb
% -------------------------------------------------------------------------
fprintf('\nSection 4: RGB visualisation of peak dF/F frame...\n');

% Find the frame with maximum mean dF/F
[~, peakFrame] = max(mean(reshape(dffSeries, [], nT)));

% Use the intensity image as background, erode an ER mask from blob union
erMask = repmat(single(blobA + blobB) > 0.05, [1 1 1 1 nT]);

% Colormap: sequential (parula or colorcet R2 if available)
if exist('colorcet','file')
    cMap = colorcet('R2');
else
    cMap = parula(256);
end

dffFrame = dffSeries(:,:,:,:,peakFrame);
intFrame = imIn(:,:,:,:,peakFrame);

rgb = physiologyConvertToRgb(dffFrame, intFrame, erMask(:,:,:,:,peakFrame), ...
    0, max(dffFrame(:))*1.05, true, cMap);

figure('Name','Section 4: dF/F RGB (peak frame)','Position',[800 400 560 480]);
imshow(squeeze(rgb(:,:,1,1,:)));
title(sprintf('ΔF/F₀ map — frame %d (peak response)', peakFrame));
colormap(cMap); c = colorbar;
c.Label.String = 'ΔF/F₀';
clim([0, max(dffFrame(:))*1.05]);

fprintf('  → Peak frame: %d; max ΔF/F₀ in frame: %.2f\n', ...
    peakFrame, max(dffFrame(:)));
fprintf('\n=== demoIntensityFluctuation complete ===\n');
