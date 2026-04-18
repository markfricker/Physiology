%DEMOFRAP  Demo: FRAP and FLIP analysis on synthetic ER network data.
%
% Demonstrates frapCompute on synthetic confocal time-series data where
% a photobleach event is applied to a defined circular region:
%   Section 1 – Generate synthetic ER network with FRAP bleach event
%   Section 2 – Run frapCompute with explicit and auto-detected bleach frames
%   Section 3 – Visualise normalised recovery series (key frames)
%   Section 4 – Region-averaged recovery curves from fitResults
%   Section 5 – Effect of diffusion rate on recovery shape
%   Section 6 – FLIP variant: fluorescence loss in unbleached region
%
% Run with:
%   demoFrap
%
% Path requirements:
%   Physiology_sandbox/src/ and Physiology_sandbox/utils/ must be on path.

close all; clc;
fprintf('=== demoFrap ===\n\n');

% =========================================================================
% Section 1 — Synthetic ER network with FRAP bleach event
% =========================================================================
fprintf('Section 1: Generating synthetic FRAP data...\n');

nY = 80; nX = 80; nT = 30; bleachFrame = 8;
Fpre    = 0.75;   % pre-bleach intensity
Fbleach = 0.05;   % post-bleach intensity in bleach zone
tau     = 6;      % recovery time constant (frames) → mobile fraction ~1.0

% Spatial structure: thin tubular ER network approximated by a grid of lines
imBg = zeros(nY, nX, 'single');
for y = 8:12:nY, imBg(y, :) = 1; end   % horizontal tubules
for x = 8:12:nX, imBg(:, x) = 1; end   % vertical tubules (junctions)
imBg = imgaussfilt(imBg, 0.8);           % smooth
imBg = imBg / max(imBg(:));             % normalise to [0,1]

% ER mask and skeleton
erMask5 = repmat(imBg > 0.15, [1 1 1 1 nT]);
erSkel  = repmat(bwmorph(imBg > 0.4, 'thin', Inf), [1 1 1 1 nT]);

% Bleach ROI: circle centred at (40,40) radius 10
[Xg, Yg] = meshgrid(1:nX, 1:nY);
bleachROI = (Xg-40).^2 + (Yg-40).^2 <= 10^2;

% Generate time series
imIn = zeros(nY, nX, 1, 1, nT, 'single');
for iT = 1:nT
    frame = Fpre * imBg;
    if iT >= bleachFrame
        t  = iT - bleachFrame;
        Ft = Fbleach + (Fpre - Fbleach) * (1 - exp(-t / tau));
        frame(bleachROI) = Ft * imBg(bleachROI);
    end
    imIn(:,:,1,1,iT) = frame + 0.008 * randn(nY, nX, 'single');
end

fprintf('  → Network: %d×%d px, %d frames, bleach at frame %d\n', ...
    nY, nX, nT, bleachFrame);

% =========================================================================
% Section 2 — Run frapCompute
% =========================================================================
fprintf('\nSection 2: Running frapCompute (explicit and auto-detected bleach)...\n');

% a) Explicit bleach frame
p_explicit = struct('channelIdx',1,'method','frap', ...
                    'bleachFrame',bleachFrame,'preFrames',0);
[recExp, recMeanExp, preMap, fitRes] = frapCompute(imIn, erSkel, ...
    false(nY,nX,1,1,nT), p_explicit);

% b) Auto-detected bleach frame (p.bleachFrame = 0)
p_auto = struct('channelIdx',1,'method','frap','bleachFrame',0,'preFrames',0);
[~, ~, ~, fitResAuto] = frapCompute(imIn, erSkel, false(nY,nX,1,1,nT), p_auto);

fprintf('  → Explicit bleach frame: %d\n', fitRes.bleachFrame);
fprintf('  → Auto-detected bleach frame: %d (true: %d)\n', ...
    fitResAuto.bleachFrame, bleachFrame);
if fitResAuto.bleachFrame == bleachFrame
    fprintf('  → Auto-detection correct!\n');
else
    fprintf('  → Auto-detection offset by %d frame(s)\n', ...
        abs(fitResAuto.bleachFrame - bleachFrame));
end

% =========================================================================
% Section 3 — Recovery series: key frames
% =========================================================================
fprintf('\nSection 3: Plotting key frames of normalised recovery...\n');

framesToShow = [bleachFrame, bleachFrame+2, bleachFrame+5, bleachFrame+10, nT];
framesToShow = framesToShow(framesToShow <= nT);
nShow = numel(framesToShow);

figure('Name','Section 3: FRAP Recovery Series','Position',[50 400 180*nShow 240]);
for k = 1:nShow
    iT = framesToShow(k);
    subplot(1, nShow, k);
    recFrame = recExp(:,:,1,1,iT);
    imagesc(recFrame, [0 1]); axis image off; colormap(hot);
    viscircles([40 40], 10, 'Color','c','LineWidth',0.8);
    if iT < bleachFrame
        title(sprintf('t=%d (pre)', iT));
    else
        title(sprintf('t=%d (+%d)', iT, iT-bleachFrame));
    end
end
colorbar('Position',[0.92 0.15 0.02 0.7]);
sgtitle('Normalised recovery F_{norm} (pre=1, bleach=0)');

% =========================================================================
% Section 4 — Region-averaged recovery curves
% =========================================================================
fprintf('\nSection 4: Region-averaged recovery curves from fitResults...\n');

tAxis = (1:nT) - bleachFrame;   % time relative to bleach (frames)
skCurve = squeeze(fitRes.skeletonCurve(1,1,:));

% Also compute whole-ROI curve manually: mean over bleach ROI pixels
roiCurve = zeros(1, nT, 'single');
for iT = 1:nT
    roiPixels = recExp(:,:,1,1,iT);
    roiCurve(iT) = mean(roiPixels(bleachROI), 'omitnan');
end

% Fit single exponential to post-bleach skeleton curve
postIdx   = bleachFrame:nT;
tFit      = tAxis(postIdx);
curveFit  = double(skCurve(postIdx));
% Simple analytical: Mf*(1-exp(-t/tau)) fitting via least squares
% Use fixed tau=6 for display; full fitting would use frapFit
MfEst = max(curveFit);
tauEst = tau;
fitCurve = MfEst * (1 - exp(-tFit / tauEst));

figure('Name','Section 4: Recovery Curves','Position',[50 80 700 320]);
subplot(1,2,1);
plot(tAxis, skCurve, 'b-o', 'LineWidth',2, 'MarkerSize',4); hold on;
plot(tAxis, roiCurve,'r--s','LineWidth',2,'MarkerSize',4);
plot(tFit, fitCurve, 'k-', 'LineWidth',1.5);
xline(0,'k:','Bleach');
xlabel('Time relative to bleach (frames)');
ylabel('Normalised fluorescence');
legend('Skeleton ROI (fitResults)','Bleach zone mean','Exp fit','Location','se');
title('Region-averaged recovery');
ylim([-0.1 1.15]); grid on;

subplot(1,2,2);
imagesc(recMeanExp(:,:,1,1),[0 1]); axis image off; colorbar;
viscircles([40 40],10,'Color','c','LineWidth',1);
colormap(gca,'hot');
title('Mean post-bleach recovery map');

fprintf('  → Skeleton curve: bleach depth %.2f, final recovery %.2f\n', ...
    min(skCurve(postIdx)), max(skCurve));
fprintf('  → Estimated mobile fraction: %.2f, τ = %d frames\n', MfEst, tauEst);

% =========================================================================
% Section 5 — Effect of diffusion rate (tau) on recovery
% =========================================================================
fprintf('\nSection 5: Recovery curves for different diffusion rates...\n');

taus   = [2, 4, 8, 16];   % recovery time constants (frames)
colors = lines(numel(taus));

figure('Name','Section 5: Diffusion Rate Effect','Position',[800 400 560 340]);
hold on;
for k = 1:numel(taus)
    imFast = zeros(nY, nX, 1, 1, nT, 'single');
    for iT = 1:nT
        frame = Fpre * imBg;
        if iT >= bleachFrame
            t = iT - bleachFrame;
            Ft = Fbleach + (Fpre - Fbleach) * (1 - exp(-t / taus(k)));
            frame(bleachROI) = Ft * imBg(bleachROI);
        end
        imFast(:,:,1,1,iT) = frame + 0.005 * randn(nY,nX,'single');
    end
    p_fast = struct('channelIdx',1,'method','frap', ...
                    'bleachFrame',bleachFrame,'preFrames',0);
    [recFast,~,~,fr_fast] = frapCompute(imFast, erSkel, false(nY,nX,1,1,nT), p_fast);
    skC_k = squeeze(fr_fast.skeletonCurve(1,1,:));
    plot(tAxis, skC_k, '-o', 'Color',colors(k,:), 'LineWidth',1.5, ...
        'MarkerSize',3, 'DisplayName', sprintf('\\tau = %d frames', taus(k)));
end
xline(0,'k:','Bleach'); yline(1,'k--','Full recovery');
xlabel('Time relative to bleach (frames)');
ylabel('Normalised fluorescence (skeleton ROI)');
title('FRAP recovery: effect of diffusion rate (τ)');
legend('Location','se'); grid on; ylim([-0.15 1.2]);

t_half_expected = taus * log(2);
fprintf('  → Expected t½ values (frames): ');
fprintf('%.1f  ', t_half_expected); fprintf('\n');

% =========================================================================
% Section 6 — FLIP variant
% =========================================================================
fprintf('\nSection 6: FLIP — fluorescence loss in remote region...\n');

% In FLIP, a remote region is repeatedly bleached; the monitored region
% shows progressive fluorescence loss.  Simulate: whole-image bleaching
% except one quadrant, which shows loss due to connectivity.
imFLIP = zeros(nY, nX, 1, 1, nT, 'single');
tauFlip = 8;
monitorROI = Xg > 60 & Yg < 40;   % top-right quadrant = monitored
for iT = 1:nT
    frame = Fpre * imBg;
    if iT >= bleachFrame
        t = iT - bleachFrame;
        % Whole image loses fluorescence (repeatedly bleached), monitor region lags
        Fbulk    = Fpre * exp(-t / tauFlip);      % bulk loss (Fpre decays to 0)
        Fmonitor = Fpre * exp(-0.4 * t / tauFlip); % monitor region lags (partial connectivity)
        frame          = Fbulk    * imBg;
        frame(monitorROI) = Fmonitor * imBg(monitorROI);
    end
    imFLIP(:,:,1,1,iT) = frame + 0.008 * randn(nY,nX,'single');
end

p_flip = struct('channelIdx',1,'method','flip', ...
                'bleachFrame',bleachFrame,'preFrames',0);
[recFlip, ~, ~, frFlip] = frapCompute(imFLIP, erSkel, false(nY,nX,1,1,nT), p_flip);

figure('Name','Section 6: FLIP','Position',[800 80 700 300]);
subplot(1,3,1);
imagesc(recFlip(:,:,1,1,bleachFrame+2),[0 1]); axis image off; colorbar;
colormap(gca,'hot'); title('FLIP: frame +2'); colorbar;

subplot(1,3,2);
imagesc(recFlip(:,:,1,1,bleachFrame+10),[0 1]); axis image off; colorbar;
colormap(gca,'hot'); title('FLIP: frame +10');

subplot(1,3,3);
skFlip = squeeze(frFlip.skeletonCurve(1,1,:));
plot(tAxis, skFlip, 'b-o','LineWidth',2,'MarkerSize',4); hold on;
xline(0,'k:','Bleach'); yline(0,'k--');
xlabel('Frames relative to bleach'); ylabel('Normalised fluorescence');
title('FLIP: skeleton recovery curve');
grid on;
fprintf('  → FLIP: fluorescence in remote region after %d frames: %.2f\n', ...
    nT - bleachFrame, min(skFlip(bleachFrame:end)));
fprintf('\n=== demoFrap complete ===\n');
