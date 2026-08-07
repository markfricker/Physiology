function [recoverySeries, recoveryMean, prebleachMap, fitResults] = ...
        frapCompute(imIn, skIn, cisternaIn, p)
%FRAPCOMPUTE  Per-pixel fluorescence recovery/loss analysis (FRAP / FLIP).
%
%   [recoverySeries, recoveryMean, prebleachMap, fitResults] = ...
%       frapCompute(imIn, skIn, cisternaIn, p)
%
% Computes normalised fluorescence dynamics following a photobleach event.
% For FRAP, fluorescence recovers into the bleached region; for FLIP,
% fluorescence is progressively lost from a remote monitored region.
%
% Per-pixel normalisation produces spatial maps suitable for display via
% the AnalyzER overlay system.  Region-averaged recovery curves over the
% ER skeleton and cisternae masks are returned in fitResults for downstream
% quantification (curve fitting is not performed here; see frapFit).
%
% INPUTS
%   imIn       – [nY nX nC nZ nT] numeric; intensity image.
%   skIn       – [nY nX nC nZ nT] logical; ER skeleton mask.
%                May have fewer C/Z dimensions (broadcast).
%   cisternaIn – [nY nX nC nZ nT] logical; cisternae mask.
%                May have fewer C/Z dimensions (broadcast).
%   p          – parameter struct:
%                  p.channelIdx  – scalar integer; channel to process.
%                  p.method      – 'frap' | 'flip' (default 'frap').
%                                  Both methods produce the same normalised
%                                  series; 'flip' annotations fitResults.
%                  p.bleachFrame – integer; frame index of the bleach event.
%                                  Set to 0 for automatic detection (default).
%                  p.preFrames   – integer; number of pre-bleach frames used
%                                  to compute F_pre baseline.  Default: all
%                                  frames before the bleach frame.
%
% OUTPUTS
%   recoverySeries – [nY nX nC nZ nT] single; per-pixel normalised
%                    fluorescence intensity.  Normalisation:
%                      F_norm(t) = (F(t) − F_bg) / (F_pre − F_bg)
%                    where F_pre = mean over pre-bleach frames,
%                    F_bg = per-pixel minimum over the full time series.
%                    Values approach 1 at full recovery (FRAP) or drop
%                    toward 0 with complete loss (FLIP).  Pre-bleach frames
%                    are set to 1; pixels with F_pre ≤ F_bg are NaN.
%   recoveryMean   – [nY nX nC nZ] single; mean of post-bleach recoverySeries
%                    (frames bleachFrame:end).
%   prebleachMap   – [nY nX nC nZ] single; F_pre baseline map (mean intensity
%                    over pre-bleach frames, before background subtraction).
%   fitResults     – struct with fields:
%                      .bleachFrame    – integer; detected or supplied bleach frame.
%                      .method         – 'frap' or 'flip'.
%                      .skeletonCurve  – [nCh×nZ×nT] single; mean recoverySeries
%                                        over skeleton ROI per channel/Z (NaN
%                                        where skeleton is empty).
%                      .cisternaCurve  – [nCh×nZ×nT] single; mean over cisternae ROI.
%                    Pass skeletonCurve / cisternaCurve to frapFit for
%                    exponential curve fitting and mobile-fraction extraction.
%
% NOTE ON BLEACH DETECTION
%   Automatic detection (p.bleachFrame = 0) identifies the bleach frame as
%   the frame with the largest drop in mean image intensity from one frame
%   to the next.  This is robust for standard FRAP acquisitions where the
%   bleach event causes a large step decrease.  For FLIP with repeated
%   low-power bleaches, supply bleachFrame explicitly.
%
% NOTE ON FRAP vs FLIP
%   The same normalisation is applied for both methods.  The distinction is
%   semantic and noted in fitResults.method.  For FLIP, the recovery curves
%   will be monotonically decreasing rather than increasing.

[nY, nX, nC, nZ, nT] = size(imIn);
[~, ~, skC, skZ, ~]  = size(skIn);
[~, ~, ciC, ciZ, ~]  = size(cisternaIn);

cIdx      = sf(p, 'channelIdx',  1);
method    = lower(sf(p, 'method',      'frap'));
bfInput   = sf(p, 'bleachFrame', 0);
preFrames = sf(p, 'preFrames',   0);       % 0 = use all frames before bleach

imIn = im2single(imIn);

% ---- detect or validate bleach frame ------------------------------------
if bfInput == 0
    bf = detectBleachFrame(imIn, cIdx);
else
    bf = round(bfInput);
end
bf = max(2, min(nT, bf));                  % clamp: need ≥1 pre-bleach frame

% ---- pre-bleach window --------------------------------------------------
if preFrames <= 0
    preWin = 1:bf-1;
else
    preWin = max(1, bf - preFrames) : bf-1;
end

nCh = 1;                                   % single channel
Cidx_arr = cIdx;

% NaN-initialised (not zero) -- only cIdx gets computed, and a stored zero
% for the other channels would be indistinguishable from a genuine
% zero-recovery measurement to every downstream consumer.
recoverySeries = nan(nY, nX, nC, nZ, nT, 'single');
recoveryMean   = nan(nY, nX, nC, nZ,     'single');
prebleachMap   = nan(nY, nX, nC, nZ,     'single');

skCurv = NaN(nCh, nZ, nT, 'single');
ciCurv = NaN(nCh, nZ, nT, 'single');

for iZ = 1:nZ
    iCsk = min(cIdx, skC);
    iCci = min(cIdx, ciC);
    iZsk = min(iZ,   skZ);
    iZci = min(iZ,   ciZ);

    im   = reshape(imIn(:,:,cIdx,iZ,:), nY, nX, nT);  % [nY nX nT]

    % ---- background: per-pixel minimum over full series -----------------
    Fbg = min(im, [], 3);                              % [nY nX]

    % ---- pre-bleach baseline --------------------------------------------
    Fpre = mean(im(:,:,preWin), 3);                    % [nY nX]

    % Pixels where Fpre ≤ Fbg have no usable signal → NaN
    span = Fpre - Fbg;                                 % [nY nX]
    validMask = span > eps('single');
    spanSafe  = span;
    spanSafe(~validMask) = 1;                          % avoid /0; NaN-ed below

    % ---- per-pixel normalised series ------------------------------------
    Fnorm = (im - Fbg) ./ spanSafe;                    % [nY nX nT]
    % Use 3-D linear indexing to set background pixels to NaN.
    % IMPORTANT: do NOT use 2-subscript indexing (e.g. Fnorm(~mask, :)) on a
    % 3-D array — MATLAB collapses it to 2-D for that call and the assignment
    % permanently reshapes Fnorm, breaking the reshape at the store-back step.
    invMask3D = repmat(~validMask, [1 1 nT]);           % [nY nX nT] logical
    Fnorm(invMask3D) = NaN;

    % Pre-bleach frames clamped to 1 (by definition fully recovered before bleach)
    Fnorm(:,:,preWin) = 1;

    % ---- region-averaged curves over ER masks ---------------------------
    sk   = squeeze(skIn(:,:,iCsk,iZsk,:));             % [nY nX nT] or [nY nX]
    cist = squeeze(cisternaIn(:,:,iCci,iZci,:));

    if ndims(sk) == 2,   sk   = repmat(sk,   [1 1 nT]); end  %#ok<ISMAT>
    if ndims(cist) == 2, cist = repmat(cist, [1 1 nT]); end  %#ok<ISMAT>

    for iT = 1:nT
        skRoi   = sk(:,:,iT);
        ciRoi   = cist(:,:,iT);
        frame   = Fnorm(:,:,iT);

        skVals  = frame(skRoi);
        ciVals  = frame(ciRoi);

        if ~isempty(skVals),  skCurv(1,iZ,iT) = mean(skVals, 'omitnan'); end
        if ~isempty(ciVals),  ciCurv(1,iZ,iT) = mean(ciVals, 'omitnan'); end
    end

    % ---- summary (mean post-bleach) -------------------------------------
    recoverySeries(:,:,cIdx,iZ,:) = reshape(Fnorm, nY, nX, 1, 1, nT);
    recoveryMean(:,:,cIdx,iZ)     = mean(Fnorm(:,:,bf:end), 3, 'omitnan');
    prebleachMap(:,:,cIdx,iZ)     = Fpre;
end

fitResults.bleachFrame   = bf;
fitResults.method        = method;
fitResults.skeletonCurve = skCurv;
fitResults.cisternaCurve = ciCurv;

end % frapCompute


% =========================================================================
function bf = detectBleachFrame(imIn, cIdx)
%DETECTBLEACHFRAME  Find the frame with the largest single-step intensity drop.
%
% Averages over the nominated channel (cIdx) and all Z planes, then finds
% the largest negative frame-to-frame difference.

[~,~,~,nZ,nT] = size(imIn);

meanPerFrame = zeros(1, nT, 'single');
for iT = 1:nT
    slice = imIn(:,:,cIdx,:,iT);
    meanPerFrame(iT) = mean(slice(:), 'omitnan');
end

diff_ = diff(meanPerFrame);           % length nT-1; negative = intensity drop
[~, idx] = min(diff_);               % most negative = steepest drop
bf = idx + 1;                        % bleach frame is the one AFTER the drop

end % detectBleachFrame


% =========================================================================
function v = sf(s, field, default)
%SF  Safe field access with default.
if isstruct(s) && isfield(s, field)
    v = s.(field);
else
    v = default;
end
end
