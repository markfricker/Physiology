function [fluctuationSeries, fluctuationMean, baselineMap] = ...
        intensityFluctuationCompute(imIn, p)
%INTENSITYFLUCTUATIONCOMPUTE  Per-pixel temporal intensity statistics for a single channel.
%
%   [fluctuationSeries, fluctuationMean, baselineMap] = ...
%       intensityFluctuationCompute(imIn, p)
%
% Computes per-pixel temporal statistics from a single fluorescence channel.
% All methods produce a per-frame series and a spatial summary map.  The
% baseline map (F₀) is always returned so that the calling context can store
% it independently for display or further analysis.
%
% INPUTS
%   imIn – [nY nX nC nZ nT] numeric; intensity image.
%   p    – parameter struct:
%            p.channelIdx     – scalar integer; channel to process.
%                               Output arrays span all nC input channels but
%                               only channelIdx is computed; others remain zero.
%            p.method         – 'variance' | 'cv' | 'dff'
%                               'variance': temporal variance per pixel.
%                               'cv':       coefficient of variation (σ/μ).
%                               'dff':      ΔF/F₀ per frame; main output for
%                                           calcium indicators.
%            p.baseline.frames – [1×2] integer; [startFrame endFrame] for the
%                               F₀ baseline window.  Default: first 10% of
%                               frames (minimum 1 frame).
%
% OUTPUTS
%   fluctuationSeries – [nY nX nC nZ nT] single; per-frame result.
%                       'dff':             ΔF/F₀ = (F(t)−F₀)/F₀ time series.
%                       'variance'/'cv':   F(t)/F₀ (baseline-normalised
%                                          intensity) at each frame.
%                       Only p.channelIdx channel is populated; others zero.
%   fluctuationMean   – [nY nX nC nZ] single; spatial summary map.
%                       'dff':      temporal mean of ΔF/F₀.
%                       'variance': temporal variance of raw intensity.
%                       'cv':       coefficient of variation (σ/μ).
%   baselineMap       – [nY nX nC nZ] single; F₀ baseline (mean intensity
%                       over p.baseline.frames window).
%
% NOTE ON METHODS
%   'dff' is the standard calcium-imaging normalisation.  The baseline window
%   should span a quiescent period before stimulation.  For variance and CV
%   the baseline window is still used to compute F₀ for the series output,
%   but the summary map is computed over the full time series.
%
%   Pixels with F₀ ≤ 0 (or very close to zero) are treated as background:
%   their series values are set to zero to avoid division artefacts.

[nY, nX, nC, nZ, nT] = size(imIn);
cIdx   = sf(p, 'channelIdx', 1);
method = lower(sf(p, 'method', 'dff'));

% Baseline window -- default to first 10 % of frames
defaultEnd    = max(1, round(nT / 10));
baselineStruct = sf(p, 'baseline', struct('frames', [1 defaultEnd]));
bf            = baselineStruct.frames;
bf(1)         = max(1, min(nT, bf(1)));
bf(2)         = max(bf(1), min(nT, bf(2)));

imIn = im2single(imIn);

fluctuationSeries = zeros(nY, nX, nC, nZ, nT, 'single');
fluctuationMean   = zeros(nY, nX, nC, nZ,     'single');
baselineMap       = zeros(nY, nX, nC, nZ,     'single');

for iZ = 1:nZ

    % Extract single-channel Z-plane: [nY nX nT]
    im = reshape(imIn(:,:,cIdx,iZ,:), nY, nX, nT);

    % ---- baseline F₀ ----------------------------------------------------
    F0   = mean(im(:,:,bf(1):bf(2)), 3);      % [nY nX]
    mask = F0 > eps('single');                 % background suppression mask
    F0safe = F0;
    F0safe(~mask) = 1;                         % avoid /0; series zeroed below

    % ---- per-frame baseline-normalised series ---------------------------
    Fnorm = im ./ F0safe;                      % F(t)/F₀,  [nY nX nT]

    % ---- method ---------------------------------------------------------
    switch method

        case 'dff'
            series = Fnorm - 1;                % ΔF/F₀
            series(:,:) = series(:,:) .* mask; % zero background pixels
            sumMap = mean(series, 3);

        case 'variance'
            series = Fnorm .* mask;            % normalised series
            sumMap = var(im, 0, 3);            % summary is raw variance

        case 'cv'
            series = Fnorm .* mask;
            mu  = mean(im, 3);
            sig = std(im,  0, 3);
            mu(mu <= eps('single')) = 1;       % avoid /0
            sumMap = sig ./ mu;
            sumMap(~mask) = 0;

        otherwise
            error('intensityFluctuationCompute:unknownMethod', ...
                'Unknown method ''%s''. Use ''dff'', ''variance'', or ''cv''.', ...
                method);
    end

    % Zero background in summary map
    sumMap(~mask) = 0;

    % ---- store back into 5-D arrays -------------------------------------
    fluctuationSeries(:,:,cIdx,iZ,:) = reshape(series, nY, nX, 1, 1, nT);
    fluctuationMean(:,:,cIdx,iZ)     = sumMap;
    baselineMap(:,:,cIdx,iZ)         = F0;
end

end % intensityFluctuationCompute


% =========================================================================
function v = sf(s, field, default)
%SF  Safe field access with default.
if isstruct(s) && isfield(s, field)
    v = s.(field);
else
    v = default;
end
end
