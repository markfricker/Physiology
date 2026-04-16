function [ratioSeries, ratioMean, calibratedSeries, calibratedMean] = ...
        ratioCompute(imIn, p)
%RATIOCOMPUTE  Per-pixel ratiometric measurement from two fluorescence channels.
%
%   [ratioSeries, ratioMean, calibratedSeries, calibratedMean] = ...
%       ratioCompute(imIn, p)
%
% Computes the per-pixel intensity ratio of a numerator channel divided by
% a denominator channel.  An optional third channel corrects for cellular
% autofluorescence before ratioing.  A calibration step converts the raw
% ratio to a physical quantity (e.g. pH) using one of four models.
%
% INPUTS
%   imIn – [nY nX nC nZ nT] numeric; intensity image.
%   p    – parameter struct:
%
%          p.numeratorChannel   – integer; channel index for numerator.
%          p.denominatorChannel – integer; channel index for denominator.
%
%          p.autofluorescence.use             – logical (default false).
%          p.autofluorescence.channel         – integer; AF channel index.
%          p.autofluorescence.scaleNumerator  – scalar ≥ 0 (default 1);
%                fraction of AF signal subtracted from numerator channel.
%          p.autofluorescence.scaleDenominator – scalar ≥ 0 (default 1);
%                fraction of AF signal subtracted from denominator channel.
%
%          p.calibration.method – 'none' | 'linear' | 'polynomial' |
%                                  'henderson-hasselbalch' (default 'none').
%
%          Linear:   calibrated = slope * ratio + intercept
%            p.calibration.slope      (default 1)
%            p.calibration.intercept  (default 0)
%
%          Polynomial:  calibrated = polyval(coefficients, ratio)
%            p.calibration.coefficients – [1×N] coefficient vector,
%                                         highest power first.
%
%          Henderson-Hasselbalch:  calibrated = pKa + log10((R-Racid)/(Rbase-R))
%            p.calibration.pKa    – indicator pKa (default 7.0)
%            p.calibration.Racid  – ratio at acid endpoint
%            p.calibration.Rbase  – ratio at base endpoint
%
% OUTPUTS
%   ratioSeries      – [nY nX 1 nZ nT] single; per-frame ratio R = num/den.
%                      Pixels where den ≤ 0 after AF correction are NaN.
%   ratioMean        – [nY nX 1 nZ] single; time-averaged ratio.
%   calibratedSeries – [nY nX 1 nZ nT] single; calibrated quantity per frame.
%                      All-NaN when p.calibration.method = 'none'.
%   calibratedMean   – [nY nX 1 nZ] single; time-averaged calibrated value.
%
% NOTE ON AUTOFLUORESCENCE CORRECTION
%   The autofluorescence channel should be one where the ratiometric probe
%   produces no signal (e.g. a spectrally distinct reference channel, or a
%   control wavelength).  The correction is:
%       corrected_num(t) = max(0,  num(t) − scaleNumerator   × AF(t))
%       corrected_den(t) = max(0,  den(t) − scaleDenominator × AF(t))
%   Scale factors can be estimated from cells expressing no probe (AF only).
%
% NOTE ON OUTPUT CHANNEL DIMENSION
%   The ratio is a single derived quantity; the output nC dimension is 1.
%   The AnalyzER run wrapper stores the result under a separate image field.

[nY, nX, ~, nZ, nT] = size(imIn);
imIn = im2single(imIn);

nCh  = sf(p, 'numeratorChannel',   1);
dCh  = sf(p, 'denominatorChannel', 2);

afStruct  = sf(p,  'autofluorescence', struct());
afUse     = sf(afStruct, 'use',                false);
afCh      = sf(afStruct, 'channel',            3);
afScaleN  = sf(afStruct, 'scaleNumerator',     1);
afScaleD  = sf(afStruct, 'scaleDenominator',   1);

calStruct = sf(p, 'calibration', struct());
calMethod = lower(sf(calStruct, 'method', 'none'));

ratioSeries      = zeros(nY, nX, 1, nZ, nT, 'single');
ratioMean        = zeros(nY, nX, 1, nZ,     'single');
calibratedSeries = NaN  (nY, nX, 1, nZ, nT, 'single');
calibratedMean   = NaN  (nY, nX, 1, nZ,     'single');

for iZ = 1:nZ

    num = reshape(imIn(:,:,nCh,iZ,:), nY, nX, nT);  % [nY nX nT]
    den = reshape(imIn(:,:,dCh,iZ,:), nY, nX, nT);

    % ---- autofluorescence correction ------------------------------------
    if afUse
        af  = reshape(imIn(:,:,afCh,iZ,:), nY, nX, nT);
        num = max(0, num - afScaleN * af);
        den = max(0, den - afScaleD * af);
    end

    % ---- ratio ----------------------------------------------------------
    % Pixels where denominator ≤ 0 after correction are undefined → NaN.
    denSafe = den;
    badPx   = den <= eps('single');
    denSafe(badPx) = 1;                           % avoid /0; NaN-ed below

    R = num ./ denSafe;                           % [nY nX nT]
    R(badPx) = NaN;

    Rmean = mean(R, 3, 'omitnan');                % [nY nX]

    ratioSeries(:,:,1,iZ,:) = reshape(R,     nY, nX, 1, 1, nT);
    ratioMean(:,:,1,iZ)     = Rmean;

    % ---- calibration ----------------------------------------------------
    if ~strcmp(calMethod, 'none')
        C     = calibrateRatio(R,     calMethod, calStruct);  % [nY nX nT]
        Cmean = calibrateRatio(Rmean, calMethod, calStruct);  % [nY nX]

        calibratedSeries(:,:,1,iZ,:) = reshape(C, nY, nX, 1, 1, nT);
        calibratedMean(:,:,1,iZ)     = Cmean;
    end

end

end % ratioCompute


% =========================================================================
function out = calibrateRatio(R, method, cal)
%CALIBRATERATIO  Apply calibration model to ratio values.
%
% R   – any shape single/double array of ratio values.
% out – same shape as R; NaN where R is NaN or calibration is undefined.

switch method

    case 'linear'
        slope     = sf(cal, 'slope',     1);
        intercept = sf(cal, 'intercept', 0);
        out = slope .* R + intercept;

    case 'polynomial'
        coeffs = sf(cal, 'coefficients', [1 0]);
        out    = polyval(coeffs, R);

    case 'henderson-hasselbalch'
        pKa   = sf(cal, 'pKa',   7.0);
        Racid = sf(cal, 'Racid', NaN);
        Rbase = sf(cal, 'Rbase', NaN);

        if isnan(Racid) || isnan(Rbase)
            error('ratioCompute:missingCalibration', ...
                'Henderson-Hasselbalch calibration requires p.calibration.Racid and Rbase.');
        end
        if Rbase <= Racid
            error('ratioCompute:calibrationRange', ...
                'Rbase (%g) must be greater than Racid (%g).', Rbase, Racid);
        end

        num_hh = R - Racid;
        den_hh = Rbase - R;
        out = pKa + log10(num_hh ./ den_hh);

        % Ratio outside calibration range is physically undefined → NaN
        out(R <= Racid | R >= Rbase) = NaN;

    otherwise
        error('ratioCompute:unknownCalibration', ...
            'Unknown calibration method ''%s''.', method);
end

end % calibrateRatio


% =========================================================================
function v = sf(s, field, default)
%SF  Safe field access with default.
if isstruct(s) && isfield(s, field)
    v = s.(field);
else
    v = default;
end
end
