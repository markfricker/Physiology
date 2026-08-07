function [oxdSeries, ehSeries] = redoxCalibrate(ratioIn, p)
%REDOXCALIBRATE  Convert ratiometric probe data to OxD and reduction potential.
%
%   [oxdSeries, ehSeries] = redoxCalibrate(ratioIn, p)
%
% Two-step calibration for genetically encoded ratiometric redox sensors
% (e.g. grx1-roGFP2, roGFP1, Hyper).  Operates on ratio images from
% ratioCompute; no channel selection is performed here.
%
% STEP 1 — Normalisation to oxidised fraction (OxD)
%   OxD = (R − Rmin) / [delta × (Rmax − R) + (R − Rmin)]
%
%   where delta (iFactor) is the instrument factor:
%     delta = (F_ex1_ox / F_ex2_ox) / (F_ex1_red / F_ex2_red)
%   accounting for the fluorescence intensities at the fully-reduced and
%   fully-oxidised probe endpoints.  OxD = 0 → fully reduced; OxD = 1 →
%   fully oxidised.
%
% STEP 2 — Nernst equation to reduction potential (mV)
%   E_hh = E0' − (RT/nF) × ln[(1 − OxD) / OxD]
%
%   The midpoint potential E0' is pH-corrected (one-electron term) from the
%   standard midpoint measured at the reference pH:
%     E0' = midpoint − (1000 × 2.303 × R × T / F) × (pH_sample − pH_ref)
%
%   This follows the standard roGFP calibration (Hanson et al. 2004;
%   Gutscher et al. 2008; Morgan et al. 2011).
%
% INPUTS
%   ratioIn – [nY nX nC nZ nT] single; ratio image from ratioCompute.
%   p       – parameter struct:
%               p.rMin     – scalar; fully-reduced endpoint ratio (R_red).
%               p.rMax     – scalar; fully-oxidised endpoint ratio (R_ox).
%               p.iFactor  – scalar ≥ 0; instrument factor (delta).  Set to
%                            1 (default) if endpoints are measured under
%                            identical acquisition conditions.
%               p.midpoint – scalar; probe E°' in mV at the reference pH
%                            (e.g. −280 mV for roGFP2; −300 mV for
%                            grx1-roGFP2).  Default −300.
%               p.temp     – scalar; temperature in K.  Default 293 (20 °C).
%               p.pH       – scalar; sample pH.  Default 7.4.
%               p.pH0      – scalar; reference pH at which midpoint was
%                            measured.  Default 7.4 (no pH correction).
%               p.n        – scalar; number of electrons (z = 2 for roGFP
%                            disulphide).  Default 2.
%
% OUTPUTS
%   oxdSeries – [nY nX nC nZ nT] single; per-pixel oxidised fraction [0,1].
%               Clipped to [0.001, 0.999] to avoid Nernst singularities at
%               the endpoints.  Pixels outside [rMin, rMax] are clipped before
%               computation.
%   ehSeries  – [nY nX nC nZ nT] single; per-pixel reduction potential (mV).
%               More negative = more reducing environment.
%
% TYPICAL PARAMETER VALUES (grx1-roGFP2 in mammalian cells)
%   rMin     = 0.1    (fully-reduced endpoint ratio, 10 mM DTT)
%   rMax     = 1.0    (fully-oxidised endpoint ratio, 200 µM aldrithiol)
%   iFactor  = 5.5    (instrument-dependent; measure from calibration images)
%   midpoint = −280   (E°' of roGFP2 vs. SHE at pH 7.0; −300 for grx1 couple)
%   temp     = 310    (37 °C for mammalian culture)
%   pH       = 7.0    (ER lumen typical value)
%   pH0      = 7.0    (reference pH used during midpoint measurement)
%
% NOTE ON pH CORRECTION
%   The one-electron pH correction (2.303 RT/F per pH unit ≈ 59 mV/pH at
%   25 °C) accounts for the proton-coupled nature of the thiol/disulphide
%   equilibrium.  If pH_sample = pH0, no correction is applied.

rMin   = sf(p, 'rMin',      0.1);
rMax   = sf(p, 'rMax',      1.0);
delta  = sf(p, 'iFactor',   1.0);
Emid   = sf(p, 'midpoint', -300);
T      = sf(p, 'temp',      293);
pH     = sf(p, 'pH',        7.4);
pH0    = sf(p, 'pH0',       7.4);
n      = sf(p, 'n',           2);

R_gas  = 8.314;    % J mol⁻¹ K⁻¹
F_fara = 96485;    % C mol⁻¹

% ---- pH-corrected midpoint potential ------------------------------------
% One-electron (z=1) pH correction, base-10 log via 2.303 factor
E0prime = Emid - (1000 * (2.303 * R_gas * T) / (1 * F_fara)) * (pH - pH0);

ratioIn = im2single(ratioIn);

% ---- clamp ratio to calibration range -----------------------------------
R = ratioIn;
R(R < rMin) = rMin;
R(R > rMax) = rMax;

% ---- step 1: OxD normalisation ------------------------------------------
num_oxd = R - rMin;
den_oxd = delta .* (rMax - R) + num_oxd;
den_oxd(den_oxd <= 0) = eps('single');   % guard: only at rMin=rMax edge case

oxdSeries = num_oxd ./ den_oxd;
% MATLAB's two-argument min/max silently DROP a NaN operand instead of
% propagating it (e.g. min(0.999, NaN) returns 0.999, not NaN) -- without
% guarding this, a NaN-padded input (e.g. ratioIn NaN for a channel ratio
% was never computed for) would get silently replaced by a fabricated
% 0.999 "oxidised" value here, indistinguishable from real data downstream.
nanMask = isnan(oxdSeries);
oxdSeries = max(0.001, min(0.999, oxdSeries));   % avoid Nernst singularities
oxdSeries(nanMask) = NaN;

% ---- step 2: Nernst equation to mV --------------------------------------
% E_hh = E0' − (RT/nF) × ln[(1−OxD)/OxD]   (natural log, result in mV)
ehSeries = E0prime - ...
    (1000 * (R_gas * T) / (n * F_fara)) .* log((1 - oxdSeries) ./ oxdSeries);

end % redoxCalibrate


% =========================================================================
function v = sf(s, field, default)
%SF  Safe field access with default.
if isstruct(s) && isfield(s, field)
    v = s.(field);
else
    v = default;
end
end
