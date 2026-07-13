function [ratioRgb, oxdRgb, ehRgb] = ...
        ratioRgb(intensityImage, maskIn, ratioIn, oxdIn, ehIn, p)
%RATIORGB  Produce colour-coded RGB images from ratiometric physiology data.
%
%   [ratioRgb, oxdRgb, ehRgb] = ...
%       ratioRgb(intensityImage, maskIn, ratioIn, oxdIn, ehIn, p)
%
% Creates up to three colour-coded visualisations from ratiometric probe data:
%
%   ratioRgb – raw channel ratio R = Ch_num / Ch_den;   sequential colormap
%   oxdRgb   – oxidised fraction OxD ∈ [0,1];           sequential colormap
%   ehRgb    – reduction potential E_hh (mV);           sequential colormap
%              (configurable independently -- see p.ehColormap)
%
% oxdRgb and ehRgb are empty ([]) when oxdIn or ehIn are empty.
%
% INPUTS
%   intensityImage – [nY nX nC nZ nT] numeric; background intensity for
%                   brightness modulation.
%   maskIn         – [nY nX ...] logical; ER/cell mask.  May have fewer
%                   C/Z/T dimensions (broadcast).
%   ratioIn        – [nY nX nC nZ nT] single; ratio from ratioCompute.
%   oxdIn          – [nY nX nC nZ nT] single; OxD from redoxCalibrate,
%                   or [] to skip oxdRgb.
%   ehIn           – [nY nX nC nZ nT] single; Eh (mV) from redoxCalibrate,
%                   or [] to skip ehRgb.
%   p              – parameter struct:
%                     p.ratioMin      – lower clip for ratio display (default 0.1)
%                     p.ratioMax      – upper clip for ratio display (default 1.0)
%                     p.ehMin         – lower clip for Eh display (mV; default −400)
%                     p.ehMax         – upper clip for Eh display (mV; default −200)
%                     p.whiteUse      – logical; white-background encoding
%                                       (default true)
%                     p.ratioColormap – colorcet name for ratio/OxD (default 'R2')
%                     p.ehColormap    – colorcet name for Eh (default 'R2')
%
% OUTPUTS
%   ratioRgb – [nY nX nC nZ nT 3] uint8; colour-coded ratio image.
%   oxdRgb   – [nY nX nC nZ nT 3] uint8; colour-coded OxD image, or [].
%   ehRgb    – [nY nX nC nZ nT 3] uint8; colour-coded Eh image, or [].
%
% COLOURMAPS
%   Default is colorcet('R2') for all three (perceptually-uniform
%   sequential/rainbow) if colorcet is on the path; caller can override via
%   p.ratioColormap/p.ehColormap with any valid colorcet name. Falls back to
%   a built-in sequential ramp if colorcet isn't available (name overrides
%   are then ignored, since there's no catalog to look them up in).
%
% DEPENDENCIES
%   physiologyConvertToRgb  (Physiology_sandbox/utils/)
%   colorcet  (optional; perceptually uniform colourmaps, Peter Kovesi).
%             Falls back to built-in colormaps if not on the MATLAB path.

ratioMin = sf(p, 'ratioMin',  0.1);
ratioMax = sf(p, 'ratioMax',  1.0);
ehMin    = sf(p, 'ehMin',   -400);
ehMax    = sf(p, 'ehMax',   -200);
whiteUse = logical(sf(p, 'whiteUse', true));

% ---- colourmaps ---------------------------------------------------------
% A diverging map (e.g. 'D7') visually centres its neutral/white point at
% the MIDPOINT of [mn,mx], not at Eh=0 -- physiologyConvertToRgb has no
% zero-crossing concept, it's a plain linear map. For an all-negative Eh
% range that wastes most of the palette on one side and can read as a flat
% single hue. Default is now the same sequential map used for ratio/OxD
% ('R2') unless the caller asks for something else.
if exist('colorcet', 'file')
    cMapSeq = colorcet(sf(p, 'ratioColormap', 'R2'));
    cMapDiv = colorcet(sf(p, 'ehColormap',    'R2'));
else
    cMapSeq = sequentialColormap(256);
    cMapDiv = sequentialColormap(256);
end

% ---- ratio RGB ----------------------------------------------------------
ratioRgb = physiologyConvertToRgb(ratioIn, intensityImage, maskIn, ...
    ratioMin, ratioMax, whiteUse, cMapSeq);

% ---- OxD RGB (OxD ∈ [0,1] always unipolar) -----------------------------
if ~isempty(oxdIn)
    oxdRgb = physiologyConvertToRgb(oxdIn, intensityImage, maskIn, ...
        0, 1, whiteUse, cMapSeq);
else
    oxdRgb = [];
end

% ---- Eh RGB (Eh in mV) ---------------------------------------------------
% Eh is always negative-valued (mV) but isn't a true bipolar quantity (no
% meaningful zero-crossing within its measured range) -- bipolar=false so
% whiteUse's white-background mode isn't overridden by the sign of ehMin.
if ~isempty(ehIn)
    ehRgb = physiologyConvertToRgb(ehIn, intensityImage, maskIn, ...
        ehMin, ehMax, whiteUse, cMapDiv, false);
else
    ehRgb = [];
end

end % ratioRgb


% =========================================================================
function cmap = sequentialColormap(n)
%SEQUENTIALCOLORMAP  Red–yellow–white sequential colormap with n levels.
% Used as fallback when colorcet is not available.
r = ones(n, 1);
g = linspace(0, 1, n)';
b = linspace(0, 0.8, n)';
cmap = [r, g, b];
end

% =========================================================================
function v = sf(s, field, default)
%SF  Safe field access with default.
if isstruct(s) && isfield(s, field)
    v = s.(field);
else
    v = default;
end
end
