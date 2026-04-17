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
%   ehRgb    – reduction potential E_hh (mV);           diverging colormap
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
%                     p.ratioMin   – lower clip for ratio display (default 0.1)
%                     p.ratioMax   – upper clip for ratio display (default 1.0)
%                     p.ehMin      – lower clip for Eh display (mV; default −400)
%                     p.ehMax      – upper clip for Eh display (mV; default −200)
%                     p.whiteUse   – logical; white-background encoding
%                                    (default true)
%
% OUTPUTS
%   ratioRgb – [nY nX nC nZ nT 3] uint8; colour-coded ratio image.
%   oxdRgb   – [nY nX nC nZ nT 3] uint8; colour-coded OxD image, or [].
%   ehRgb    – [nY nX nC nZ nT 3] uint8; colour-coded Eh image, or [].
%
% COLOURMAPS
%   Ratio and OxD: colorcet('R2') if available; else a red-yellow-white
%                  sequential ramp.
%   Eh (diverging): colorcet('D7') if available; else a red-white-blue ramp.
%   For other probe types, the caller can replace intensityImage/maskIn but
%   the colourmap assignment matches the standard grx1-roGFP2 convention.
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
if exist('colorcet', 'file')
    cMapSeq  = colorcet('R2');   % ratio and OxD: perceptually uniform seq.
    cMapDiv  = colorcet('D7');   % Eh: diverging blue-white-red
else
    cMapSeq  = sequentialColormap(256);
    cMapDiv  = divergingColormap(256);
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

% ---- Eh RGB (Eh in mV; bipolar → diverging colormap, dark background) ---
if ~isempty(ehIn)
    ehRgb = physiologyConvertToRgb(ehIn, intensityImage, maskIn, ...
        ehMin, ehMax, false, cMapDiv);
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
function cmap = divergingColormap(n)
%DIVERGINGCOLORMAP  Blue–white–red diverging colormap with n levels.
% Centre (white) maps to the midpoint value; same as opticalFlowRgb fallback.
half  = floor(n / 2);
red   = [linspace(0.7, 1, half), linspace(1, 1, half)]';
green = [linspace(0.0, 1, half), linspace(1, 0, half)]';
blue  = [linspace(0.0, 1, half), linspace(1, 0.7, half)]';
cmap  = flipud([red, green, blue]);
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
