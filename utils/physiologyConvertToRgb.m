function rgb = physiologyConvertToRgb(scalarField, intensityImage, maskIn, mn, mx, whiteUse, cMap, bipolar, intensityMin, intensityMax)
%PHYSIOLOGYCONVERTTORGB  Colour-encode a scalar physiology field using HSV.
%
%   rgb = physiologyConvertToRgb(scalarField, intensityImage, maskIn, ...
%                                mn, mx, whiteUse, cMap, bipolar, ...
%                                intensityMin, intensityMax)
%
% Maps a per-pixel scalar quantity (e.g. dF/F, ratio, CV, normalised recovery)
% to a colour-coded RGB image.  The colour (hue) encodes the scalar field;
% brightness (value) is modulated by the underlying intensity image so that
% the ER network structure remains visible within the colour map.
%
% ENCODING SCHEME
%   Hue   : scalar field rescaled to [mn, mx] → index into cMap.
%   Saturation and Value: depend on whiteUse flag:
%     whiteUse = true, mn >= 0 (unipolar, e.g. dF/F, CV, recovery):
%       S = sqrt(V_image)   (full saturation at bright pixels, white at dim)
%       V = 1               (always bright)
%     whiteUse = true, mn < 0 (bipolar, e.g. ΔF/F with negative dips):
%       S = 1               (always saturated)
%       V = sqrt(V_image) × 0.625 + 0.375   (dark background, not black)
%     whiteUse = false:
%       S = 1               (always saturated)
%       V = sqrt(V_image)   (pure darkness at zero intensity)
%
% INPUTS
%   scalarField    – [nY nX nC nZ nT] single/double; quantity to colour-encode.
%                   NaN values are rendered as background (black or white
%                   depending on whiteUse).
%   intensityImage – [nY nX nC nZ nT] numeric; background intensity for
%                   brightness modulation.  Must match scalarField size.
%   maskIn         – [nY nX ...] logical; pixels outside mask rendered as
%                   background.  May have fewer C/Z/T dimensions (broadcast).
%   mn, mx         – scalar; range of scalarField mapped to full colourmap.
%                   Values outside [mn, mx] are clamped.
%   whiteUse       – logical; true = white background encoding (see above).
%   cMap           – [N×3] double colourmap matrix.
%   bipolar        – (optional) logical; selects which whiteUse=true sub-mode
%                   applies (see above). Defaults to (mn < 0) for backward
%                   compatibility. Pass explicitly to decouple the choice from
%                   the sign of mn -- e.g. Eh (mV) is always negative-valued
%                   but isn't a true bipolar quantity, so callers rendering it
%                   should pass bipolar=false to get the white-background mode.
%   intensityMin, intensityMax – (optional) scalar in [0,1]; explicit
%                   low/high input range the brightness (Value) channel is
%                   stretched against (via imadjust) after mat2gray's own
%                   per-plane [0,1] normalisation. Defaults to [0 1] (no
%                   extra stretch) if omitted -- pixels at or below
%                   intensityMin render black, at or above intensityMax
%                   render full brightness. Previously this stage used
%                   imadjust's own single-argument auto-stretch
%                   (effectively ignoring any caller-specified range,
%                   which is why the RatioRgb panel's Intensity Min/Max
%                   fields never had any visible effect).
%
% OUTPUT
%   rgb – [nY nX nC nZ nT 3] uint8; colour-coded image.

if nargin < 8 || isempty(bipolar)
    bipolar = mn < 0;
end
if nargin < 9 || isempty(intensityMin)
    intensityMin = 0;
end
if nargin < 10 || isempty(intensityMax)
    intensityMax = 1;
end
if intensityMin >= intensityMax
    intensityMin = 0;
    intensityMax = 1;
end

[nY, nX, nC, nZ, nT] = size(scalarField);
[~,  ~,  mC, mZ, mT] = size(maskIn);

rgb = zeros(nY, nX, nC, nZ, nT, 3, 'uint8');

for iT = 1:nT
    for iZ = 1:nZ
        for iC = 1:nC
            sf_  = scalarField(:,:,iC,iZ,iT);
            img  = intensityImage(:,:,iC,iZ,iT);
            msk  = maskIn(:,:, min(iC,mC), min(iZ,mZ), min(iT,mT));

            % NaN pixels treated as background (outside mask)
            msk  = msk & ~isnan(sf_);

            % ---- hue from scalar field ----------------------------------
            im = rescale(double(sf_), 0, 1, 'InputMin', mn, 'InputMax', mx);
            im(~msk) = 0;

            hsvIm = rgb2hsv(ind2rgb(uint8(255 .* im), cMap));

            % ---- value channel from intensity image ---------------------
            V = imadjust(im2single(mat2gray(double(img))), [intensityMin intensityMax]);
            V(~msk) = 0;

            % ---- compose HSV --------------------------------------------
            if whiteUse
                if bipolar
                    % bipolar (e.g. ΔF/F with baseline negative): dark bg
                    hsvIm(:,:,3) = sqrt(V) * 0.625 + 0.375;
                else
                    % unipolar (e.g. dF/F ≥ 0, CV, recovery): white bg
                    hsvIm(:,:,2) = sqrt(V);
                    hsvIm(:,:,3) = ones(nY, nX, 'single');
                end
            else
                hsvIm(:,:,3) = sqrt(V);
            end

            hsvIm(~msk) = 0;
            rgb(:,:,iC,iZ,iT,1:3) = im2uint8(hsv2rgb(hsvIm));
        end
    end
end

end
