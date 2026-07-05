function [dal,raw,cvd] = daltonize(rgb,typ,sev)
% Daltonize/recolor an image to improve contrast for CVD observers (protanomaly, deuteranomaly, tritanomaly)
%
% (c) 2026 Stephen Cobeldick
%
% Adjusts the colors of an image so that an observer with color vision
% deficiency (CVD) can better distinguish details that would otherwise be
% lost. This is the counterpart of CVDSIM: CVDSIM shows a normal-sighted
% viewer what a CVD observer sees; DALTONIZE instead modifies a (normal)
% image so that a CVD observer sees more of the original information.
%
%%% Syntax %%%
%
%   dal = daltonize(rgb,typ)
%   dal = daltonize(rgb,typ,sev)
%   [dal,raw,cvd] = daltonize(...)
%
%% Algorithm %%
%
% Method (Fidaner, Lin & Ozguven, 2005, "Analysis of Color Blindness"):
%  1. Simulate the image as seen by the CVD observer using CVDSIM().
%  2. Compute the error: the difference between the original and simulated
%     colors, i.e. the information the CVD observer's cones cannot convey.
%  3. Redistribute that error into the color channels the deficiency does
%     not compromise, using a fixed, type-specific 3x3 matrix.
%  4. Add the redistributed error back onto the original image.
%
% Unlike a common (and questionable) shortcut, steps 1-4 here are all
% performed on linearized RGB, with a single sRGB decode at the start and
% a single re-encode at the end (matching CVDSIM's own linear-light
% pipeline): doing the error/redistribution arithmetic directly on
% gamma-encoded values is not physically meaningful, since color mixing
% is additive in linear light, not in gamma-encoded space.
%
% The error-redistribution matrix is type-specific: many public
% implementations reuse the protanopia matrix for all three types,
% which Simon-Liedtke and Farup explains is incorrect.
%
% References:
%  Fidaner O, Lin P, Ozguven N: "Analysis of Color Blindness", 2005.
%  <http://scien.stanford.edu/pages/labsite/2005/psych221/projects/05/ofidaner/colorblindness_project.htm>
%  Type-specific matrices and their derivation:
%  <http://ixora.io/projects/colorblindness/color-blindness-simulation-research.html>
%  Simon-Liedtke J T, Farup I: "Evaluating color vision deficiency
%  daltonization methods using a behavioral visual-search method", J. Vis.
%  Commun. Image Represent., 2016. <https://ivarfa.folk.ntnu.no/publications/journal/Simon_16_jvci.pdf>
%
%% Examples %%
%
%%% Enhance a single RGB triple for protanopia %%%
%
%   >> daltonize([1,0,0], 'protan')
%   ans = [1.0000  0.7213  0.7961]
%
%%% Daltonize an image for moderate deuteranopy %%%
%
%   >> rgb = imread("peppers.png");
%   >> imshow(daltonize(rgb, 'deutan',0.6))
%
%%% View Parula enhanced for all three CVD types %%%
%
%   >> I = permute(parula(17),[3,1,2]);
%   >> imshow([I;daltonize(I,'p');daltonize(I,'d');daltonize(I,'t')])
%
%% Notes %%
%
% * The 0.7 coefficients in the redistribution matrices are heuristic
%   constants from the original Fidaner/Lin/Ozguven algorithm, not derived
%   from first-principles color science; using type-specific matrices and
%   a single, consistent (linear-light) color space is more internally
%   consistent than the common alternative, but the underlying 0.7 weights
%   remain an empirical choice, not a rigorously optimal one.
% * This is a fixed, content-independent transform: it can overcorrect
%   (introducing an unwanted color cast) for images that did not need much
%   correction, and does not adapt to the image's actual color distribution
%   the way newer naturalness-preserving methods do, e.g.:
%   Kuhn, Oliveira & Fernandes (2008), IEEE TVCG 14(6):1747-1754;
%   Machado & Oliveira (2010), Computer Graphics Forum 29(3):933-942;
%   Milic et al. (2015), J. Imaging Sci. Technol. 59(1):10504.
%
%% Input Arguments (**=default) %%
%
%   rgb = Numeric array of sRGB values to convert. Floating point values
%         must be 0<=rgb<=1, integer must be 0<=rgb<=intmax(class(rgb)).
%         Size Nx3 or RxCx3, the last dimension encodes the R,G,B values.
%   typ = CharRowVector or StringScalar, the type of CVD to correct for:
%         'p' / 'protan' / 'protanomaly'   / 'protanopia'   (L-cone deficiency).
%         'd' / 'deutan' / 'deuteranomaly' / 'deuteranopia' (M-cone deficiency).
%         't' / 'tritan' / 'tritanomaly'   / 'tritanopia'   (S-cone deficiency).
%   sev = NumericScalar, 0<=sev<=1, the severity of the deficiency to
%         correct for. Passed through to CVDSIM: see CVDSIM for details.
%         1** = complete dichromacy (protanopia/deuteranopia) or the most
%               severe tabulated case of tritanomaly.
%
%% Output Arguments %%
%
%   dal = NumericArray, the same size and class as <rgb>, the daltonized
%         (recolored) image. Float values are clipped to 0..1.
%   raw = FloatArray, the same size as <rgb>, the daltonized image without
%         clipping (i.e. values may be outside 0..1).
%   cvd = NumericArray, the same size and class as <rgb>, containing the
%         simulated CVD colors. Float values are clipped to 0..1.
%
%% Dependencies %%
%
% * MATLAB R2009b or later.
% * cvdsim().
%
% See also CVDSIM MACHADO2010
% COLORMAP COLORORDER BREWERMAP MAXDISTCOLOR
% SRGB_TO_CAM02UCS CAM02UCS_TO_SRGB SRGB_TO_CAM16UCS CAM16UCS_TO_SRGB

% Release | Feature
% --------|--------
% R2016b  | string class                                 [only if supplied]
% R2009b  | tilde argument placeholder
% R2008a  | assert: message-identifier
%
%% Input Wrangling %%
%
% CVDSIM performs input checking, LIN & SIM are Nx3:
if nargin<3
	[cvd,~,lin,sim] = cvdsim(rgb,typ);
else
	[cvd,~,lin,sim] = cvdsim(rgb,typ,sev);
end
%
isz = size(rgb);
icl = class(rgb);
if isfloat(rgb)
	mxv = 1;
elseif isinteger(rgb)
	mxv = double(intmax(icl));
end
%
% Matrices from Ixora.io's ColorBlindness library (see References).
% Each type's own row is zero (no self-correction for the lost channel);
% its error is split 0.7/0.7 into the two channels the deficiency spares.
%
switch lower(typ)
	case {'p','protan','protanopia','protanomaly'}
		e2m = [0.0,0.0,0.0; 0.7,1.0,0.0; 0.7,0.0,1.0];
	case {'d','deutan','deuteranopia','deuteranomaly'}
		e2m = [1.0,0.7,0.0; 0.0,0.0,0.0; 0.0,0.7,1.0];
	case {'t','tritan','tritanopia','tritanomaly'}
		e2m = [1.0,0.0,0.7; 0.0,1.0,0.7; 0.0,0.0,0.0];
	otherwise
		error('SC:daltonize:typ:NotSupported',...
			'Second input <typ> "%s" is not supported: use "protan"/"deutan"/"tritan" or their full names or their initials.',typ)
end
%
%% Simulate, Compute Error, Redistribute, and Recombine (in linear light) %%
%
err = lin - sim;
lin = lin + err*e2m.';
raw = reshape(sGammaCor(lin),isz);
%
if mxv>1
	dal = cast(mxv*raw,icl);
else
	dal = min(1,max(0,raw));
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%daltonize
function out = sGammaCor(inp)
% Gamma correction: Nx3 linear RGB -> Nx3 sRGB.
idx = inp > 0.0031308;
out = 12.92 * inp;
out(idx) = real(1.055 * inp(idx) .^ (1./2.4) - 0.055);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sGammaCor
%
% Code and Implementation:
% Copyright (c) 2026 Stephen Cobeldick
% Error-redistribution Matrix Only:
% Copyright (c) 2005 Onur Fidaner, Poliang Lin, and Nevran Ozguven.
% Type-specific Matrix Values Only:
% Source: Ixora.io ColorBlindness library research page (see References above).
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%license