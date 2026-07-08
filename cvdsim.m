function [cvd,raw,lin,sim] = cvdsim(rgb,typ,sev,gamma)
% Simulate dichromat color vision deficiency (CVD) for a trichromat observer.
%
% (c) 2026 Stephen Cobeldick
%
% Simulates the perceived colors seen by an observer with color vision
% deficiency (CVD), using the physiologically-based model of Machado,
% Oliveira & Fernandes (2009). Anomalous trichromacy and dichromacy are
% handled in a single unified way, controlled by a continuously-variable
% severity parameter <sev> (0 = normal vision, 1 = complete dichromacy).
%
%%% Syntax %%%
%
%   cvd = cvdsim(rgb,typ)
%   cvd = cvdsim(rgb,typ,sev)
%   cvd = cvdsim(rgb,typ,sev,gamma)
%   [cvd,raw] = cvdsim(...)
%
%% Algorithm %%
%
% Reference:
%  Machado G M, Oliveira M M, Fernandes L A F: "A Physiologically-based Model
%  for Simulation of Color Vision Deficiency" IEEE TVCG 15(6):1291-1298, 2009.
%  <https://doi.org/10.1109/TVCG.2009.113>
%  Author-hosted matrices: <https://www.inf.ufrgs.br/~oliveira/pubs_files/CVD_Simulation/CVD_Simulation.html>
%
% The simulated color is obtained by multiplying the linearized input sRGB
% values with a 3x3 matrix <mat>, selected according to <typ> and <sev>.
% Marix <mat> is obtained by linearly interpolating between the two nearest
% of the eleven matrices published by the paper's authors (tabulated at
% severity steps of 0.1); this is the fast approximation recommended by
% the authors themselves for intermediate severities.
%
%% Examples %%
%
%%% Simulate a single RGB triple for protanopia %%%
%
%   >> cvdsim([1,0,0], 'protan')
%   ans = [0.4266  0.3727  0]
%
%%% View an image as seen with moderate deuteranomaly %%%
%
%   >> I = imread("peppers.png");
%   >> imshow(cvdsim(I, 'deutan', 0.6))
%
%%% View Parula colormap under all three CVD types %%%
%
%   >> I = permute(parula(17),[3,1,2]);
%   >> imshow([I;cvdsim(I,'p');cvdsim(I,'d');cvdsim(I,'t')])
%
%% Input Arguments (**=default) %%
%
%   rgb = Numeric array of sRGB values to convert. Floating point values
%         must be 0<=rgb<=1, integer must be 0<=rgb<=intmax(class(rgb)).
%         Size Nx3 or RxCx3, the last dimension encodes the R,G,B values.
%   typ = CharRowVector or StringScalar, the type of CVD to simulate:
%         'p' / 'protan' / 'protanomaly'   / 'protanopia'   (L-cone deficiency).
%         'd' / 'deutan' / 'deuteranomaly' / 'deuteranopia' (M-cone deficiency).
%         't' / 'tritan' / 'tritanomaly'   / 'tritanopia'   (S-cone deficiency):
%         True tritanopia (i.e. complete S-cone loss) is NOT modeled: the
%         matrices approximate tritanomaly only, via a spectral-shift
%         method, which is informally extended up to severity 1.
%   sev = NumericScalar, the severity of the deficiency, with value range:
%         0   = normal color vision (<gam> is the identity matrix).
%         1** = complete dichromacy (protanopia/deuteranopia)
%               or the most severe tabulated case of tritanomaly.
%   gamma = LogicalScalar controlling whether the simulation matrices
%         are applied to linearized RGB values, where:
%         true** = apply inverse sRGB gamma correction, simulate in
%                  linear RGB, then reapply sRGB gamma correction.
%                  This is the theoretically correct interpretation.
%         false  = apply the simulation matrices directly to sRGB values
%                  values without gamma correction. This more closely
%                  reproduces the images embedded in Machado (2009) PDF.
%
%% Output Arguments %%
%
%   cvd = NumericArray, the same size and class as <rgb>, containing the
%         simulated CVD colors. Float values are clipped to 0..1.
%   raw = FloatArray, the same size as <rgb>, containing the simulated
%         CVD colors without clipping (i.e. values may be outside 0..1).
%
%% Dependencies %%
%
% * MATLAB R2008a or later.
%
% See also BRETTEL1997 DALTONIZER MACHADO2010 MILIC2015
% IMSHOW PARULA LINES COLORMAP COLORORDER BREWERMAP MAXDISTCOLOR

% Release | Feature
% --------|--------
% R2016b  | string class                                 [only if supplied]
% R2008a  | assert: message-identifier
%
%% Input Wrangling %%
%
isz = size(rgb);
icl = class(rgb);
if isfloat(rgb)
	mxv = 1;
elseif isinteger(rgb)
	mxv = double(intmax(icl));
	rgb = double(rgb)./mxv;
else
	error('SC:cvdsim:rgb:NotNumeric',...
		'1st input <rgb> must be a numeric array, not %s',class(rgb))
end
assert(isreal(rgb),...
	'SC:cvdsim:rgb:NotReal',...
	'1st input <rgb> must be a real array (not complex).')
assert(isz(end)==3 || isequal(isz,[3,1]),...
	'SC:cvdsim:rgb:InvalidSize',...
	'1st input <rgb> last dimension must have size 3 (e.g. Nx3 or RxCx3).')
assert(all(0<=rgb(:)&rgb(:)<=1),'SC:cvdsim:rgb:OutOfRange',...
	'1st input <rgb> values must be 0<=rgb<=%d',mxv)
%
typ = svcSS2C(typ);
assert(ischar(typ)&&ndims(typ)==2&&size(typ,1)==1,...
	'SC:cvdsim:typ:NotText',...
	'Second input <typ> must be a character vector or a string scalar.') %#ok<ISMAT>
%
switch lower(typ)
	case {'p','protan','protanopia','protanomaly'}
		arr = cat(3, eye(3), [+0.856167,+0.182038,-0.038205;+0.029342,+0.955115,+0.015544;-0.002880,-0.001563,+1.004443], [+0.734766,+0.334872,-0.069637;+0.051840,+0.919198,+0.028963;-0.004928,-0.004209,+1.009137], [+0.630323,+0.465641,-0.095964;+0.069181,+0.890046,+0.040773;-0.006308,-0.007724,+1.014032], [+0.539009,+0.579343,-0.118352;+0.082546,+0.866121,+0.051332;-0.007136,-0.011959,+1.019095], [+0.458064,+0.679578,-0.137642;+0.092785,+0.846313,+0.060902;-0.007494,-0.016807,+1.024301], [+0.385450,+0.769005,-0.154455;+0.100526,+0.829802,+0.069673;-0.007442,-0.022190,+1.029632], [+0.319627,+0.849633,-0.169261;+0.106241,+0.815969,+0.077790;-0.007025,-0.028051,+1.035076], [+0.259411,+0.923008,-0.182420;+0.110296,+0.804340,+0.085364;-0.006276,-0.034346,+1.040622], [+0.203876,+0.990338,-0.194214;+0.112975,+0.794542,+0.092483;-0.005222,-0.041043,+1.046265], [+0.152286,+1.052583,-0.204868;+0.114503,+0.786281,+0.099216;-0.003882,-0.048116,+1.051998]);
	case {'d','deutan','deuteranopia','deuteranomaly'}
		arr = cat(3, eye(3), [+0.866435,+0.177704,-0.044139;+0.049567,+0.939063,+0.011370;-0.003453,+0.007233,+0.996220], [+0.760729,+0.319078,-0.079807;+0.090568,+0.889315,+0.020117;-0.006027,+0.013325,+0.992702], [+0.675425,+0.433850,-0.109275;+0.125303,+0.847755,+0.026942;-0.007950,+0.018572,+0.989378], [+0.605511,+0.528560,-0.134071;+0.155318,+0.812366,+0.032316;-0.009376,+0.023176,+0.986200], [+0.547494,+0.607765,-0.155259;+0.181692,+0.781742,+0.036566;-0.010410,+0.027275,+0.983136], [+0.498864,+0.674741,-0.173604;+0.205199,+0.754872,+0.039929;-0.011131,+0.030969,+0.980162], [+0.457771,+0.731899,-0.189670;+0.226409,+0.731012,+0.042579;-0.011595,+0.034333,+0.977261], [+0.422823,+0.781057,-0.203881;+0.245752,+0.709602,+0.044646;-0.011843,+0.037423,+0.974421], [+0.392952,+0.823610,-0.216562;+0.263559,+0.690210,+0.046232;-0.011910,+0.040281,+0.971630], [+0.367322,+0.860646,-0.227968;+0.280085,+0.672501,+0.047413;-0.011820,+0.042940,+0.968881]);
	case {'t','tritan','tritanopia','tritanomaly'}
		arr = cat(3, eye(3), [+0.926670,+0.092514,-0.019184;+0.021191,+0.964503,+0.014306;+0.008437,+0.054813,+0.936750], [+0.895720,+0.133330,-0.029050;+0.029997,+0.945400,+0.024603;+0.013027,+0.104707,+0.882266], [+0.905871,+0.127791,-0.033662;+0.026856,+0.941251,+0.031893;+0.013410,+0.148296,+0.838294], [+0.948035,+0.089490,-0.037526;+0.014364,+0.946792,+0.038844;+0.010853,+0.193991,+0.795156], [+1.017277,+0.027029,-0.044306;-0.006113,+0.958479,+0.047634;+0.006379,+0.248708,+0.744913], [+1.104996,-0.046633,-0.058363;-0.032137,+0.971635,+0.060503;+0.001336,+0.317922,+0.680742], [+1.193214,-0.109812,-0.083402;-0.058496,+0.979410,+0.079086;-0.002346,+0.403492,+0.598854], [+1.257728,-0.139648,-0.118081;-0.078003,+0.975409,+0.102594;-0.003316,+0.501214,+0.502102], [+1.278864,-0.125333,-0.153531;-0.084748,+0.957674,+0.127074;-0.000989,+0.601151,+0.399838], [+1.255528,-0.076749,-0.178779;-0.078411,+0.930809,+0.147602;+0.004733,+0.691367,+0.303900]);
	otherwise
		error('SC:cvdsim:typ:NotSupported',...
			'Second input <typ> "%s" is not supported: use "protan"/"deutan"/"tritan" or their full names or their initials.',typ)
end
%
if nargin<3
	sev = 10;
else
	assert(isnumeric(sev)&&isscalar(sev)&&isreal(sev),...
		'SC:cvdsim:sev:NotUnitScalar',...
		'Third input <sev> must be a real numeric scalar.')
	assert(sev>=0&&sev<=1,...
		'SC:cvdsim:sev:OutOfRange',...
		'Third input <sev> must have a value between 0 and 1 (inclusive).')
	sev = 10*double(sev);
end
%
if nargin<4
	gamma = true;
else
	assert(isequal(gamma,false)||isequal(gamma,true),...
		'SC:cvdsim:gamma:NotScalarLogical',...
		'Fourth input <sev> must be true/1 or false/0.')
	gamma = logical(gamma);
end
%
%% Interpolate the Simulation Matrix %%
%
id0 = min(floor(sev),9);
adj = sev-id0;
mat = (1-adj).*arr(:,:,id0+1) + adj.*arr(:,:,id0+2);
%
%% Apply the Transformation %%
%
if gamma
	lin = reshape(sGammaInv(rgb),[],3);
else
	lin = reshape(rgb,[],3);
end
%
sim = lin*mat.';
%
if gamma
	raw = reshape(sGammaCor(sim),isz);
else
	raw = reshape(sim,isz);
end
%
if mxv>1
	cvd = cast(mxv*raw,icl);
else
	cvd = min(1,max(0,raw));
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%cvdsim
function out = sGammaCor(inp)
% Gamma correction: Nx3 linear RGB -> Nx3 sRGB.
idx = inp > 0.0031308;
out = 12.92 * inp;
out(idx) = real(1.055 * inp(idx) .^ (1./2.4) - 0.055);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sGammaCor
function out = sGammaInv(inp)
% Inverse gamma correction: Nx3 sRGB -> Nx3 linear RGB.
idx = inp > 0.04045;
out = inp / 12.92;
out(idx) = real(((inp(idx) + 0.055) ./ 1.055) .^ 2.4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sGammaInv
function arr = svcSS2C(arr)
% If scalar string then extract the character vector, otherwise data is unchanged.
if isa(arr,'string') && isscalar(arr)
	arr = arr{1};
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%svcSS2C
%
% Code and Implementation:
% Copyright (c) 2026 Stephen Cobeldick
% Simulation Matrices Only:
% Copyright (c) 2009 Gustavo M. Machado, Manuel M. Oliveira, and Leandro A. F. Fernandes.
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%license