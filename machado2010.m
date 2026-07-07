function [rec,raw,cd1,pm1] = machado2010(rgb,typ,exg,cd0,pm0)
% Image color contrast enhancement for CVD/dichromats (Machado et al, 2010)
%
% (c) 2026 Stephen Cobeldick
%
% Recolors an RGB image to improve local color contrast for dichromats
% (protanopes, deuteranopes, and tritanopes), using the projection-based
% method of Machado & Oliveira (2010). The method estimates the direction
% of maximum local color-contrast loss in the CIELab L*a*b* chromaticity
% plane, projects the image colors onto that direction, and rotates the
% projected colors onto an approximate dichromat gamut plane.
%
%%% Syntax %%%
%
%   rec = machado2010(rgb,typ)
%   rec = machado2010(rgb,typ,exg)
%   rec = machado2010(rgb,typ,exg,cd0)
%   rec = machado2010(rgb,typ,exg,cd0,pm0)
%   [rec,raw,cd1,pm1] = machado2010(...)
%
%% Algorithm %%
%
% Method (Machado & Oliveira, 2010):
%  1. Convert the input sRGB image colors to CIE L*a*b*.
%  2. Estimate local contrast loss using one randomly chosen nearby pixel
%     for each pixel (Gaussian pairing in the image plane).
%  3. Compute the chromatic direction <cd1> maximizing the weighted loss of
%     local contrast, as the main eigenvector of a 2x2 scatter matrix.
%  4. Project the original colors onto the plane defined by the L* axis and <cd1>.
%  5. Rotate those projected chromatic coordinates onto the approximate
%     dichromat gamut plane.
%  6. Convert the recolored L*a*b* colors back to sRGB.
%
% This implementation is intentionally restricted to RGB images of size
% RxCx3. Unlike CVDSIM and many color-space conversions, this algorithm is
% not purely pointwise: the recoloring direction is estimated from local
% spatial color-contrast losses. The paper defines the sampling step using
% pixel neighborhoods, image width/height, and Gaussian-distributed row and
% column offsets. Inputs such as Nx3 colormaps therefore have no defined
% image-space neighborhood in this implementation.
%
% This implementation leaves the random number generator untouched: users
% who require reproducible recoloring, e.g. for tests or video processing,
% should set and restore the RNG state outside this function using whatever
% seed-control API is appropriate for their MATLAB release.
%
% References:
%  Machado G M, Oliveira M M: "Real-Time Temporal-Coherent Color Contrast
%  Enhancement for Dichromats", Computer Graphics Forum 29(3):933-942, 2010.
%  <https://doi.org/10.1111/j.1467-8659.2009.01586.x>
%  Gamut-plane angles from:
%  Kuhn G R, Oliveira M M, Fernandes L A F: "An efficient naturalness-
%  preserving image-recoloring method for dichromats", IEEE TVCG 14(6):1747-1754, 2008.
%
%% Examples %%
%
%%% Recolor an image %%%
%
%   >> I = imread("peppers.png");
%   >> imshow(machado2010(I,'deutan'))
%   >> imshow(machado2010(I,'protan',true)) % exagerated
%
%%% Preserve temporal coherence across frames %%%
%
%   >> cdt = [];
%   >> pmt = [];
%   >> for k = 1:numel(frm)
%   ..     [out,~,cdt,pmt] = machado2010(frm{k},'d',false,cdt,pmt);
%   ..     imshow(out)
%   ..     drawnow
%   .. end
%
%% Notes %%
%
% * This algorithm is described for dichromats. It does not implement a
%   severity parameter for anomalous trichromacy.
% * The input must be RxCx3. Nx3 colormaps/color lists are intentionally
%   not supported because the paper-defined Gaussian pairing step requires
%   a two-dimensional image-space neighborhood.
% * Degenerate images with R==1 or C==1 are accepted but produce a warning:
%   the Gaussian pairing then samples along a one-pixel-wide image domain,
%   which is a weak approximation of the intended 2D local-neighborhood
%   estimator.
% * <cd0>/<cd1> control only the sign of the chromatic direction between
%   frames. <pm0>/<pm1> additionally allow reuse of the same pixel-pair
%   map across frames, matching the paper's Section 3.1 statement that
%   pairs are pre-computed once and "used during the entire sequence". If
%   <pm0> does not match the current image size, a new pair map is silently
%   (re)computed and returned as <pm1>: this is expected to be rare, and a
%   caller intentionally varying frame size mid-sequence is assumed to be
%   doing something unusual enough to not need a warning here.
% * If every sampled pixel-pair has zero (or non-finite) measured contrast
%   loss, the image is left unrecolored for that call (see the code comment
%   above the relevant branch for the paper's justification for this choice).
% * For exaggerated contrast, chromaticity coordinates are scaled so that
%   the maximum output chroma is 148, following the value stated by the
%   paper. This is not the preferred/default use of the method.
%
%% Input Arguments (**=default) %%
%
%   rgb = NumericArray of sRGB values to convert, size RxCx3. Floating
%         point values must be 0<=rgb<=1, integer values must be >=0.
%         Dimensions 1 and 2 are interpreted as rows and columns
%         respectively, dimension 3 encodes the R,G,B values.
%   typ = CharRowVector or StringScalar, the type of dichromacy to correct for:
%         'p' / 'protan' / 'protanopia'   / 'protanomaly'   (L-cone deficiency).
%         'd' / 'deutan' / 'deuteranopia' / 'deuteranomaly' (M-cone deficiency).
%         't' / 'tritan' / 'tritanopia'   / 'tritanomaly'   (S-cone deficiency).
%   exg = LogicalScalar for selecting exaggerated contrast, where:
%         true    => exaggerated contrast.
%         false** => regular recoloring.
%   cd0 = NumericVector, size 1x2, the previous chromatic direction <cd1>.
%         If supplied, the current chromatic direction may be flipped to
%         avoid abrupt sign changes between consecutive frames.
%   pm0 = NumericVector, optional previous pixel-pair map, size R*Cx1, as
%         returned by a prior call as <pm1>. If supplied and its number of
%         elements matches R*C for the current <rgb>, the same pixel pairing
%         is reused (as required for exact video/sequence temporal coherence
%         per the paper); otherwise a new pair map is computed and returned.
%
%% Output Arguments %%
%
%   rec = NumericArray, the same size and class as <rgb>, the recolored
%         image. Float values are clipped to 0..1.
%   raw = FloatArray, the same size as <rgb>, the recolored image without
%         clipping (i.e. values may be outside 0..1, depending on LAB2RGB).
%   cd1 = FloatVector, size 1x2, the chromatic direction used for recoloring.
%   pm1 = NumericVector, size R*Cx1, the pixel-pair map used for this call.
%
%% Dependencies %%
%
% * MATLAB R2009b or later.
%
% See also CVDSIM DALTONIZER MILIC2015
% COLORMAP COLORORDER MAXDISTCOLOR BREWERMAP
% SRGB_TO_CAM02UCS CAM02UCS_TO_SRGB SRGB_TO_CAM16UCS CAM16UCS_TO_SRGB

% Release | Feature
% --------|--------
% R2016b  | string class                                 [only if supplied]
% R2009b  | tilde argument placeholder
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
	error('SC:machado2010:rgb:NotNumeric',...
		'1st input <rgb> must be a numeric array, not %s',class(rgb))
end
assert(isreal(rgb),...
	'SC:machado2010:rgb:NotReal',...
	'1st input <rgb> must be a real array (not complex).')
assert(ndims(rgb)==3 && isz(3)==3,...
	'SC:machado2010:rgb:InvalidSize',...
	'1st input <rgb> must be an RxCx3 RGB image. Nx3 color lists and N-D arrays are not supported.') %#ok<ISMAT>
assert(all(0<=rgb(:)&rgb(:)<=1),'SC:machado2010:rgb:OutOfRange',...
	'1st input <rgb> values must be 0<=rgb<=%d',mxv)
if isz(1)==1 || isz(2)==1
	warning('SC:machado2010:rgb:DegenerateImage',...
		'1st input <rgb> has R==1 or C==1: Gaussian pairing is intended for 2D images, so this is a degenerate image case.')
end
%
typ = svcSS2C(typ);
assert(ischar(typ)&&ndims(typ)==2&&size(typ,1)==1,...
	'SC:machado2010:typ:NotText',...
	'Second input <typ> must be a character vector or a string scalar.') %#ok<ISMAT>
%
switch lower(typ)
	case {'p','protan','protanopia','protanomaly'}
		ang = -11.48;
	case {'d','deutan','deuteranopia','deuteranomaly'}
		ang = -8.11;
	case {'t','tritan','tritanopia','tritanomaly'}
		ang = 46.37;
	otherwise
		error('SC:machado2010:typ:NotSupported',...
			'Second input <typ> "%s" is not supported: use "protan"/"deutan"/"tritan" or their full names or their initials.',typ)
end
%
if nargin<3 || isempty(exg)
	exg = false;
else
	assert(isequal(exg,0)||isequal(exg,1),...
		'SC:machado2010:exg:NotScalar',...
		'Third input <exg> must be a true/1 or false/0.')
	exg = logical(exg);
end
%
if nargin<4 || isempty(cd0)
	tcd = [];
else
	assert(isnumeric(cd0)&&isreal(cd0)&&numel(cd0)==2,...
		'SC:machado2010:cd0:InvalidSize',...
		'Fourth input <cd0> must be empty or a real numeric vector with two elements.')
	tcd = reshape(double(cd0),1,2);
end
%
if nargin<5 || isempty(pm0)
	pm0 = [];
else
	assert(isnumeric(pm0)&&isreal(pm0)&&isvector(pm0),...
		'SC:machado2010:pm0:InvalidSize',...
		'Fifth input <pm0> must be empty or a real numeric vector.')
	pm0 = double(pm0(:));
end
%
%% Convert to CIE L*a*b* and Estimate Lost Contrast %%
%
lab = sRGB2Lab(reshape(rgb,[],3));
ab0 = lab(:,2:3);
%
idx = localPairs(isz(1),isz(2),pm0);
pm1 = idx;
lab1 = lab;
lab2 = lab(idx,:);
%
the = pi*ang/180;
dch = [sin(the),cos(the)]; % chromatic axis of the approximate dichromat plane.
pr1 = lab1(:,2:3)*dch.';
pr2 = lab2(:,2:3)*dch.';
lab1(:,2:3) = pr1*dch;
lab2(:,2:3) = pr2*dch;
%
dlt = lab - lab(idx,:);
dsm = lab1 - lab2;
% dlt/dsm are 3-column (L*,a*,b*) differences: Eq. (2)'s norms are full Lab
% distances, so all three components must be included (nested HYPOT avoids
% intermediate overflow/underflow better than SQRT(SUM(X.^2,2))).
len = hypot(hypot(dlt(:,1),dlt(:,2)),dlt(:,3));
los = zeros(size(len));
idn = len>0;
% Eq. (2): relative loss of contrast. Pairs whose dichromat-perceived
% distance exceeds the original distance may contribute a negative
% (contrast-gain) term.
los(idn) = (len(idn)-hypot(hypot(dsm(idn,1),dsm(idn,2)),dsm(idn,3)))./len(idn);
vec = dlt(:,2:3) .* los(:,[1,1]);
%
%% Compute the Dominant Chromatic Loss Direction %%
%
A = vec.'*vec;
isdeg = ~(any(isfinite(A(:))) && any(A(:)));
if isdeg
	% Degenerate/no-contrast input: every sampled pixel pair had zero (or
	% non-finite) measured contrast loss, so there is no eigenvector to
	% extract. This case is not covered by an equation in the paper, but
	% Section 4.1 discusses the same underlying situation (Figure 7,
	% "Pink Head"): "Note that deuteranopes (and protanopes) already
	% perceive the reference image as having sufficient contrast, and no
	% recoloring is necessary." We extrapolate that guidance here: leave
	% the image unrecolored (cd1 has no meaningful direction) rather than
	% projecting onto an arbitrary fallback axis.
	cd1 = [0,0];
else
	[V,D] = eig(A);
	[~,idv] = max(abs(diag(D)));
	cd1 = V(:,idv).';
	cd1 = cd1 ./ max(eps,norm(cd1));
	%
	if ~isempty(tcd) && all(isfinite(tcd)) && any(tcd)
		tcd = tcd ./ max(eps,norm(tcd));
		if cd1*tcd.' < 0
			cd1 = -cd1;
		end
	end
end
%
%% Project, Rotate, and Convert Back to sRGB %%
%
if ~isdeg
	lab(:,2:3) = (ab0*cd1.') * dch;
	if exg
		chr = hypot(lab(:,2),lab(:,3));
		mxc = max(chr);
		if mxc>0
			lab(:,2:3) = lab(:,2:3) * (148/mxc);
		end
	end
end
%
raw = reshape(Lab2sRGB(lab),isz);
%
if mxv>1
	rec = cast(mxv*raw,icl);
else
	rec = min(1,max(0,raw));
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%machado2010
function idx = localPairs(hgt,wid,pm0)
% Random Gaussian pairing: choose one spatial neighbor per pixel, unless a
% valid precomputed pair map <pm0> from a previous frame is supplied, in
% which case it is reused as-is (paper Section 3.1: pairs are pre-computed
% once and "used during the entire sequence"). A mismatched <pm0> (e.g. the
% image size changed between calls) is assumed to be rare/exceptional, so a
% new pair map is silently (re)computed rather than raising a warning.
n = hgt*wid;
if numel(pm0)==n && all(pm0>=1 & pm0<=n & pm0==round(pm0))
	idx = pm0;
	return
end
[y,x] = ndgrid(1:hgt,1:wid);
% Per the paper (Section 3.1): horizontal/vertical offsets are drawn from
% a zero-mean univariate Gaussian with variance (2/pi)*sig2, where
% sig2 = sqrt(2*min(width,height)). randn() draws unit-variance samples,
% so the offsets must be scaled by the standard deviation, i.e. sqrt of
% that variance expression.
sig2 = sqrt(2*min(hgt,wid));
sig = max(1,sqrt((2/pi)*sig2));
x = round(x + sig*randn(hgt,wid));
y = round(y + sig*randn(hgt,wid));
x = max(1,min(wid,x));
y = max(1,min(hgt,y));
idx = sub2ind([hgt,wid],y(:),x(:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%localPairs
function Lab = sRGB2Lab(rgb)
%% RGB2XYZ
M = [... IEC 61966-2-1:1999, used for compatibility with other implementations.
	0.4124,0.3576,0.1805;...
	0.2126,0.7152,0.0722;...
	0.0193,0.1192,0.9505];
XYZ = sGammaInv(rgb) * M.';
%% XYZ2LAB
wpt = [0.95047,1,1.08883]; % (D65)
epsilon = 216/24389; % (6/29)^3
kappa   = 24389/27;  % (29/3)^3
% source: <http://www.brucelindbloom.com/index.html?LContinuity.html>
xyzr = bsxfun(@rdivide,XYZ,wpt);
idx  = xyzr>epsilon;
fxyz = (kappa*xyzr+16)/116;
fxyz(idx) = nthroot(xyzr(idx),3);
Lab = [116*fxyz(:,2)-16,...
	500*(fxyz(:,1)-fxyz(:,2)),...
	200*(fxyz(:,2)-fxyz(:,3))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sRGB2Lab
function rgb = Lab2sRGB(Lab)
%% LAB2XYZ
wpt = [0.95047,1,1.08883]; % (D65)
epsilon = 216/24389; % (6/29)^3
kappa   = 24389/27;  % (29/3)^3
% source: <http://www.brucelindbloom.com/index.html?Eqn_Lab_to_XYZ.html>
fxyz = bsxfun(@rdivide,Lab(:,[2,1,3]),[500,Inf,-200]);
fxyz = bsxfun(@plus,fxyz,(Lab(:,1)+16)/116);
tmp  = fxyz.^3;
idx  = tmp>epsilon;
idx(:,2) = true;
xyzr = idx.*tmp + ~idx.*(116*fxyz-16)/kappa;
idl  = Lab(:,1)>(kappa*epsilon);
xyzr(~idl,2) = Lab(~idl,1)/kappa;
XYZ = bsxfun(@times,xyzr,wpt);
%% YXZ2RGB
M = [... IEC 61966-2-1:1999 (for compatibility)
	0.4124,0.3576,0.1805;...
	0.2126,0.7152,0.0722;...
	0.0193,0.1192,0.9505];
rgb = sGammaCor(XYZ / M.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Lab2sRGB
function out = sGammaCor(inp)
% Forward Gamma correction: Nx3 linear RGB -> Nx3 sRGB.
idx = inp > 0.0031308;
out = 12.92 * inp;
out(idx) = real(1.055 * inp(idx).^(1./2.4) - 0.055);
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
% Algorithm:
% Copyright (c) 2010 Gustavo M. Machado and Manuel M. Oliveira.
% Gamut-plane Angle Values Only:
% Source: Kuhn, Oliveira, and Fernandes (2008), as cited by Machado & Oliveira (2010).
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%license