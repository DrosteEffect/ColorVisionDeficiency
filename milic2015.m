function [rec,raw,ctr0,ctr1,idx,opts] = milic2015(rgb,typ,opts,varargin)
% Naturalness-preserving image rgb color enhancement for CVD/dichromats (Milić et al, 2015)
%
% (c) 2026 Stephen Cobeldick
%
% Recolors an RGB image to improve color distinction for dichromats
% (protanopes, deuteranopes, and tritanopes), using the content-dependent
% naturalness-preserving method of Milić, Hoffmann, Tómács, Novaković &
% Milosavljević (2015). The method segments the image by chromaticity in
% the CIE L*u'v' space, remaps representative segment centers by rotating
% them around the relevant dichromatic confusion point, and then applies
% the resulting center displacement to all pixels in each segment.
%
%%% Syntax %%%
%
% rec = milic2015(rgb,typ)
% rec = milic2015(rgb,typ,opts)
% rec = milic2015(rgb,typ,<name-value pairs>)
% [rec,raw] = milic2015(...)
% [rec,raw,ctr0,ctr1,idx,opts] = milic2015(...)
%
%% Algorithm %%
%
% Reference:
% Milić N, Hoffmann M, Tómács T, Novaković D, Milosavljević B:
% "A Content-Dependent Naturalness-Preserving Daltonization Method for
% Dichromatic and Anomalous Trichromatic Color Vision Deficiencies",
% Journal of Imaging Science and Technology 59(1):010504, 2015.
% <https://doi.org/10.2352/J.ImagingSci.Technol.2015.59.1.010504>
%
% Method:
% 1. Convert the input sRGB image colors to CIE L*u'v'.
% 2. Segment pixels by k-means clustering their u'v' chromaticities.
% 3. Compute one representative L*u'v' center for each non-empty segment.
% 4. Compute angular intervals around each center, measured from the
% selected dichromatic confusion point.
% 5. Redistribute center angles inside their admissible intervals using
% the stable-position equation described by the paper.
% 6. Rotate each center to its new angle, preserving its distance from the
% confusion point in the u'v' chromaticity diagram.
% 7. Recolor all pixels by preserving their original L*u'v' offset from
% their segment center.
% 8. Convert the recolored L*u'v' colors back to sRGB.
%
% The paper does not define automatic rules for choosing the number of
% segments, k-means settings, or the admissible angular half-width. These
% are therefore user options in this implementation. This implementation
% supports only the dichromatic confusion-point geometry explicitly given
% by the paper. Anomalous trichromacy is not supported.
%
%% Options %%
%
% The options may be supplied either
% 1) in a scalar structure, or
% 2) as a comma-separated list of name-value pairs.
%
% Field names are case-insensitive. The following field names are permitted
% as options (**=default value):
%
% Field | Permitted |
% Name: | Values:   | Description:
% ======|===========|======================================================
% nseg  | scalar    | Number of chromatic k-means segments.     (6**)
% ------|-----------|------------------------------------------------------
% ang   | scalar    | Admissible angular half-width in degrees. (5**)
% ------|-----------|------------------------------------------------------
% kmi   | scalar    | Maximum k-means iterations per replicate. (100**)
% ------|-----------|------------------------------------------------------
% kmr   | scalar    | Number of k-means replicates.             (3**)
% ------|-----------|------------------------------------------------------
% tol   | scalar    | Stable-angle solver tolerance.            (1e-12**)
% ------|-----------|------------------------------------------------------
% nit   | scalar    | Stable-angle solver maximum iterations.
%       | []**      | []** uses max(100,20*M), where M is the number
%       |           | of non-empty segment centers.
% ------|-----------|------------------------------------------------------
% wpt   | 1x3 float | XYZ reference white point, scaled Y==1.
%
%% Examples %%
%
%%% Recolor an image %%%
%
% >> I = imread("peppers.png");
% >> imshow(milic2015(I,'protan'))
% >> imshow(milic2015(I,'deutan', 'nseg',4, 'ang',4))
%
%%% Use a scalar structure of options %%%
%
% >> opt.nseg = 4;
% >> opt.ang = 4;
% >> imshow(milic2015(I,'p',opt))
%
%% Notes %%
%
% * This implementation is intentionally restricted to RGB images of size
%   RxCx3. Nx3 colormaps are not supported because the method is defined
%   using image segmentation.
% * This implementation performs chromatic remapping only. The optional
%   lightness remapping idea mentioned by the paper is not implemented
%   because the paper does not specify an algorithm for it.
% * Anomalous trichromacy is not supported. Although the paper states that
%   the concept can be adapted to anomalous trichromacy, it does not
%   define the required anomalous confusion-line geometry.
% * Empty k-means clusters are silently removed. This avoids arbitrary
%   reseeding rules and keeps the effective set of segment centers equal
%   to the non-empty image segments.
% * This implementation leaves the random number generator untouched:
%   users who require reproducible segmentation should set and restore
%   RNG state outside this function using whatever seed-control API is
%   appropriate for their MATLAB release.
%
%% Input Arguments %%
%
%   rgb = NumericArray of sRGB values to convert, size RxCx3. Floating
%         point values must be 0<=rgb<=1, integer values must be >=0.
%         Dimensions 1 and 2 are interpreted as rows and columns
%         respectively, dimension 3 encodes the R,G,B values.
%   typ = CharRowVector or StringScalar, the type of dichromacy to correct for:
%        'p' / 'protan' / 'protanopia'   (L-cone absence).
%        'd' / 'deutan' / 'deuteranopia' (M-cone absence).
%        't' / 'tritan' / 'tritanopia'   (S-cone absence).
%   opts = StructureScalar, optional parameter values as per 'Options' above.
%   <name-value pairs> = a comma-separated list of names and corresponding values.
%
%% Output Arguments %%
%
%   rec  = NumericArray, the same size and class as <rgb>, the recolored
%          image. Float values are clipped to 0..1.
%   raw  = FloatArray, the same size as <rgb>, the recolored image without
%          clipping (i.e. values may be outside 0..1, depending on LUV2RGB).
%   ctr0 = FloatArray, size Mx3, the original non-empty segment centers in
%          CIE L*u'v', where M is the effective number of non-empty segments.
%   ctr1 = FloatArray, size Mx3, the remapped segment centers in CIE L*u'v'.
%   idx  = NumericArray, size RxC, the effective segment index of each pixel.
%   opts = StructureScalar, the used parameter values as per 'Options' above.
%
%% Dependencies %%
%
% * MATLAB R2009b or later.
% * No toolboxes are required.
%
% See also CVDSIM DALTONIZER MACHADO2010
% COLORMAP COLORORDER BREWERMAP MAXDISTCOLOR
% SRGB_TO_CAM02UCS CAM02UCS_TO_SRGB SRGB_TO_CAM16UCS CAM16UCS_TO_SRGB

% Release | Feature
% --------|--------
% R2016b  | string class [only if supplied]
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
	error('SC:milic2015:rgb:NotNumeric',...
		'1st input <rgb> must be a numeric array, not %s',class(rgb))
end
assert(isreal(rgb),...
	'SC:milic2015:rgb:NotReal',...
	'1st input <rgb> must be a real array (not complex).')
assert(ndims(rgb)==3 && isz(3)==3,...
	'SC:milic2015:rgb:InvalidSize',...
	'1st input <rgb> must be an RxCx3 RGB image. Nx3 color lists and N-D arrays are not supported.') %#ok<ISMAT>
assert(all(0<=rgb(:)&rgb(:)<=1),'SC:milic2015:rgb:OutOfRange',...
	'1st input <rgb> values must be 0<=rgb<=%d',mxv)
%
typ = mSS2C(typ);
assert(ischar(typ)&&ndims(typ)==2&&size(typ,1)==1,...
	'SC:milic2015:typ:NotText',...
	'Second input <typ> must be a character vector or a string scalar.') %#ok<ISMAT>
%
switch lower(typ)
	case {'p','protan','protanopia'}
		cnf = [+0.68,+0.50];
	case {'d','deutan','deuteranopia'}
		cnf = [-1.22,+0.78];
	case {'t','tritan','tritanopia'}
		cnf = [+0.26,+0.00];
	case {'protanomaly','deuteranomaly','tritanomaly'}
		error('SC:milic2015:typ:AnomalousNotSupported',...
			'Second input <typ> "%s" is not supported: this implementation uses dichromatic confusion points only.',typ)
	otherwise
		error('SC:milic2015:typ:NotSupported',...
			'Second input <typ> "%s" is not supported: use "protan"/"deutan"/"tritan" or their dichromatic full names or initials.',typ)
end
%
stpo = struct(... Default option values.
	'nseg',6, 'ang',5, 'kmi',100, 'kmr',3,...
	'tol',1e-12, 'nit',[], 'wpt',[0.95047,1,1.08883]);
%
% Check any supplied option field names and values:
switch nargin
	case 2 % no user-supplied options
		% Use defaults.
	case 3 % options in a struct
		assert(isstruct(opts)&&isscalar(opts),...
			'SC:milic2015:options:NotScalarStruct',...
			'Third input <opts> must be a scalar structure, or options must be supplied as name-value pairs.')
		opts = structfun(@mSS2C,opts,'UniformOutput',false);
		stpo = mOptions(stpo,opts);
	otherwise % options as <name-value> pairs
		temp = cellfun(@mSS2C,[{opts},varargin],'UniformOutput',false);
		opts = cell2struct(temp(2:2:end),temp(1:2:end),2);
		stpo = mOptions(stpo,opts);
end
opts = stpo;
%
assert(stpo.nseg<=isz(1)*isz(2),...
	'SC:milic2015:options:nseg:TooLarge',...
	'The <nseg> value must not exceed the number of image pixels.')
assert(stpo.ang<=90,...
	'SC:milic2015:options:ang:TooLarge',...
	'The angular half-width <ang> must not exceed 90 degrees.')
%
%% Convert to CIE L*u''v'' and Segment by Chromaticity %%
%
luv = sRGB2Lupvp(reshape(rgb,[],3),stpo.wpt);
uv0 = luv(:,2:3);
%
% The paper specifies k-means segmentation by chromaticity in the u'v'
% plane, excluding lightness. It does not specify an initialization method,
% empty-cluster handling, or replicate heuristic. This implementation uses
% a small local k-means routine, keeps the RNG under user control, and
% removes empty clusters after clustering.
%
idx = mKMeans(uv0,stpo.nseg,stpo.kmi,stpo.kmr);
%
ctr0 = zeros(stpo.nseg,3);
use = false(stpo.nseg,1);
for k = 1:stpo.nseg
	idk = idx==k;
	if any(idk)
		ctr0(k,:) = mean(luv(idk,:),1);
		use(k) = true;
	end
end
%
map = zeros(stpo.nseg,1);
map(use) = 1:nnz(use);
idx = map(idx);
ctr0 = ctr0(use,:);
nctr = size(ctr0,1);
%
%% Remap Segment Centers %%
%
ctr1 = ctr0;
if nctr>1 && stpo.ang>0
	% The paper defines angular remapping intervals but does not define a
	% baseline-selection algorithm. We therefore choose the largest angular
	% gap as the branch cut, giving a compact unwrapped angular sequence.
	% For real sRGB colors and sane angular half-widths this avoids interval
	% wraparound; the option check above rejects extreme artificial widths.
	uvc = ctr0(:,2:3);
	th0 = atan2(uvc(:,2)-cnf(2),uvc(:,1)-cnf(1));
	rad = hypot(uvc(:,1)-cnf(1),uvc(:,2)-cnf(2));
	[ths,ord0] = sort(th0);
	gap = diff([ths;ths(1)+2*pi]);
	[~,idg] = max(gap);
	ord = ord0([idg+1:nctr,1:idg]);
	p0 = unwrap(th0(ord));
	%
	% Milić et al. demonstrate fixed admissible angles (e.g. 4 and 5
	% degrees) but do not define an automatic rule for choosing them.
	% The user option <ang> is therefore used directly as the symmetric
	% angular half-width around every original center.
	%
	% The paper describes an elliptical admissible region related to MacAdam
	% ellipses, but does not provide the ellipse axes or an automatic rule for
	% deriving them. This implementation therefore uses a constant-radius
	% angular approximation: each center is rotated around the confusion point
	% while preserving its original distance from that point.
	dth = pi*stpo.ang/180;
	aa = p0 - dth;
	bb = p0 + dth;
	%
	p1 = mStableAngles(p0,aa,bb,stpo.tol,stpo.nit);
	%
	% Apply the angular displacement on circular arcs centered at the
	% confusion point. This is the explicit geometric simplification relative
	% to the paper's under-specified elliptical admissible regions.
	ctr1(ord,2) = cnf(1) + rad(ord).*cos(p1);
	ctr1(ord,3) = cnf(2) + rad(ord).*sin(p1);
end
%
%% Recolor Pixels and Convert Back to sRGB %%
%
luv1 = luv;
dlt = ctr1 - ctr0;
for k = 1:nctr
	idk = idx==k;
	% Preserve each pixel's original offset from its segment center, as
	% specified by the paper. No clipping is applied in L*u'v' space.
	luv1(idk,:) = bsxfun(@plus,luv(idk,:),dlt(k,:));
end
%
raw = reshape(Lupvp2sRGB(luv1,stpo.wpt),isz);
idx = reshape(idx,isz(1),isz(2));
%
if mxv>1
	rec = cast(mxv*raw,icl);
else
	rec = min(1,max(0,raw));
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%milic2015
function p = mStableAngles(p0,aa,bb,tol,nit)
% Solve the stable angular redistribution from Milić et al. Eq. (5).
% The paper defines the fixed-point equation but not the numerical solver.
% We use bounded Gauss-Seidel iteration because it is simple, robust,
% directly applies the equation, and cannot hang due to the iteration
% limit. Because this updates p(k) in place, it is not a literal
% simultaneous Jacobi update; in pathological overlapping-interval cases
% this could affect the exact numerical equilibrium, although the solution
% remains bounded by the admissible intervals.
%
% Monotonicity is structurally preserved: p0 is sorted before this solver
% is called, and the admissible intervals are formed with non-negative
% half-widths, so adjacent intervals can overlap but cannot invert.
n = numel(p0);
p = min(bb,max(aa,p0));
p(1) = aa(1);
p(n) = bb(n);
if isempty(nit)
	nit = max(100,20*n);
end
for it = 1:nit
	old = p;
	for k = 2:n-1
		p(k) = min(bb(k),max(aa(k),0.5*(p(k-1)+p(k+1))));
	end
	if max(abs(p-old))<=tol
		break
	end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mStableAngles
function idx = mKMeans(dat,k,kmi,kmr)
% Small toolbox-free k-means implementation for Nx2 chromaticity data.
% Empty clusters keep their previous center during iteration; they are
% removed by the caller after clustering. This avoids arbitrary reseeding.
n = size(dat,1);
bst = Inf;
idx = ones(n,1);
for r = 1:kmr
	prm = randperm(n);
	ctr = dat(prm(1:k),:);
	tmp = zeros(n,1);
	for it = 1:kmi
		old = tmp;
		dst = zeros(n,k);
		for j = 1:k
			dxy = bsxfun(@minus,dat,ctr(j,:));
			dst(:,j) = dxy(:,1).^2 + dxy(:,2).^2;
		end
		[~,tmp] = min(dst,[],2);
		oldctr = ctr;
		for j = 1:k
			idj = tmp==j;
			if any(idj)
				ctr(j,:) = mean(dat(idj,:),1);
			end
		end
		if isequal(tmp,old) || max(abs(ctr(:)-oldctr(:)))<=eps
			break
		end
	end
	dst = zeros(n,k);
	for j = 1:k
		dxy = bsxfun(@minus,dat,ctr(j,:));
		dst(:,j) = dxy(:,1).^2 + dxy(:,2).^2;
	end
	[val,tmp] = min(dst,[],2);
	sse = sum(val);
	if sse<bst
		bst = sse;
		idx = tmp;
	end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mKMeans
function luv = sRGB2Lupvp(rgb,wpt)
%% RGB2XYZ
M = [... IEC 61966-2-1:1999, used for compatibility with other implementations.
	0.4124,0.3576,0.1805;...
	0.2126,0.7152,0.0722;...
	0.0193,0.1192,0.9505];
XYZ = sGammaInv(rgb) * M.';
%% XYZ2LUPVP
epsilon = 216/24389; % (6/29)^3
kappa = 24389/27; % (29/3)^3
yr = XYZ(:,2) ./ wpt(2);
idx = yr>epsilon;
L = kappa * yr;
L(idx) = 116*nthroot(yr(idx),3) - 16;
%
dnw = wpt(1) + 15*wpt(2) + 3*wpt(3);
un = 4*wpt(1) / dnw;
vn = 9*wpt(2) / dnw;
%
den = XYZ(:,1) + 15*XYZ(:,2) + 3*XYZ(:,3);
idu = abs(den)>eps;
up = un + zeros(size(den));
vp = vn + zeros(size(den));
up(idu) = 4*XYZ(idu,1)./den(idu);
vp(idu) = 9*XYZ(idu,2)./den(idu);
%
luv = [L,up,vp];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sRGB2Lupvp
function rgb = Lupvp2sRGB(luv,wpt)
%% LUPVP2XYZ
epsilon = 216/24389; % (6/29)^3
kappa = 24389/27; % (29/3)^3
L = luv(:,1);
up = luv(:,2);
vp = luv(:,3);
%
fy = (L+16)/116;
Y = wpt(2) * fy.^3;
idl = L<=(kappa*epsilon);
Y(idl) = wpt(2) * L(idl) / kappa;
%
% Inverting u'=4X/(X+15Y+3Z), v'=9Y/(X+15Y+3Z) gives:
% X = Y*9*u'/(4*v'), Z = Y*(12-3*u'-20*v')/(4*v').
% The denominator guard is purely defensive: valid sRGB-derived colors have
% positive v', but remapping and raw intermediate values are intentionally
% allowed to leave the display gamut before final clipping.
den = 4*vp;
den(abs(den)<eps) = eps;
X = Y .* 9 .* up ./ den;
Z = Y .* (12 - 3*up - 20*vp) ./ den;
%
idz = L<=0;
X(idz) = 0;
Y(idz) = 0;
Z(idz) = 0;
XYZ = [X,Y,Z];
%% XYZ2RGB
M = [... IEC 61966-2-1:1999 (for compatibility)
	0.4124,0.3576,0.1805;...
	0.2126,0.7152,0.0722;...
	0.0193,0.1192,0.9505];
rgb = sGammaCor(XYZ / M.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Lupvp2sRGB
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
function stpo = mOptions(stpo,opts)
% Options check: only supported fieldnames with suitable option values.
dfc = fieldnames(stpo);
ofc = fieldnames(opts);
%
for k = 1:numel(ofc)
	ofn = ofc{k};
	dix = strcmpi(ofn,dfc);
	oix = strcmpi(ofn,ofc);
	if ~any(dix)
		dfs = sort(dfc);
		ont = sprintf(', <%s>',dfs{:});
		error('SC:milic2015:options:UnknownOptionName',...
			'Unknown option: <%s>.\nOptions are:%s.',ofn,ont(2:end))
	elseif nnz(oix)>1
		dnt = sprintf(', <%s>',ofc{oix});
		error('SC:milic2015:options:DuplicateOptionNames',...
			'Duplicate option names:%s.',dnt(2:end))
	end
	arg = opts.(ofn);
	dfn = dfc{dix};
	switch dfn
		case 'nseg'
			mInteger(0)
		case {'kmi','kmr'}
			mInteger(0)
		case 'nit'
			mInteger(1)
		case 'ang'
			mScalar()
		case 'tol'
			mScalar()
		case 'wpt'
			mWhitePoint()
		otherwise
			error('SC:milic2015:options:MissingCase','Please report this bug.')
	end
	stpo.(dfn) = arg;
end
%
%% Nested Functions %%
%
	function mInteger(ise) % positive integer scalar.
		if ise && isnumeric(arg) && isempty(arg)
			arg = [];
			return
		end
		assert(isnumeric(arg)&&isscalar(arg),...
			sprintf('SC:milic2015:%s:NotScalarNumeric',dfn),...
			'The <%s> value must be a scalar numeric.',dfn)
		assert(isreal(arg),...
			sprintf('SC:milic2015:%s:NotRealNumeric',dfn),...
			'The <%s> value cannot be complex. Input: %g%+gi',dfn,real(arg),imag(arg))
		assert(isfinite(arg),...
			sprintf('SC:milic2015:%s:NotFiniteNumeric',dfn),...
			'The <%s> value must be finite. Input: %g',dfn,arg)
		assert(fix(arg)==arg && arg>=1,...
			sprintf('SC:milic2015:%s:NotPositiveInteger',dfn),...
			'The <%s> value must be a positive integer. Input: %g',dfn,arg)
		arg = double(arg);
	end
	function mScalar() % finite non-negative scalar.
		assert(isnumeric(arg)&&isscalar(arg),...
			sprintf('SC:milic2015:%s:NotScalarNumeric',dfn),...
			'The <%s> value must be a scalar numeric.',dfn)
		assert(isreal(arg),...
			sprintf('SC:milic2015:%s:NotRealNumeric',dfn),...
			'The <%s> value cannot be complex. Input: %g%+gi',dfn,real(arg),imag(arg))
		assert(isfinite(arg),...
			sprintf('SC:milic2015:%s:NotFiniteNumeric',dfn),...
			'The <%s> value must be finite. Input: %g',dfn,arg)
		assert(arg>=0,...
			sprintf('SC:milic2015:%s:NegativeNumeric',dfn),...
			'The <%s> value must be non-negative. Input: %g',dfn,arg)
		arg = double(arg);
	end
	function mWhitePoint() % positive real 1x3 vector.
		assert(isnumeric(arg),...
			'SC:milic2015:wpt:NotNumeric',...
			'The <wpt> value must be numeric.')
		assert(isreal(arg),...
			'SC:milic2015:wpt:NotRealNumeric',...
			'The <wpt> value cannot be complex.')
		assert(numel(arg)==3 && all(arg(:)>0) && arg(2)==1,...
			'SC:milic2015:wpt:InvalidValue',...
			'The <wpt> value must be a 1x3 numeric vector where 0<=wpt and wpt(2)==1.')
		arg = reshape(double(arg),1,3);
	end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mOptions
function arr = mSS2C(arr)
% If scalar string then extract the character vector, otherwise data is unchanged.
if isa(arr,'string') && isscalar(arr)
	arr = arr{1};
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mSS2C
%
% Code and Implementation:
% Copyright (c) 2026 Stephen Cobeldick
% Algorithm:
% Copyright (c) 2015 Neda Milić, Miklós Hoffmann, Tibor Tómács,
% Dragoljub Novaković, and Branko Milosavljević.
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%license