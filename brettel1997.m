function [cvd,raw,lms,sim,opts] = brettel1997(rgb,typ,opts,varargin)
% Simulate dichromat color vision deficiency (CVD) using Brettel et al. (1997).
%
% (c) 2026 Stephen Cobeldick
%
% Simulates the perceived colors seen by a dichromat observer using the
% colorimetric projection method of Brettel, Viénot & Mollon (1997). Each
% input color is converted to an LMS-like cone-response space, projected
% parallel to the missing-cone axis onto one of two reduced-stimulus
% half-planes, and converted back to sRGB.
%
%%% Syntax %%%
%
%   cvd = brettel1997(rgb,typ)
%   cvd = brettel1997(rgb,typ,opts)
%   cvd = brettel1997(rgb,typ,<name-value pairs>)
%   [cvd,raw] = brettel1997(...)
%   [cvd,raw,lms,sim,opts] = brettel1997(...)
%
%% Algorithm %%
%
% Reference:
%  Brettel H, Viénot F, Mollon J D: "Computerized simulation of color
%  appearance for dichromats", JOSA A 14(10):2647-2655, 1997.
%  <https://vision.psychol.cam.ac.uk/jdmollon/papers/Dichromatsimulation.pdf>
%
% The original Brettel et al. implementation used a calibrated CRT display:
% monitor RGB values were converted to LMS values using measured spectral
% power distributions of the CRT primaries, transformed in LMS space, and
% converted back to monitor RGB. This implementation accepts and returns
% sRGB by default. This is a practical deviation from the reference paper:
% sRGB values are linearized, converted to CIE XYZ using the IEC 61966-2-1
% matrix, converted to an LMS-like space using the default HPE matrix,
% then projected using Brettel et al.'s reduced-stimulus geometry.
%
% Method:
% 1. Convert the input sRGB colors to linear sRGB.
% 2. Convert linear sRGB to CIE XYZ.
% 3. Convert CIE XYZ to LMS-like cone responses.
% 4. Select one of two reduced-stimulus half-planes according to the
%    position of the stimulus relative to the neutral axis:
%      protan: if S/M < S_E/M_E, use the 575 nm anchor, otherwise 475 nm.
%      deutan: if S/L < S_E/L_E, use the 575 nm anchor, otherwise 475 nm.
%      tritan: if M/L < M_E/L_E, use the 660 nm anchor, otherwise 485 nm.
% 5. Project each stimulus parallel to the missing-cone axis onto the selected plane.
% 6. Convert the simulated LMS values back to CIE XYZ, linear RGB, and sRGB.
%
% Given neutral-axis stimulus E and a monochromatic anchor A, each reduced-
% stimulus half-plane contains the origin, E, and A. The plane equation is
%
%      a*L + b*M + c*S = 0
%
% where [a,b,c] = cross(E,A). The missing cone coordinate is replaced by
% the value which satisfies that plane equation, while the two remaining
% cone coordinates are preserved.
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
% Field   | Permitted | Description:
% Name:   | Values:   |
% ========|===========|====================================================
% gamma   | logical   | Apply inverse/forward sRGB gamma correction. (true**)
% --------|-----------|----------------------------------------------------
% wpt     | 1x3 float | XYZ reference white point, scaled Y==1.
%         |           | Default is D65: [0.95047,1,1.08883].
% --------|-----------|----------------------------------------------------
% RGB2XYZ | 3x3 float | Matrix converting linear RGB to XYZ.
%         |           | Default is IEC 61966-2-1 sRGB.
% --------|-----------|----------------------------------------------------
% XYZ2LMS | 3x3 float | Matrix converting XYZ to LMS-like responses.
%         |           | Default is the Hunt-Pointer-Estevez matrix used
%         |           | by CIECAM02.
% --------|-----------|----------------------------------------------------
% anchXYZ | 4x3 float | Monochromatic anchor XYZ values. Rows correspond
%         |           | to wavelengths [575;475;660;485] nm. These values
%         |           | are converted to LMS using <XYZ2LMS>. Only their
%         |           | chromatic directions are relevant to the planes.
% --------|-----------|----------------------------------------------------
% anchLMS | 4x3 float | Monochromatic anchor LMS values. Rows correspond
%         | []**      | to wavelengths [575;475;660;485] nm. If non-empty,
%         |           | these override <anchXYZ> and are used directly.
%
%% Examples %%
%
%%% Simulate a single RGB triple for protanopia %%%
%
%   >> brettel1997([1,0,0], 'protan')
%   ans = [0.5122   0.4356   0]
%
%%% View an image as seen by a deuteranope %%%
%
%   >> I = imread("peppers.png");
%   >> imshow(brettel1997(I, 'd'))
%
%%% View Parula colormap under all three CVD types %%%
%
%   >> I = permute(parula(17),[3,1,2]);
%   >> imshow([I;brettel1997(I,'p');brettel1997(I,'d');brettel1997(I,'t')])
%
%%% Use custom LMS anchors %%%
%
%   >> opt.anchLMS = [...
%   ..     0.9594,0.8884,0.0018;... 575nm
%   ..     0.0330,0.1006,1.0419;... 475nm
%   ..     0.1061,0.0344,0.0000;... 660nm
%   ..     0.0953,0.2130,0.6162]; % 485nm
%   >> imshow(brettel1997(I, 'deutan', opt))
%
%% Notes %%
%
% * This implementation models dichromacy only. It does not implement a
%   severity parameter for anomalous trichromacy.
% * Unlike the reference paper, this implementation defaults to sRGB input
%   and output, rather than a specifically calibrated CRT monitor. Users
%   requiring a different display model may supply the conversion matrices
%   <RGB2XYZ>, <XYZ2LMS>, the whitepoint <wpt>, and/or <anchLMS>.
% * The default <XYZ2LMS> matrix is the Hunt-Pointer-Estevez matrix defined
%   by CIECAM02. It is used here as a practical LMS-like cone-response
%   space, not as a full CIECAM02 appearance model.
% * The default monochromatic anchor XYZ values are tabulated CIE 1931
%   2-degree color-matching values at 575, 475, 660, and 485 nm. The scale
%   of each anchor is irrelevant because each anchor defines a plane
%   through the origin and the neutral axis.
% * Branch selection is implemented by cross-multiplication rather than
%   explicit division, e.g. S*ME < SE*M instead of S/M < SE/ME. This avoids
%   division-by-zero edge cases on artificial or out-of-gamut values.
% * The 2nd output <raw> is not clipped to the display gamut. The 1st
%   output is clipped to 0..1 for floating-point input, or cast back to
%   the input integer class for integer input.
%
%% Input Arguments %%
%
%   rgb = NumericArray of sRGB values to convert. Floating point values
%         must be 0<=rgb<=1, integer must be 0<=rgb<=intmax(class(rgb)).
%         Size Nx3 or RxCx3, the last dimension encodes the R,G,B values.
%   typ = CharRowVector or StringScalar, the type of dichromacy to simulate:
%        'p' / 'protan' / 'protanopia'     (L-cone absence).
%        'd' / 'deutan' / 'deuteranopia'   (M-cone absence).
%        't' / 'tritan' / 'tritanopia'     (S-cone absence).
%   opts = StructureScalar, optional parameter values as per 'Options' above.
%   <name-value pairs> = a comma-separated list of names and corresponding values.
%
%% Output Arguments %%
%
%   cvd  = NumericArray, the same size and class as <rgb>, containing the
%          simulated dichromat colors. Float values are clipped to 0..1.
%   raw  = FloatArray, the same size as <rgb>, containing the simulated
%          dichromat colors without clipping (i.e. values may be outside 0..1).
%   lms  = FloatArray, the same size as <rgb>, containing the input LMS
%          values used by the simulation.
%   sim  = FloatArray, the same size as <rgb>, containing the simulated LMS
%          values after Brettel projection.
%   opts = StructureScalar, the used parameter values as per 'Options' above.
%
%% Dependencies %%
%
% * MATLAB R2009b or later.
% * No toolboxes are required.
%
% See also CVDSIM DALTONIZER MACHADO2010 MILIC2015
% IMSHOW PARULA LINES COLORMAP COLORORDER BREWERMAP MAXDISTCOLOR

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
	error('SC:brettel1997:rgb:NotNumeric',...
		'1st input <rgb> must be a numeric array, not %s',class(rgb))
end
assert(isreal(rgb),...
	'SC:brettel1997:rgb:NotReal',...
	'1st input <rgb> must be a real array (not complex).')
assert(isz(end)==3 || isequal(isz,[3,1]),...
	'SC:brettel1997:rgb:InvalidSize',...
	'1st input <rgb> last dimension must have size 3 (e.g. Nx3 or RxCx3).')
assert(all(0<=rgb(:)&rgb(:)<=1),'SC:brettel1997:rgb:OutOfRange',...
	'1st input <rgb> values must be 0<=rgb<=%d',mxv)
rgb = reshape(rgb,[],3);
%
typ = mSS2C(typ);
assert(ischar(typ)&&ndims(typ)==2&&size(typ,1)==1,...
	'SC:brettel1997:typ:NotText',...
	'Second input <typ> must be a character vector or a string scalar.') %#ok<ISMAT>
%
switch lower(typ)
	case {'p','protan','protanopia'}
		dst = 'p';
	case {'d','deutan','deuteranopia'}
		dst = 'd';
	case {'t','tritan','tritanopia'}
		dst = 't';
	case {'protanomaly','deuteranomaly','tritanomaly'}
		error('SC:brettel1997:typ:AnomalousNotSupported',...
			'Second input <typ> "%s" is not supported: Brettel et al. (1997) defines dichromatic simulations only.',typ)
	otherwise
		error('SC:brettel1997:typ:NotSupported',...
			'Second input <typ> "%s" is not supported: use "protan"/"deutan"/"tritan" or their dichromatic full names or initials.',typ)
end
%
stpo = struct(... Default option values.
	'gamma',true,...
	'wpt',[0.95047,1,1.08883],...
	'RGB2XYZ',[... IEC 61966-2-1:1999, used for compatibility with other implementations.
	0.4124,0.3576,0.1805;...
	0.2126,0.7152,0.0722;...
	0.0193,0.1192,0.9505],...
	'XYZ2LMS',[... Hunt-Pointer-Estevez, as used by CIECAM02.
	+0.38971,+0.68898,-0.07868;...
	-0.22981,+1.18340,+0.04641;...
	+0      ,+0      ,+1      ],...
	'anchXYZ',[... CIE 1931 2-degree CMF values. Rows: 575, 475, 660, 485 nm.
	0.842500,0.915400,0.001800;...
	0.142100,0.112600,1.041900;...
	0.164900,0.061000,0.000000;...
	0.057950,0.169300,0.616200],...
	'anchLMS',[]);
%
% Check any supplied option field names and values:
switch nargin
	case 2 % no user-supplied options
		% Use defaults.
	case 3 % options in a struct
		assert(isstruct(opts)&&isscalar(opts),...
			'SC:brettel1997:options:NotScalarStruct',...
			'Third input <opts> must be a scalar structure, or options must be supplied as name-value pairs.')
		opts = structfun(@mSS2C,opts,'UniformOutput',false);
		stpo = mOptions(stpo,opts);
	otherwise % options as <name-value> pairs
		tmp = [{opts},varargin];
		assert(mod(numel(tmp),2)==0,...
			'SC:brettel1997:options:NameValuePairsNotPaired',...
			'Options supplied as name-value pairs must have one value for every name.')
		tmp = cellfun(@mSS2C,tmp,'UniformOutput',false);
		opts = cell2struct(tmp(2:2:end),tmp(1:2:end),2);
		stpo = mOptions(stpo,opts);
end
opts = stpo;
%
%% Prepare Simulation Parameters %%
%
E = stpo.wpt * stpo.XYZ2LMS.';
%
if isempty(stpo.anchLMS)
	anc = stpo.anchXYZ * stpo.XYZ2LMS.';
else
	anc = stpo.anchLMS;
end
%
% Plane normals [a,b,c] = cross(E,A) for anchors [575;475;660;485] nm.
plane = [...
	E(2).*anc(:,3) - E(3).*anc(:,2),...
	E(3).*anc(:,1) - E(1).*anc(:,3),...
	E(1).*anc(:,2) - E(2).*anc(:,1)];
%
% Check that the anchors and neutral axis define usable Brettel planes.
% Each anchor is used only to define one plane with the neutral axis, so
% the anchor matrix itself need not be full rank. What matters is that
% each row defines a non-degenerate plane, and that when replacing the
% missing cone response the denominator is not zero.
anm = hypot(hypot(anc(:,1),anc(:,2)),anc(:,3));
pnm = hypot(hypot(plane(:,1),plane(:,2)),plane(:,3));
enm = hypot(hypot(E(1),E(2)),E(3));
tol = 100 * eps(max([1;enm;anm;pnm]));
%
assert(all(anm>0),...
	'SC:brettel1997:options:anchors:ZeroAnchor',...
	'Each anchor row must contain at least one non-zero value.')
assert(all(pnm > tol*enm.*anm),...
	'SC:brettel1997:options:anchors:DegeneratePlane',...
	'Each anchor must define a non-degenerate plane with the white point.')
assert(all(abs(plane(1:2,1)) > tol*pnm(1:2)),...
	'SC:brettel1997:options:anchors:DegenerateProtanPlane',...
	'The 575 nm and 475 nm anchors must define protan planes with non-zero L coefficients.')
assert(all(abs(plane(1:2,2)) > tol*pnm(1:2)),...
	'SC:brettel1997:options:anchors:DegenerateDeutanPlane',...
	'The 575 nm and 475 nm anchors must define deutan planes with non-zero M coefficients.')
assert(all(abs(plane(3:4,3)) > tol*pnm(3:4)),...
	'SC:brettel1997:options:anchors:DegenerateTritanPlane',...
	'The 660 nm and 485 nm anchors must define tritan planes with non-zero S coefficients.')
%
%% Convert to LMS %%
%
if stpo.gamma
	lin = sGammaInv(rgb);
else
	lin = rgb;
end
%
XYZ = lin * stpo.RGB2XYZ.';
lms = XYZ * stpo.XYZ2LMS.';
sim = lms;
%
%% Apply Brettel et al. Projection %%
%
switch dst
	case 'p' % protanopia: replace L, preserve M and S.
		idx = lms(:,3).*E(2) < E(3).*lms(:,2); % S/M < S_E/M_E
		idy = 1 + ~idx; % true -> 575 nm, false -> 475 nm.
		a = plane(idy,1);
		b = plane(idy,2);
		c = plane(idy,3);
		sim(:,1) = -(b.*lms(:,2) + c.*lms(:,3)) ./ a;
	case 'd' % deuteranopia: replace M, preserve L and S.
		idx = lms(:,3).*E(1) < E(3).*lms(:,1); % S/L < S_E/L_E
		idy = 1 + ~idx; % true -> 575 nm, false -> 475 nm.
		a = plane(idy,1);
		b = plane(idy,2);
		c = plane(idy,3);
		sim(:,2) = -(a.*lms(:,1) + c.*lms(:,3)) ./ b;
	case 't' % tritanopia: replace S, preserve L and M.
		idx = lms(:,2).*E(1) < E(2).*lms(:,1); % M/L < M_E/L_E
		idy = 3 + ~idx; % true -> 660 nm, false -> 485 nm.
		a = plane(idy,1);
		b = plane(idy,2);
		c = plane(idy,3);
		sim(:,3) = -(a.*lms(:,1) + b.*lms(:,2)) ./ c;
	otherwise
		error('SC:brettel1997:typ:MissingCase','Please report this bug: missing case.')
end
%
%% Convert Back to sRGB %%
%
XYZ = sim / stpo.XYZ2LMS.';
lin = XYZ / stpo.RGB2XYZ.';
%
if stpo.gamma
	raw = reshape(sGammaCor(lin),isz);
else
	raw = reshape(lin,isz);
end
%
lms = reshape(lms,isz);
sim = reshape(sim,isz);
%
if mxv>1
	cvd = cast(mxv*raw,icl);
else
	cvd = min(1,max(0,raw));
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%brettel1997
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
		error('SC:brettel1997:options:UnknownOptionName',...
			'Unknown option: <%s>.\nOptions are:%s.',ofn,ont(2:end))
	elseif nnz(oix)>1
		dnt = sprintf(', <%s>',ofc{oix});
		error('SC:brettel1997:options:DuplicateOptionNames',...
			'Duplicate option names:%s.',dnt(2:end))
	end
	arg = opts.(ofn);
	dfn = dfc{dix};
	switch dfn
		case 'gamma'
			mLogical()
		case {'RGB2XYZ','XYZ2LMS'}
			mFullRank(false)
		case 'anchXYZ'
			mMatrix(false,[4,3])
		case 'anchLMS'
			mMatrix(true,[4,3])
		case 'wpt'
			mWhitePoint(false)
		otherwise
			error('SC:brettel1997:options:MissingCase','Please report this bug.')
	end
	stpo.(dfn) = arg;
end
%
%% Nested Functions %%
%
	function mLogical() % scalar logical.
		assert(isequal(arg,false)||isequal(arg,true),...
			sprintf('SC:brettel1997:%s:NotScalarLogical',dfn),...
			'The <%s> value must be true/1 or false/0.',dfn)
		arg = logical(arg);
	end
%
	function mWhitePoint(imt) % positive real 1x3 white point vector.
		mMatrix(imt,3)
		assert(all(arg(:)>0) && arg(2)==1,...
			'SC:brettel1997:wpt:InvalidValue',...
			'The <wpt> value must contain three positive numeric values and wpt(2)==1.')
		arg = reshape(arg,1,3);
	end
%
	function mFullRank(imt) % real finite full-rank 3x3 matrix.
		mMatrix(imt,[3,3])
		assert(rank(arg)==3,...
			sprintf('SC:brettel1997:%s:NotFullRank',dfn),...
			'The <%s> value must be a full-rank 3x3 matrix.',dfn)
	end
%
	function mMatrix(imt,siz) % real finite numeric array of specified size, optionally [].
		if imt && isnumeric(arg) && isempty(arg)
			arg = [];
			return
		end
		assert(isnumeric(arg),...
			sprintf('SC:brettel1997:%s:NotNumeric',dfn),...
			'The <%s> value must be numeric.',dfn)
		assert(isreal(arg),...
			sprintf('SC:brettel1997:%s:NotRealNumeric',dfn),...
			'The <%s> value cannot be complex.',dfn)
		if isscalar(siz)
			assert(numel(arg)==siz,...
				sprintf('SC:brettel1997:%s:InvalidSize',dfn),...
				'The <%s> value must contain %d numeric values.',dfn,siz)
		else
			assert(isequal(size(arg),siz),...
				sprintf('SC:brettel1997:%s:InvalidSize',dfn),...
				'The <%s> value must have size %s.',dfn,mat2str(siz))
		end
		assert(all(isfinite(arg(:))),...
			sprintf('SC:brettel1997:%s:NotFiniteNumeric',dfn),...
			'The <%s> value must contain only finite values.',dfn)
		arg = double(arg);
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
% Copyright (c) 1997 Hans Brettel, Francoise Vienot, and John D. Mollon.
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%license