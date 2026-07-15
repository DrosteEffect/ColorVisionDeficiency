function [cvd,raw,lms,sim,opts] = vienot1999(rgb,typ,opts,varargin)
% Simulate protanopic/deuteranopic CVD using a modernized Vienot et al. (1999) method.
%
% (c) 2026 Stephen Cobeldick
%
% Simulates the perceived colors seen by protanopic or deuteranopic
% observers using the single-plane reduction idea of Vienot, Brettel and
% Mollon (1999), modified to use a modern sRGB/XYZ/LMS workflow (just like
% BRETTEL1997). Each input color is converted to an LMS-like cone-response
% space, projected parallel to the missing-cone axis onto the reduced plane
% containing [black, the reference white, the display blue primary] and
% is then converted back to sRGB.
%
%%% Syntax %%%
%
%   cvd = vienot1999(rgb,typ)
%   cvd = vienot1999(rgb,typ,opts)
%   cvd = vienot1999(rgb,typ,<name-value pairs>)
%   [cvd,raw,lms,sim,opts] = vienot1999(...)
%
%% Algorithm %%
%
% Reference:
%  Vienot F, Brettel H, Mollon J D: "Digital Video Colourmaps for Checking
%  the Legibility of Displays by Dichromats", Color Research and
%  Application, 24(4):243-252, 1999.
%
% The original Vienot et al. implementation was a practical digital-video
% colourmap method for protanopes and deuteranopes. It used ITU-R BT.709
% CRT primaries, D65 white, a CRT-style transfer function, Judd-Vos
% modified colorimetry, Smith-Pokorny fundamentals, and a fixed RGB-to-LMS
% matrix derived from those assumptions. This implementation instead uses
% an explicit modern color pipeline. By default, sRGB values are linearized,
% converted to CIE XYZ using the IEC 61966-2-1 matrix, converted to an
% LMS-like space using the default HPE matrix, projected using Vienot
% et al.'s single reduced-plane geometry, and converted back to sRGB.
%
% Method:
% 1. Convert the input sRGB colors to linear sRGB.
% 2. Convert linear sRGB to CIE XYZ.
% 3. Convert CIE XYZ to LMS-like cone responses.
% 4. Define the reduced plane as the plane through the origin, the reference
%    white, and the display blue primary. With neutral-axis stimulus E and
%    blue-primary stimulus B, the plane equation is
%
%       a*L + b*M + c*S = 0
%
%    where [a,b,c] = cross(E,B).
% 5. Project each stimulus parallel to the missing-cone axis onto this one
%    reduced plane:
%      protan: replace L, preserving M and S.
%      deutan: replace M, preserving L and S.
% 6. Convert the simulated LMS values back to CIE XYZ, linear RGB, and sRGB.
%
% This is a modernization of the 1999 algorithm, not a byte-for-byte
% reproduction of the paper's published replacement colourmaps. In
% particular, the paper's printed RGB-to-LMS matrix is treated as an
% implementation artifact of its assumed CRT display model, not as a
% universal conversion matrix.
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
% gamma   | true**    | Use inverse/forward sRGB gamma correction.
%         | false     | Apply the simulation directly to the <rgb> values.
% --------|-----------|----------------------------------------------------
% wpt     | 1x3 float | XYZ reference white point, scaled Y==1.
%         |           | Default is D65: [0.95047,1,1.08883].
% --------|-----------|----------------------------------------------------
% RGB2XYZ | 3x3 float | Matrix converting linear RGB to XYZ.
%         |           | Default is IEC 61966-2-1 sRGB.
% --------|-----------|----------------------------------------------------
% XYZ2LMS | 3x3 float | Matrix converting XYZ to LMS-like responses. Default
%         |           | is the Hunt-Pointer-Estevez matrix used by CIECAM02.
% --------|-----------|----------------------------------------------------
% bluXYZ  | 1x3 float | XYZ value of the display blue primary.
%         | []**      | If empty, it is derived from <RGB2XYZ> as the XYZ
%         |           | value of linear RGB [0,0,1]. Only its chromatic
%         |           | direction is relevant to the reduced plane.
% --------|-----------|----------------------------------------------------
% bluLMS  | 1x3 float | LMS value of the display blue primary.
%         | []**      | If non-empty, this overrides <bluXYZ> and is used
%         |           | directly. Only its chromatic direction is relevant.
%
%% Examples %%
%
%%% Simulate a single RGB triple for protanopia %%%
%
%   >> vienot1999([1,0,0], 'protan')
%   ans = [0.4498   0.4498   0]
%
%%% View an image as seen by a deuteranope %%%
%
%   >> I = imread("peppers.png");
%   >> imshow(vienot1999(I, 'd'))
%
%%% View Parula colormap under protanopia and deuteranopia %%%
%
%   >> I = permute(parula(17),[3,1,2]);
%   >> imshow([I;vienot1999(I,'p');vienot1999(I,'d')])
%
%%% Use a custom blue-primary direction %%%
%
%   >> opt.bluLMS = [0.0330,0.1006,1.0419];
%   >> imshow(vienot1999(I, 'deutan', opt))
%
%% Notes %%
%
% * This implementation models the protanopic and deuteranopic cases
%   defined by Vienot et al. (1999). Tritan simulation and anomalous
%   trichromacy are not implemented.
% * Unlike the reference paper, this implementation defaults to sRGB input
%   and output, rather than a specifically assumed CRT monitor and printed
%   RGB-to-LMS matrix. Users requiring a different display model may supply
%   the conversion matrices <RGB2XYZ>, <XYZ2LMS>, the whitepoint <wpt>,
%   and/or the blue-primary direction <bluLMS>.
% * The default <XYZ2LMS> matrix is the Hunt-Pointer-Estevez matrix defined
%   by CIECAM02. It is used here as a practical LMS-like cone-response
%   space, not as a full CIECAM02 appearance model.
% * The reduced plane is the plane containing the origin, white point, and
%   display blue primary. This corresponds to the KBWY diagonal plane used
%   by Vienot et al., generalized to the selected RGB/XYZ/LMS display model.
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
%          values after Vienot projection.
%   opts = StructureScalar, the used parameter values as per 'Options' above.
%
%% Dependencies %%
%
% * MATLAB R2009b or later.
% * No toolboxes are required.
%
% See also BRETTEL1997 CVDSIM DALTONIZER MACHADO2010 MILIC2015
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
	error('SC:vienot1999:rgb:NotNumeric',...
		'1st input <rgb> must be a numeric array, not %s',class(rgb))
end
assert(isreal(rgb),...
	'SC:vienot1999:rgb:NotReal',...
	'1st input <rgb> must be a real array (not complex).')
assert(isz(end)==3 || isequal(isz,[3,1]),...
	'SC:vienot1999:rgb:InvalidSize',...
	'1st input <rgb> last dimension must have size 3 (e.g. Nx3 or RxCx3).')
assert(all(0<=rgb(:)&rgb(:)<=1),'SC:vienot1999:rgb:OutOfRange',...
	'1st input <rgb> values must be 0<=rgb<=%d',mxv)
rgb = reshape(rgb,[],3);
%
typ = mSS2C(typ);
assert(ischar(typ)&&ndims(typ)==2&&size(typ,1)==1,...
	'SC:vienot1999:typ:NotText',...
	'Second input <typ> must be a character vector or a string scalar.') %#ok<ISMAT>
%
switch lower(typ)
	case {'p','protan','protanopia'}
		dst = 'p';
	case {'d','deutan','deuteranopia'}
		dst = 'd';
	case {'t','tritan','tritanopia'}
		error('SC:vienot1999:typ:TritanNotSupported',...
			'Second input <typ> "%s" is not supported: Vienot et al. (1999) defines protanopic and deuteranopic simulations only.',typ)
	case {'protanomaly','deuteranomaly','tritanomaly'}
		error('SC:vienot1999:typ:AnomalousNotSupported',...
			'Second input <typ> "%s" is not supported: Vienot et al. (1999) defines dichromatic simulations only.',typ)
	otherwise
		error('SC:vienot1999:typ:NotSupported',...
			'Second input <typ> "%s" is not supported: use "protan"/"deutan" or their dichromatic full names or initials.',typ)
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
	'bluXYZ',[],...
	'bluLMS',[]);
%
% Check any supplied option field names and values:
switch nargin
	case 2 % no user-supplied options
		% Use defaults.
	case 3 % options in a struct
		assert(isstruct(opts)&&isscalar(opts),...
			'SC:vienot1999:options:NotScalarStruct',...
			'Third input <opts> must be a scalar structure, or options must be supplied as name-value pairs.')
		opts = structfun(@mSS2C,opts,'UniformOutput',false);
		stpo = mOptions(stpo,opts);
	otherwise % options as <name-value> pairs
		tmp = cellfun(@mSS2C,[{opts},varargin],'UniformOutput',false);
		opts = cell2struct(tmp(2:2:end),tmp(1:2:end),2);
		stpo = mOptions(stpo,opts);
end
opts = stpo;
%
%% Prepare Simulation Parameters %%
%
E = stpo.wpt * stpo.XYZ2LMS.';
%
if isempty(stpo.bluLMS)
	if isempty(stpo.bluXYZ)
		bluXYZ = [0,0,1] * stpo.RGB2XYZ.';
	else
		bluXYZ = stpo.bluXYZ;
	end
	B = bluXYZ * stpo.XYZ2LMS.';
else
	B = stpo.bluLMS;
end
%
% Plane normal [a,b,c] = cross(E,B), where B is the display blue primary.
plane = [...
	E(2).*B(3) - E(3).*B(2),...
	E(3).*B(1) - E(1).*B(3),...
	E(1).*B(2) - E(2).*B(1)];
%
% Check that the blue primary and neutral axis define a usable Vienot plane.
bnm = hypot(hypot(B(1),B(2)),B(3));
pnm = hypot(hypot(plane(1),plane(2)),plane(3));
enm = hypot(hypot(E(1),E(2)),E(3));
tol = 100 * eps(max([1;enm;bnm;pnm]));
%
assert(bnm>0,...
	'SC:vienot1999:options:blue:ZeroBlue',...
	'The blue-primary direction must contain at least one non-zero value.')
assert(pnm > tol*enm.*bnm,...
	'SC:vienot1999:options:blue:DegeneratePlane',...
	'The blue-primary direction must define a non-degenerate plane with the white point.')
assert(abs(plane(1)) > tol*pnm,...
	'SC:vienot1999:options:blue:DegenerateProtanPlane',...
	'The blue-primary direction must define a protan plane with a non-zero L coefficient.')
assert(abs(plane(2)) > tol*pnm,...
	'SC:vienot1999:options:blue:DegenerateDeutanPlane',...
	'The blue-primary direction must define a deutan plane with a non-zero M coefficient.')
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
%% Apply Vienot et al. Projection %%
%
a = plane(1);
b = plane(2);
c = plane(3);
switch dst
	case 'p' % protanopia: replace L, preserve M and S.
		sim(:,1) = -(b.*lms(:,2) + c.*lms(:,3)) ./ a;
	case 'd' % deuteranopia: replace M, preserve L and S.
		sim(:,2) = -(a.*lms(:,1) + c.*lms(:,3)) ./ b;
	otherwise
		error('SC:vienot1999:typ:MissingCase','Please report this bug: missing case.')
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%vienot1999
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
		error('SC:vienot1999:options:UnknownOptionName',...
			'Unknown option: <%s>.\nOptions are:%s.',ofn,ont(2:end))
	elseif nnz(oix)>1
		dnt = sprintf(', <%s>',ofc{oix});
		error('SC:vienot1999:options:DuplicateOptionNames',...
			'Duplicate option names:%s.',dnt(2:end))
	end
	arg = opts.(ofn);
	dfn = dfc{dix};
	switch dfn
		case 'gamma'
			mLogical()
		case {'RGB2XYZ','XYZ2LMS'}
			mFullRank(false)
		case {'bluXYZ','bluLMS'}
			mMatrix(true,3)
		case 'wpt'
			mWhitePoint(false)
		otherwise
			error('SC:vienot1999:options:MissingCase','Please report this bug.')
	end
	stpo.(dfn) = arg;
end
%
%% Nested Functions %%
%
	function mLogical() % scalar logical.
		assert(isequal(arg,false)||isequal(arg,true),...
			sprintf('SC:vienot1999:%s:NotScalarLogical',dfn),...
			'The <%s> value must be true/1 or false/0.',dfn)
		arg = logical(arg);
	end
%
	function mWhitePoint(imt) % positive real 1x3 white point vector.
		mMatrix(imt,3)
		assert(all(arg(:)>0) && arg(2)==1,...
			'SC:vienot1999:wpt:InvalidValue',...
			'The <wpt> value must contain three positive numeric values and wpt(2)==1.')
		arg = reshape(arg,1,3);
	end
%
	function mFullRank(imt) % real finite full-rank 3x3 matrix.
		mMatrix(imt,[3,3])
		assert(rank(arg)==3,...
			sprintf('SC:vienot1999:%s:NotFullRank',dfn),...
			'The <%s> value must be a full-rank 3x3 matrix.',dfn)
	end
%
	function mMatrix(imt,siz) % real finite numeric array of specified size, optionally [].
		if imt && isnumeric(arg) && isempty(arg)
			arg = [];
			return
		end
		assert(isnumeric(arg),...
			sprintf('SC:vienot1999:%s:NotNumeric',dfn),...
			'The <%s> value must be numeric.',dfn)
		assert(isreal(arg),...
			sprintf('SC:vienot1999:%s:NotRealNumeric',dfn),...
			'The <%s> value cannot be complex.',dfn)
		if isscalar(siz)
			assert(numel(arg)==siz,...
				sprintf('SC:vienot1999:%s:InvalidSize',dfn),...
				'The <%s> value must contain %d numeric values.',dfn,siz)
		else
			assert(isequal(size(arg),siz),...
				sprintf('SC:vienot1999:%s:InvalidSize',dfn),...
				'The <%s> value must have size %s.',dfn,mat2str(siz))
		end
		assert(all(isfinite(arg(:))),...
			sprintf('SC:vienot1999:%s:NotFiniteNumeric',dfn),...
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
% Copyright (c) 1999 Francoise Vienot, Hans Brettel, and John D. Mollon.
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%license