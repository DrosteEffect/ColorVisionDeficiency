# ColorVisionDeficiency #

MATLAB functions for simulating color vision deficiency (CVD/colorblindness) and recoloring images (daltonize) to improve visual contrast for CVD observers.

## Overview ##

This repository contains:

| Function | Description | Algorithm |
| --- | --- | --- |
| `cvdsim` | Simulate how RGB colors, colormaps, or images may appear to observers with protan, deutan, or tritan color vision deficiency | Machado, Oliveira and Fernandes (2009) |
| `daltonize` | Recolor RGB images using a fixed error-redistribution method. Because it is pointwise, it also works with RGB triples and colormaps | Fidaner, Lin & Ozguven (2005) |
| `machado2010` | Recolor RGB images using a projection-based contrast enhancement method | Machado and Oliveira (2010) |
| `milic2015` | Recolor RGB images using a content-dependent naturalness-preserving daltonization method | Milić, Hoffmann, Tómács, Novaković and Milosavljević (2015) |

The three recoloring functions all have the same broad purpose: they modify RGB images to make information more distinguishable to observers with color vision deficiency. They differ mainly in how much image context they use and what trade-off they make between visibility, naturalness, simplicity, and repeatability.

`cvdsim` does the complementary task: it does not improve an image, but simulates how colors may be perceived by a CVD observer. This is useful when checking figures, colormaps, plots, diagrams, GUIs, and images for colorblind accessibility.

The functions are self-contained MATLAB code and do not require any toolboxes.

---

## Which Function Should I Use? ##

Use `cvdsim` when you want to *check* how colors may look to someone with CVD.

Use `daltonize` when you want a simple, fast image recoloring method. It is the easiest option to apply to images, RGB values, and MATLAB colormaps, but because it is fixed and content-independent it can sometimes overcorrect.

Use `machado2010` when you want an image-dependent recoloring method based on local color-contrast loss. It is intended for dichromats and can preserve temporal coherence across image sequences when previous-state outputs are reused.

Use `milic2015` when you want an image-dependent recoloring method that aims to preserve naturalness by segmenting image chromaticities before recoloring. It exposes several options because the underlying paper leaves some implementation choices to the user.

---

## What These Functions Can Do ##

These functions support the main families of color vision deficiency:

- L-cone deficiency: protan / protanomaly / protanopia
- M-cone deficiency: deutan / deuteranomaly / deuteranopia
- S-cone deficiency: tritan / tritanomaly / tritanopia

The functions accept common MATLAB numeric RGB data. Floating-point inputs use the range `0..1`; integer inputs use up to intmax. The main output preserves the input class.

`cvdsim` and `daltonize` work with both `Nx3` colormaps and `RxCx3` images. `machado2010` and `milic2015` work with `RxCx3` images only, because their algorithms depend on spatial image structure.

---

## What These Functions Cannot Do ##

These functions are not a universal solution to color vision deficiency. Recoloring is always a compromise: it can improve some distinctions while changing the appearance of the image, reducing naturalness, or introducing unwanted color shifts.

`cvdsim` provides a model-based simulation, not a guarantee of how any particular person will perceive an image. Real perception varies between observers, displays, viewing conditions, and adaptation states.

`daltonize` is intentionally simple and content-independent. It can be very convenient, especially for colormaps, but it does not know whether an image actually needs strong correction.

`machado2010` and `milic2015` are intended for image recoloring, not arbitrary colormaps. They also target dichromatic recoloring rather than providing a continuous anomalous-trichromacy severity parameter.

For tritan simulation, `cvdsim` follows the Machado et al. matrices, which approximate tritanomaly; true complete S-cone loss is not explicitly modeled by those source matrices.

---

## Basic Syntax ##

```matlab
cvd = cvdsim(rgb,typ)
cvd = cvdsim(rgb,typ,sev)
cvd = cvdsim(rgb,typ,sev,gamma)

rgb = daltonize(rgb,typ)
rgb = daltonize(rgb,typ,sev)

rgb = machado2010(rgb,typ)
rgb = machado2010(rgb,typ,exaggerate)

rgb = milic2015(rgb,typ)
rgb = milic2015(rgb,typ,options)
rgb = milic2015(rgb,typ,'name',value,...)
```

Accepted type names include:

```matlab
'p', 'protan', 'protanomaly', 'protanopia'
'd', 'deutan', 'deuteranomaly', 'deuteranopia'
't', 'tritan', 'tritanomaly', 'tritanopia'
```

See the help text in each M-file for the full calling syntax, optional outputs, and implementation notes.

---

## Examples ##

Simulate a single RGB triple for protanopia:

```matlab
rgb = cvdsim([1,0,0], 'protan')
```

Simulate an image as seen with moderate deuteranomaly:

```matlab
I = imread("peppers.png");
imshow(cvdsim(I, 'deutan', 0.6))
```

Simulate a MATLAB colormap under the three CVD families:

```matlab
I = permute(parula(17),[3,1,2]);
imshow([I; cvdsim(I,'p'); cvdsim(I,'d'); cvdsim(I,'t')])
```

Recolor an image with the fixed daltonization method:

```matlab
I = imread("peppers.png");
imshow(daltonize(I, 'deutan', 0.6))
```

Recolor an image with `machado2010`:

```matlab
I = imread("peppers.png");
imshow(machado2010(I,'deutan'))
```

Recolor an image with `milic2015`:

```matlab
I = imread("peppers.png");
imshow(milic2015(I,'protan'))
```

Use custom options with `milic2015`:

```matlab
I = imread("peppers.png");
imshow(milic2015(I,'deutan', 'nseg',4, 'ang',4))
```

---

## References ##

- Machado G. M., Oliveira M. M., Fernandes L. A. F.
  "A Physiologically-based Model for Simulation of Color Vision Deficiency", IEEE TVCG 15(6):1291–1298, 2009.
  <https://doi.org/10.1109/TVCG.2009.113>

- Fidaner O., Lin P., Ozguven N.
  "Analysis of Color Blindness", 2005.
  <https://acorn.stanford.edu/psych221/projects/2005/ofidaner/project_report.pdf>

- Machado G. M., Oliveira M. M.
  "Real-Time Temporal-Coherent Color Contrast Enhancement for Dichromats", Computer Graphics Forum 29(3):933–942, 2010.
  <https://www.inf.ufrgs.br/~oliveira/pubs_files/CVD_PCA/Machado_Oliveira_EuroVis2010.pdf>

- Milić N., Hoffmann M., Tómács T., Novaković D., Milosavljević B.
  "A Content-Dependent Naturalness-Preserving Daltonization Method for Dichromatic and Anomalous Trichromatic Color Vision Deficiencies", Journal of Imaging Science and Technology 59(1):010504, 2015.
  <https://www.researchgate.net/publication/276455941_A_Content-Dependent_Naturalness-Preserving_Daltonization_Method_for_Dichromatic_and_Anomalous_Trichromatic_Color_Vision_Deficiencies>

---