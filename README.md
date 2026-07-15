# Color Vision Deficiency Simulation and Daltonization Algorithms for MATLAB #

[![View ColorVisionDeficiency on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/184209)
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=DrosteEffect/ColorVisionDeficiency)

MATLAB functions for simulating dichromatic color vision deficiency (CVD/colorblindness) and enhancing images to improve visual contrast for dichromatic CVD observers (daltonization/recoloring).

## Overview ##

This repository contains:

| Function | Description | Algorithm |
| --- | --- | --- |
| `brettel1997` | Simulate how RGB colors, colormaps, or images may appear to observers with dichromatic color vision deficiency. | Brettel, Viénot and Mollon (1997) |
| `cvdsim` | Simulate how RGB colors, colormaps, or images may appear to observers with dichromatic color vision deficiency. | Machado, Oliveira and Fernandes (2009) |
| `daltonizer` | Enhance RGB images using a fixed error-redistribution enhancement. Because it is pointwise, it also works with RGB triples and colormaps. Calls `cvdsim`. | Fidaner, Lin & Ozguven (2005) |
| `machado2010` | Enhance RGB images using a projection-based contrast enhancement. | Machado and Oliveira (2010) |
| `milic2015` | Enhance RGB images using a content-dependent naturalness-preserving enhancement. | Milić, Hoffmann, Tómács, Novaković and Milosavljević (2015) |

The enhancement functions all have the same broad purpose: they modify RGB images to make information more distinguishable to observers with dichromatic color vision deficiency. They differ mainly in how much image context they use and what trade-off they make between visibility, naturalness, simplicity, and repeatability.

The simulation functions perform the complementary task: they simulate for a normal trichromat observer how RGB colors may be perceived by a dichromatic CVD observer. This is useful when checking colors, figures, colormaps, plots, diagrams, GUIs, and images for colorblind accessibility.

The functions are self-contained MATLAB code and do not require any toolboxes.

---

## Which Function Should I Use? ##

I recommend starting with `cvdsim` and `daltonizer`.

Use `brettel1997` when you want to *check* how colors may look to someone with dichromatic CVD. This is a classic algorithm which has been modified to work with sRGB images and colormaps.

Use `cvdsim` when you want to *check* how colors may look to someone with dichromatic CVD. This is a simple, fast, and easy function to apply to images, RGB values, and colormaps.

Use `daltonizer` when you want a simple, fast image recoloring method. It is the easiest function to apply to images, RGB values, and MATLAB colormaps, but because it is fixed and content-independent it can sometimes overcorrect.

Use `machado2010` when you want an image-dependent recoloring method based on local color-contrast loss. It can preserve temporal coherence across image sequences when previous-state outputs are reused (i.e. it can be used for processing video data).

Use `milic2015` when you want an image-dependent recoloring method that aims to preserve naturalness by segmenting image chromaticities before recoloring. It exposes several options because the underlying paper leaves some implementation choices to the user.

---

## What These Functions Can Do ##

These functions support the main types of dichromatic color vision deficiency:

- L-cone deficiency: protan / protanomaly / protanopia
- M-cone deficiency: deutan / deuteranomaly / deuteranopia
- S-cone deficiency: tritan / tritanomaly / tritanopia

The functions accept common MATLAB numeric RGB data. Floating-point inputs use the range `0..1`; integer inputs use up to intmax. The main output preserves the input class.

`brettel1997`, `cvdsim` and `daltonizer` work with both `Nx3` colormaps and `RxCx3` images. `machado2010` and `milic2015` work with `RxCx3` images only, because their algorithms depend on spatial image structure.

---

## What These Functions Cannot Do ##

These functions are not a universal solution to color vision deficiency. Enhancement is always a compromise: it can improve some distinctions while changing the appearance of the image, reducing naturalness, or introducing unwanted color shifts.

`cvdsim` provides a model-based simulation, not a guarantee of how any particular person will perceive an image. Real perception varies between observers, displays, viewing conditions, and adaptation states.

`daltonizer` is intentionally simple and content-independent. It can be very convenient, especially for colormaps, but it does not know whether any particular image content actually needs strong correction.

`machado2010` and `milic2015` are intended for image recoloring, not arbitrary colormaps. They also target dichromatic recoloring rather than providing a continuous anomalous-trichromacy severity parameter.

For tritan simulation, `cvdsim` follows the Machado et al. matrices, which approximate tritanomaly; true complete S-cone loss is not explicitly modeled by those source matrices.

---

## Basic Syntax ##

```matlab
cvd = brettel1997(rgb,typ)
cvd = brettel1997(rgb,typ,options)
cvd = brettel1997(rgb,typ,'name',value,...)

cvd = cvdsim(rgb,typ)
cvd = cvdsim(rgb,typ,sev)
cvd = cvdsim(rgb,typ,sev,gamma)

rgb = daltonizer(rgb,typ)
rgb = daltonizer(rgb,typ,sev)

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

Simulate a single RGB triple for protanopia, using both `brettel1997` and `cvdsim`:

```matlab
rgb = brettel1997([1,0,0], 'protan')
rgb = cvdsim([1,0,0], 'protan')
```

Simulate an image as seen with moderate deuteranomaly:

```matlab
I = imread("peppers.png");
imshow(cvdsim(I, 'deutan', 0.6))
```

Simulate a MATLAB colormap under the three CVD types:

```matlab
I = permute(parula(17),[3,1,2]);
imshow([I; cvdsim(I,'p'); cvdsim(I,'d'); cvdsim(I,'t')])
```

Enhance an image with the fixed daltonization method:

```matlab
I = imread("peppers.png");
imshow(daltonizer(I, 'deutan', 0.6))
```

Enhance an image with `machado2010`:

```matlab
I = imread("peppers.png");
imshow(machado2010(I,'deutan'))
```

Enhance an image with `milic2015`:

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

- Brettel, H., Viénot, F., & Mollon, J. D. (1997).
  "Computerized simulation of color appearance for dichromats", Journal of the Optical Society of America A, 14(10), 2647-2655.
  <https://vision.psychol.cam.ac.uk/jdmollon/papers/Dichromatsimulation.pdf>

- Machado G. M., Oliveira M. M., Fernandes L. A. F.
  "A Physiologically-based Model for Simulation of Color Vision Deficiency", IEEE TVCG 15(6):1291–1298, 2009.
  <https://doi.org/10.1109/TVCG.2009.113>
  <https://www.inf.ufrgs.br/~oliveira/pubs_files/CVD_Simulation/CVD_Simulation.html>

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