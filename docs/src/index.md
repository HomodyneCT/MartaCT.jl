# Marta

_CT images reconstruction and analysis_

A package for image reconstruction. The main field of
application for this package is medical imaging, but it can
also serve as a base for Quantum Tomography applications.

## Contents

```@contents
```

## Getting started

### Installation

Marta can be installed with the following commands:

```julia
import Pkg
Pkg.add(url="https://gitlab.com/homodyne-ct/Marta.git")
```

## Tutorial

Let's see a brief example on how Marta can be used.

As a starting point we prepare the session importing some
packages:

```@example ex_img
using Marta, Plots, Plots.Measures
nothing # hide
```

```@setup ex_img
ENV["GKS_WSTYPE"] = "svg"
using Plots.Measures
gr()
Plots.reset_defaults()
Plots.default(size=(500,300))
```

Now we can generate our first image from the available test
images (see [`Test Images`](@ref)) for a full list of test
images.

```@example ex_img
img = GrayScalePyramid()
pbg = ParallelBeamGeometry(img, nϕ=800)
nothing # hide
```

We can plot the image as

```@example ex_img
plot(img, rightmargin=0.5cm)
savefig("img.svg"); nothing # hide
```

![](img.svg)

The gray scale values of the image are in Hounsfield units,
typically used in medical applications.

Now we can compute the sinogram of `img` using the function
[`project_image`](@ref):

```@example ex_img
sinog = project_image(img, RadonInfo(pbg), progress=false, rescaled=true)
plot(sinog, rightmargin=1cm)
savefig("sinog.svg"); nothing # hide
```

![](sinog.svg)

One can also use the low level function [`radon`](@ref).

The reconstruction of the image can be done using the
function [`reconstruct_image`](@ref) or the low level
function [`iradon`](@ref):

```@example ex_img
tomog = reconstruct_image(sinog, FBPInfo(pbg), progress=false)
plot(tomog, rightmargin=1cm)
savefig("tomog.svg"); nothing # hide
```

![](tomog.svg)

## Index

```@index
```
