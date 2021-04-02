# Marta

_CT images reconstruction and analysis_

A package for image reconstruction. The main field of application for this package is medical imaging, but it can also serve as a base for Quantum Tomography applications.

## Getting started

### Installation

Marta can be installed with the following commands:

```julia
import Pkg
Pkg.add(url="https://gitlab.com/homodyne-ct/Marta.git")
```

### Tutorial

Let's see a brief example on how Marta can be used.

As a starting point we prepare the session importing some
packages:

```@example ex_img
using Marta, Plots
nothing # hide
```

```@setup ex_img
using Plots.Measures
gr()
Plots.reset_defaults()
Plots.default(size=(400,200), margin=1cm)
```

Now we can generate our first image from the available test
images (see [`Test Images`](@ref)) for a full list of test
images.

```@example ex_img
img = GrayScalePyramid()
pbg = ParallelBeamGeometry(img)
nothing # hide
```

We can plot the image as

```@example ex_img
plot(img)
savefig("img.svg"); nothing # hide
```

![](img.svg)

Now we can compute the sinogram of `img` using the function
[`project_image`](@ref):

```@example ex_img
sinog = project_image(img, RadonInfo(pbg), progress=false, rescaled=true)
plot(sinog)
savefig("sinog.svg"); nothing # hide
```

![](sinog.svg)

One can also use the low level function [`radon`](@ref).

The reconstruction of the image can be done using the
function [`reconstruct_image`](@ref) or the low level
function [`iradon`](@ref):

```@example ex_img
tomog = reconstruct_image(sinog, FBPInfo(pbg), progress=false)
plot(tomog)
savefig("tomog.svg"); nothing # hide
```

![](tomog.svg)
