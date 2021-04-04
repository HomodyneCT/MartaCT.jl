# Tutorial

```@setup ex_img
ENV["GKS_WSTYPE"] = "svg"
using Marta, Plots, Plots.Measures, IntervalSets
gr()
Plots.reset_defaults()
Plots.default(size=(500,300), rightmargin=1cm, color=:grays)
```

Let's see a brief example on how Marta can be used.

As a starting point we prepare the session importing some
packages:

```julia
using Marta, Plots, IntervalSets
```

Now we can generate our first image from the available test
images (see [Test Images](@ref)) for a full list of test
images.

```@example ex_img
img = GrayScalePyramid()
pbg = ParallelBeamGeometry(img, nϕ=800)
nothing # hide
```

We can plot the image as

```@example ex_img
heatmap(img)
savefig("img.svg"); nothing # hide
```

which produces the image below.

![](img.svg)

The gray scale values of the image are in Hounsfield units,
typically used in medical applications. In this case, where
the image is just used to test the algorithm, there is, of
course, no meaning in the actual values. We use this scale
by default so that it can be used to compare with standard
imaging procedures.

Now we can compute the sinogram of `img` using the function
[`project_image`](@ref):

```@example ex_img
sinog = project_image(img, RadonInfo(pbg), progress=false, rescaled=true)
heatmap(sinog)
savefig("sinog.svg"); nothing # hide
```

![](sinog.svg)

One can also use the low level function [`radon`](@ref).

The reconstruction of the image can be done using the
function [`reconstruct_image`](@ref) or the low level
function [`iradon`](@ref):

```@example ex_img
tomog = reconstruct_image(sinog, FBPInfo(pbg), progress=false)
heatmap(tomog)
savefig("tomog.svg"); nothing # hide
```

![](tomog.svg)

Here we can see that the final scale is not recovered
properly: we need to calibrate the image! Marta provides a
convenient framework for image calibration and analysis.
Since we are studying a well defined image, Marta knows how
to calibrate the resulting image according to the original
image. In the general case, this would require to provide at
least two ROIs (regions of interest) which are used to
calibrate the final image.

We can just do the following

```@example ex_img
calibrate_tomogram!(tomog, img; window=-1000..1000)
heatmap(tomog)
savefig("tomog_calib.svg"); nothing # hide
```

obtaining the final result

![](tomog_calib.svg)

The `window` keyword argument can be used to apply a window
filter on the image so that the final interval of values
specified through `window` are displayed.
