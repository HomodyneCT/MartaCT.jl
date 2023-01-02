var documenterSearchIndex = {"docs":
[{"location":"api/Applicative/#Applicative","page":"Applicative","title":"Applicative","text":"","category":"section"},{"location":"api/Applicative/","page":"Applicative","title":"Applicative","text":"Modules = [MartaCT.Applicative]","category":"page"},{"location":"api/CTTestImages/#Test-images","page":"Test images","title":"Test images","text":"","category":"section"},{"location":"api/CTTestImages/","page":"Test images","title":"Test images","text":"Modules = [MartaCT.CTTestImages]","category":"page"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.ImageParams","page":"Test images","title":"MartaCT.CTTestImages.ImageParams","text":"struct ImageParams{T}\n\nHold information to construct the input test image.\n\n\n\n\n\n","category":"type"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.ImageParams-Union{Tuple{}, Tuple{Type{T}}, Tuple{T}} where T","page":"Test images","title":"MartaCT.CTTestImages.ImageParams","text":"ImageParams([T=Float32]; <keyword arguments>) where {T}\n\nConstruct ImageParams object.\n\nArguments:\n\nwidth=200: width of the rectangle gray scale image.\nheight=40: height of the rectangle gray scale image.\npad=30: padding around image.\ndist=10: distance between images.\nradius=20: radius of calibration circle.\ngray_scale=1000..1000: interval of gray scale values.\ncalibration_value=nothing: value of the calibration circle. If not specified, is gray_scale[2].\nbackground=nothing: value to be used as background.\nhounsfield=false: whether the gray scale should be in Hounsfield units.\n\n\n\n\n\n","category":"method"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.background_position","page":"Test images","title":"MartaCT.CTTestImages.background_position","text":"background_position(imp::AbstractImageParams)\n\nReturn suitable position to calibrate background.\n\n\n\n\n\n","category":"function"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.calibration_position","page":"Test images","title":"MartaCT.CTTestImages.calibration_position","text":"calibration_position(imp::AbstractImageParams)\n\nGet the position for calibration (max value).\n\n\n\n\n\n","category":"function"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.circle_image","page":"Test images","title":"MartaCT.CTTestImages.circle_image","text":"circle_image(imp::AbstractImageParams)\n\nCreate a square image with a circle of given value from parameters imp.\n\n\n\n\n\n","category":"function"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.circle_image-Union{Tuple{}, Tuple{Type{T}}, Tuple{T}} where T","page":"Test images","title":"MartaCT.CTTestImages.circle_image","text":"circle_image([T=Float32]; radius=20, value=1) where {T}\n\nCreate a square image with a circle of given value.\n\nExamples\n\njulia> circle_image(30, 0.8)\n30×30 Array{Float32,2}:\n[...]\n\nSee also: gray_scale_image, combine_images\n\n\n\n\n\n","category":"method"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.circle_polar_image-Union{Tuple{T}, Tuple{Type{T}, Integer, Integer, Integer}} where T","page":"Test images","title":"MartaCT.CTTestImages.circle_polar_image","text":"circle_polar_image([T=Float32] nr, nϕ, radius; value=1) where {T}\n\nCreate a circle of radius radius inside a nr×nϕ image in polar coordinates .\n\nExamples\n\njulia> circle_polar_image(30, 30, 10)\n30×30 Array{Float32,2}:\n[...]\n\n\n\n\n\n","category":"method"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.circle_position","page":"Test images","title":"MartaCT.CTTestImages.circle_position","text":"circle_position(imp::AbstractImageParams)\n\nReturn circle position inside image given image parameters.\n\n\n\n\n\n","category":"function"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.combine_images-Union{Tuple{T}, Tuple{ImageParams{T}, AbstractMatrix{T}, AbstractMatrix{T}}} where T","page":"Test images","title":"MartaCT.CTTestImages.combine_images","text":"combine_images(imp::ImageParams{T}, rect::AbstractMatrix{T}, circle::AbstractMatrix{T}) where {T}\n\nHelper function to create a gray scale image.\n\nCombine the rect image with the circle image for calibration. The circle image should be smaller than the gray scale image!\n\nSee also: gray_scale_image, circle_image\n\n\n\n\n\n","category":"method"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.create_image-Tuple{ImageParams}","page":"Test images","title":"MartaCT.CTTestImages.create_image","text":"create_image(par::ImageParams)\n\nCreate gray scale image.\n\n\n\n\n\n","category":"method"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.create_pyramid_image-Tuple{ImageParams}","page":"Test images","title":"MartaCT.CTTestImages.create_pyramid_image","text":"create_pyramid_image(par::ImageParams)\n\nCreate pyramid gray scale image.\n\n\n\n\n\n","category":"method"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.gray_scale_image-Union{Tuple{ImageParams{T}}, Tuple{T}} where T","page":"Test images","title":"MartaCT.CTTestImages.gray_scale_image","text":"gray_scale_image(imp::ImageParams)\n\nCreate an image with a gray scale rectangle with given scale gray_scale from parameters imp.\n\n\n\n\n\n","category":"method"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.gray_scale_image-Union{Tuple{}, Tuple{Type{T}}, Tuple{T}} where T","page":"Test images","title":"MartaCT.CTTestImages.gray_scale_image","text":"gray_scale_image(\n    [T=Float32];\n    rows=40,\n    cols=200,\n    swidth=nothing,\n    sheight=nothing,\n    gray_scale=-1000..1000,\n) where {T}\n\nCreate an image with a gray scale rectangle with given scale gray_scale.\n\nExamples\n\njulia> gray_scale_image(swidth=80, sheight=40, gray_scale=-1000..1000)\n40×80 Array{Float32,2}:\n[...]\n\nSee also: circle_image, combine_images\n\n\n\n\n\n","category":"method"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.gray_scale_indices-Tuple{ImageParams}","page":"Test images","title":"MartaCT.CTTestImages.gray_scale_indices","text":"gray_scale_indices(imp::ImageParams)\n\nReturn a tuple (row, column_range) where row is the row index of the scale and column_range is the range of columns.\n\n\n\n\n\n","category":"method"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.pyramid_gray_scale_image-Union{Tuple{ImageParams{T}}, Tuple{T}} where T","page":"Test images","title":"MartaCT.CTTestImages.pyramid_gray_scale_image","text":"pyramid_gray_scale_image(imp::ImageParams)\n\nCreate an image with a pyramid gray scale rectangle with given scale gray_scale from parameters imp.\n\nCreate an image with a pyramid gray scale rectangle with given scale gray_scale from parameters imp.\n\n\n\n\n\n","category":"method"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.pyramid_gray_scale_image-Union{Tuple{}, Tuple{Type{T}}, Tuple{T}} where T","page":"Test images","title":"MartaCT.CTTestImages.pyramid_gray_scale_image","text":"pyramidgrayscaleimage([T=Float32]; swidth=200, sheight=40, grayscale=0..1, plateau = 0) where {T}\n\nCreate an image with a pyramid gray scale rectangle with given scale grayscale. The plateau parameter denotes the percentage of the rectangle width that has the maximum value of `grayscale`.\n\nExamples\n\njulia> pyramid_gray_scale_image(swidth=80, sheight=40, gray_scale=-1000..1000)\n40×80 Array{Float32,2}:\n[...]\n\nSee also: gray_scale_image, circle_image, combine_images\n\n\n\n\n\n","category":"method"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.square_image-Union{Tuple{T}, Tuple{Type{T}, Integer, Integer}, Tuple{Type{T}, Integer, Integer, Union{Nothing, Integer}}} where T","page":"Test images","title":"MartaCT.CTTestImages.square_image","text":"square_image([T=Float32] r, c; l=nothing) where {T}\n\nCreate a l×l square inside a r×c image.\n\nExamples\n\njulia> square_image(30, 30; l=10)\n30×30 Array{Float32,2}:\n[...]\n\n\n\n\n\n","category":"method"},{"location":"api/CTTestImages/#MartaCT.CTTestImages.square_position-Tuple{SquareParams}","page":"Test images","title":"MartaCT.CTTestImages.square_position","text":"square_position(imp::SquareParams)\n\nReturn square position inside image given image parameters.\n\n\n\n\n\n","category":"method"},{"location":"api/Monads/#Monads","page":"Monads","title":"Monads","text":"","category":"section"},{"location":"api/Monads/","page":"Monads","title":"Monads","text":"Modules = [MartaCT.Monads]","category":"page"},{"location":"api/Algorithms/#Algorithms-documentation","page":"Algorithms documentation","title":"Algorithms documentation","text":"","category":"section"},{"location":"api/Algorithms/","page":"Algorithms documentation","title":"Algorithms documentation","text":"Modules = [MartaCT.AbstractAlgorithms, MartaCT.RadonAlgorithm]","category":"page"},{"location":"api/Algorithms/#MartaCT.AbstractAlgorithms.iradon","page":"Algorithms documentation","title":"MartaCT.AbstractAlgorithms.iradon","text":"iradon(sinog::AbstractMatrix, xs[, ys[, algorithm::AbstractProjectionAlgorithm[, coo::AbstractCoordinates]]]; <keyword arguments>)\n\nCompute the inverse Radon transform of sinog on the points given by the vectors xs and ys. If ys is omitted, the reconstruction is performed on the square with xs == ys. The default reconstruction algorithm is FBP. Please refer to the respective documentation for additional parameters.\n\nSee Also: reconstruct_image\n\n\n\n\n\n","category":"function"},{"location":"api/Algorithms/#MartaCT.AbstractAlgorithms.iradon-Tuple{AbstractMatrix{T} where T, AbstractParallelBeamGeometry, AbstractReconstructionAlgorithm}","page":"Algorithms documentation","title":"MartaCT.AbstractAlgorithms.iradon","text":"iradon(sinog::AbstractMatrix[, geometry::AbstractGeometry[, params::AbstractParams[, algorithm::AbstractProjectionAlgorithm]]]; <keyword arguments>)\n\nCompute the inverse Radon transform of sinog with explicit geometry. If the algorithm needs specific parameters, these can be passed with params. Additional parameters can be passed to the algorithm through keyword arguments. Please see the respective documentation for more details. The default algorithm is FBP.\n\nSee Also: reconstruct_image\n\n\n\n\n\n","category":"method"},{"location":"api/Algorithms/#MartaCT.AbstractAlgorithms.iradon-Tuple{AbstractMatrix{T} where T, AbstractReconstructionAlgorithm}","page":"Algorithms documentation","title":"MartaCT.AbstractAlgorithms.iradon","text":"iradon(sinog::AbstractMatrix[, algorithm::AbstractProjectionAlgorithm]; <keyword arguments>)\n\nCompute the inverse Radon transform of sinog with parameters given as keyword arguments. The parameters depend on the algorithm, please see the relative documentation. The default algorithm is FBP.\n\nSee Also: reconstruct_image\n\n\n\n\n\n","category":"method"},{"location":"api/Algorithms/#MartaCT.AbstractAlgorithms.radon","page":"Algorithms documentation","title":"MartaCT.AbstractAlgorithms.radon","text":"radon(image::AbstractMatrix[, algorithm::AbstractProjectionAlgorithm]; <keyword arguments>)\n\nCompute the Radon transform of image with parameters given as keyword arguments. The parameters depend on the algorithm, please see the relative documentation. The default algorithm is Radon.\n\nSee Also: project_image\n\n\n\n\n\n","category":"function"},{"location":"api/Algorithms/#MartaCT.AbstractAlgorithms.radon-2","page":"Algorithms documentation","title":"MartaCT.AbstractAlgorithms.radon","text":"radon(image::AbstractMatrix, ts, ϕs[, algorithm::AbstractProjectionAlgorithm]; <keyword arguments>)\n\nCompute the Radon transform of image with integration points given by the vector ts and projection angles given by the vector ϕs. Additional parameters can be passed to the algorithm through keyword arguments.Please see the relative documentation. The default algorithm is Radon.\n\nSee Also: project_image\n\n\n\n\n\n","category":"function"},{"location":"api/Algorithms/#MartaCT.AbstractAlgorithms.radon-3","page":"Algorithms documentation","title":"MartaCT.AbstractAlgorithms.radon","text":"radon(image::AbstractMatrix, geometry::AbstractGeometry[, algorithm::AbstractProjectionAlgorithm]; <keyword arguments>)\n\nCompute the Radon transform of image with explicit geometry. Additional parameters can be passed to the algorithm through keyword arguments. Please see the relative documentation. The default algorithm is Radon.\n\nSee Also: project_image\n\n\n\n\n\n","category":"function"},{"location":"api/Algorithms/#MartaCT.RadonAlgorithm.radon_diag","page":"Algorithms documentation","title":"MartaCT.RadonAlgorithm.radon_diag","text":"radon_diag(image::AbstractMatrix, ts::AbstractVector, ϕs::AbstractVector; <keyword arguments>)\n\nCompute the Radon transform of image inside a circle of radius hypot(rows,cols)/2 where rows and cols are the dimensions of image.\n\nSee Also: radon_square\n\n\n\n\n\n","category":"function"},{"location":"api/Algorithms/#MartaCT.RadonAlgorithm.radon_square","page":"Algorithms documentation","title":"MartaCT.RadonAlgorithm.radon_square","text":"radon_square(image::AbstractMatrix, ts::AbstractVector, ϕs::AbstractVector; <keyword arguments>)\n\nCompute the Radon transform of image inside the circle contained in the square of side min(rows,cols) where rows and cols are the dimensions of image.\n\nSee also: radon_diag\n\n\n\n\n\n","category":"function"},{"location":"api/Simulations/#Simulations","page":"Simulations","title":"Simulations","text":"","category":"section"},{"location":"api/Simulations/","page":"Simulations","title":"Simulations","text":"Modules = [MartaCT.Simulations]","category":"page"},{"location":"api/Simulations/#MartaCT.Simulations.generate_photons-Tuple{Integer, Integer, Integer}","page":"Simulations","title":"MartaCT.Simulations.generate_photons","text":"generate_photons(n::Integer, nx::Integer, nϕ::Integer)\n\nGenerate n Poisson distributed photons per projection angle over an array of nx detectors.\n\nReturn a nd × nϕ matrix with the generated photons. ```\n\n\n\n\n\n","category":"method"},{"location":"api/Simulations/#MartaCT.Simulations.sample_sinogram_external-Tuple{AbstractMatrix{T} where T}","page":"Simulations","title":"MartaCT.Simulations.sample_sinogram_external","text":"sample_sinogram_external(sinog::AbstractMatrix; <keyword arguments>)\n\nSimulate a low dose CT scan. This Resample sinog with n photons.\n\nArguments\n\nsinog: sinogram data.\nsinog_path='_tmp_sinog.dat': path to a file where to write the sinogram to.\nresampled_path='_tmp_resampled.dat': path to a file where to read the resampled sinogram from.\nn=10000: mean number of photons.\nϵ=1: detectors quantum efficiency.\ntake_log=true: whether to take the logarithm of the resampled intensities to obtain the corresponding sinogram.\nverbosity=0: set verbosity level.\nprogress=false: show progress bar.\noptions=[]: Additional options to be passed to the program.\n\n\n\n\n\n","category":"method"},{"location":"api/Simulations/#MartaCT.Simulations.simulate-Tuple{MartaCT.Simulations.CTSimulation, AbstractMatrix{T} where T}","page":"Simulations","title":"MartaCT.Simulations.simulate","text":"simulate(sim::CTSimulation, sinog::AbstractMatrix; <keyword arguments>)\n\nSimulate a CT scan by generating sim.nphotons Poisson-distributed photons for each angle, which are then detected with a probability given by p_ϕ_j(x_i) = exp(-s_ij), where s_ij is sinog[i,j].\n\nThe keyword arguments are the same as those for simulate_ct.\n\nSee also: simulate_ct\n\n\n\n\n\n","category":"method"},{"location":"api/Simulations/#MartaCT.Simulations.simulate_ct-Tuple{AbstractMatrix{T} where T}","page":"Simulations","title":"MartaCT.Simulations.simulate_ct","text":"simulate_ct(sinog::AbstractMatrix; <keyword arguments>) where {T}\n\nSimulate a low dose CT scan. This samples sinog with n random photons per projection angle.\n\nArguments\n\nsinog: sinogram data.\nnphotons::Integer=10000: mean number of photons.\nϵ::Real=1: detectors quantum efficiency.\ntake_log::Bool=false: whether to take the logarithm of the resampled intensities to obtain the corresponding sinogram.\n\n\n\n\n\n","category":"method"},{"location":"api/Simulations/#StatsBase.sample-Tuple{AbstractMatrix{T} where T, AbstractVector{T} where T}","page":"Simulations","title":"StatsBase.sample","text":"StatsBase.sample(data::AbstractMatrix, xs::AbstractVector[; nsamples=1000, nblks=1, nbins=nothing])\n\nCompute nsamples from xs according to the distribution given by each column of data. Optionally, if nblks > 1, nblks × nsamples samples are computed. The sampled data are collected into a histogram of length nbins for each column in data. The resulting data have size (nbins, size(data, 2), nblks).\n\n\n\n\n\n","category":"method"},{"location":"QuantumOptics/#Quantum-Optics-example","page":"Quantum Optics","title":"Quantum Optics example","text":"","category":"section"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"ENV[\"GKS_WSTYPE\"] = \"svg\"\nusing MartaCT, Plots, QuantumOptics, IntervalSets, LinearAlgebra, Plots.Measures\ngr()\nPlots.reset_defaults()\ndefault(size=(500,300), rightmargin=0.5cm)","category":"page"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"MartaCT can also be used for Quantum Tomography as the basic filtered backprojection algorithm can be used to reconstruct a quantum state (or more precisely the Wigner function) from the measurements of the position.","category":"page"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"First we need to import the packages for this example. We use the QuantumOptics Julia package for the definition of the quantum objects.","category":"page"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"using MartaCT, Plots, QuantumOptics, IntervalSets, LinearAlgebra","category":"page"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"At this point, let's construct an optical cat state:","category":"page"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"N = 64 # density matrix dimension\nα = 5\nbs = FockBasis(N-1)\nν = inv(√(2*(1+exp(-2abs2(α)))))\nψ = ν * (coherentstate(bs, α) + coherentstate(bs, -α))\nρ = ψ ⊗ ψ'\n@show tr(ρ)\nnothing # hide","category":"page"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"bar(diag(ρ.data)|>real; c=:red, leg=:none, lc=:red)\nsavefig(\"rho-diag.svg\"); nothing # hide","category":"page"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"(Image: )","category":"page"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"We need to get the Wigner function representation in order to compute the marginal distributions of the position.","category":"page"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"ζ = 10\nxs = linspace(-ζ..ζ, 200)\nW = wigner(ρ, xs, xs) |> permutedims\nδW = 4ζ^2 / length(W)\n@show sum(W) * δW\nheatmap(xs, xs, W)\nsavefig(\"Wcat.svg\"); nothing # hide","category":"page"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"(Image: )","category":"page"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"The marginal distributions can be computed as a Radon transform of the Wigner function:","category":"page"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"ϕs = linspace(ORI(0..2π), 800)\nmarg = radon(W, xs, ϕs, RadonSquare())\n@show sum(marg[:,1])\nheatmap(ϕs, xs, marg)\nsavefig(\"marg.svg\"); nothing # hide","category":"page"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"(Image: )","category":"page"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"Now we can employ the standard FBP algorithm to recover the Wigner distribution:","category":"page"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"Wrec = iradon(marg, xs, xs, FBPFFTSquare())\nheatmap(xs, xs, Wrec)\nsavefig(\"Wrec.svg\"); nothing # hide","category":"page"},{"location":"QuantumOptics/","page":"Quantum Optics","title":"Quantum Optics","text":"(Image: )","category":"page"},{"location":"Tutorial/#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"ENV[\"GKS_WSTYPE\"] = \"svg\"\nusing MartaCT, Plots, Plots.Measures, IntervalSets\ngr()\nPlots.reset_defaults()\nPlots.default(size=(500,300), rightmargin=1cm, color=:grays)","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"Let's see a brief example on how MartaCT can be used.","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"As a starting point we prepare the session importing some packages:","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"using MartaCT, Plots, IntervalSets","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"Now we can generate our first image from the available test images (see Test Images) for a full list of test images.","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"img = GrayScalePyramid()\npbg = ParallelBeamGeometry(img, nϕ=800)\nnothing # hide","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"We can plot the image as","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"heatmap(img)\nsavefig(\"img.svg\"); nothing # hide","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"which produces the image below.","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: )","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"The gray scale values of the image are in Hounsfield units, typically used in medical applications. In this case, where the image is just used to test the algorithm, there is, of course, no meaning in the actual values. We use this scale by default so that it can be used to compare with standard imaging procedures.","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"Now we can compute the sinogram of img using the function project_image:","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"sinog = project_image(img, RadonInfo(pbg), progress=false, rescaled=true)\nheatmap(sinog)\nsavefig(\"sinog.svg\"); nothing # hide","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: )","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"One can also use the low level function radon.","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"The reconstruction of the image can be done using the function reconstruct_image or the low level function iradon:","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"tomog = reconstruct_image(sinog, FBPInfo(pbg), progress=false)\nheatmap(tomog)\nsavefig(\"tomog.svg\"); nothing # hide","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: )","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"Here we can see that the final scale is not recovered properly: we need to calibrate the image! MartaCT provides a convenient framework for image calibration and analysis. Since we are studying a well defined image, MartaCT knows how to calibrate the resulting image according to the original image. In the general case, this would require to provide at least two ROIs (regions of interest) which are used to calibrate the final image.","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"We can just do the following","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"calibrate_tomogram!(tomog, img; window=-1000..1000)\nheatmap(tomog)\nsavefig(\"tomog_calib.svg\"); nothing # hide","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"obtaining the final result","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: )","category":"page"},{"location":"Tutorial/","page":"Tutorial","title":"Tutorial","text":"The window keyword argument can be used to apply a window filter on the image so that the final interval of values specified through window are displayed.","category":"page"},{"location":"api/Interpolation/#Interpolation","page":"Interpolation","title":"Interpolation","text":"","category":"section"},{"location":"api/Interpolation/","page":"Interpolation","title":"Interpolation","text":"Modules = [MartaCT.Interpolation]","category":"page"},{"location":"api/Geometry/#Geometry","page":"Geometry","title":"Geometry","text":"","category":"section"},{"location":"api/Geometry/","page":"Geometry","title":"Geometry","text":"Modules = [MartaCT.Geometry, MartaCT.FanBeam]","category":"page"},{"location":"api/Geometry/#MartaCT.Geometry.FanBeamGeometry-Union{Tuple{}, Tuple{Type{T}}, Tuple{T}, Tuple{Type{T}, AbstractTomograph}} where T<:Real","page":"Geometry","title":"MartaCT.Geometry.FanBeamGeometry","text":"FanBeamGeometry([T=Float32]; <keyword arguments>) where {T}\n\nCreate a new FanBeamGeometry object with given parameters.\n\nIt is assumed that the detectors are arranged on a circle arc, not flat.\n\nArguments\n\n`nϕ: number of scan angles.\nnd=nothing: number of detectors.\nrows=nothing: number of rows in the reconstructed image.\ncols=nothing: number of columns in the reconstructed image.\nwidth=nothing: alias for cols.\nheight=nothing: alias for rows.\nD: focal spot to ISO distance.\nD′=nothing: focal spot to detectors distance; if not specified   and also γ is not specified, is defined so that the total   fan angle is 1 radian.\nγ=nothing: total fan angle in degrees.\nδ=1: detectors spacing (cell size).\nα=360: scan angle in degrees.\nα₀=0: scan starting angle in degrees.\ncenter=nothing: virtual center channel, defaults to (nd-1)/2.\n\n\n\n\n\n","category":"method"},{"location":"api/Geometry/#MartaCT.Geometry.ParallelBeamGeometry-Union{Tuple{}, Tuple{Type{T}}, Tuple{T}, Tuple{Type{T}, AbstractTomograph}} where T<:Real","page":"Geometry","title":"MartaCT.Geometry.ParallelBeamGeometry","text":"ParallelBeamGeometry([T=Float32]; <keyword arguments>) where {T}\n\nConstruct geometry for the simulation.\n\nArguments\n\nnϕ: number of projections.\nnd=nothing: number of detectors. If not specified, same as nϕ.\nrows=nothing: number of rows in the reconstructed image. If not specified, same as nd.\ncols=nothing: number of columns in the reconstructed image. If not specified, same as rows.\nwidth=nothing: alias for cols.\nheight=nothing: alias for rows.\nα=360: scan angle in degrees.\nα₀=0: starting scan angle in degrees.\ncenter=nothing: virtual center channel, defaults to (nd-1)/2.\n\n\n\n\n\n","category":"method"},{"location":"api/Geometry/#MartaCT.FanBeam.fan2para-Union{Tuple{Interp}, Tuple{U}, Tuple{T}, Tuple{AbstractMatrix{T}, FanBeamGeometry{U, DefaultTomograph}}} where {T, U, Interp<:MartaCT.Interpolation.AbstractInterp2DOrNone}","page":"Geometry","title":"MartaCT.FanBeam.fan2para","text":"fan2para(\n    sinog_fan::AbstractMatrix{T},\n    fbg::FanBeamGeometry;\n    <keyword arguments>\n) where {T} -> pbg, sinog_para\n\nConvert given sinogram sinog_fan from fan beam geometry to parallel beam projections.\n\nIt is assumed that the detectors are arranged on a circle arc, not flat.\n\nThe function returns a tuple where fbg is the new geometry with fan beam projections parameters and sinog_fan is the converted sinogram.\n\nArguments\n\nsinog_para: the input sinogram.\npbg: the parallel geometry parameters.\nbackground=nothing: background value to be used, defaults to 0.\ninterpolation: interpolation strategy; it should be a function   taking a matrix as input and returning a function of the indices to   get the interpolated value. By default it is a bilinear interpolation.\n\nSee Also: para2fan\n\n\n\n\n\n","category":"method"},{"location":"api/Geometry/#MartaCT.FanBeam.para2fan-Union{Tuple{Interp}, Tuple{U}, Tuple{T}, Tuple{AbstractMatrix{T}, FanBeamGeometry{U, CT} where CT<:AbstractTomograph}} where {T, U, Interp<:MartaCT.Interpolation.AbstractInterp2DOrNone}","page":"Geometry","title":"MartaCT.FanBeam.para2fan","text":"para2fan(\n    sinog_para::AbstractMatrix{T},\n    fbg::FanBeamGeometry;\n    <keyword arguments>\n) where {T} -> fbg, sinog_fan\n\nConvert given sinogram sinog_para from parallel beam geometry to fan beam projections.\n\nThe function returns a tuple where fbg is the new geometry with fan beam projections parameters and sinog_fan is the converted sinogram.\n\nArguments\n\nsinog_para: the input sinogram.\npbg: the geometry parameters.\nD: focal spot to ISO distance.\nD′=nothing: focal spot to detectors distance; if not specified   and also γ is not specified, is defined so that the total   fan angle is 1 radian.\nγ=nothing: total fan angle in degrees.\nδ=1: detectors spacing (cell size).\nbackground=nothing: background value to be used, defaults to 0.\ninterpolation: interpolation strategy; it should be a function   taking a matrix as input and returning a function of the indices to   get the interpolated value. By default it is a bilinear interpolation.\n\nSee Also: fan2para\n\n\n\n\n\n","category":"method"},{"location":"GeneralIndex/#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"GeneralIndex/","page":"Index","title":"Index","text":"","category":"page"},{"location":"api/Calibration/#Calibration","page":"Calibration","title":"Calibration","text":"","category":"section"},{"location":"api/Calibration/","page":"Calibration","title":"Calibration","text":"Modules = [MartaCT.CalibrationBase, MartaCT.Calibration]","category":"page"},{"location":"api/Calibration/#MartaCT.CalibrationBase.calibrate_image!-Tuple{AbstractMatrix{T} where T, AbstractImageParams}","page":"Calibration","title":"MartaCT.CalibrationBase.calibrate_image!","text":"calibrate_image!(image::AbstractMatrix, imp::AbstractImageParams; interval=nothing, window=nothing)\n\nPerform calibration of image using image parameters as reference.\n\n\n\n\n\n","category":"method"},{"location":"api/Calibration/#MartaCT.CalibrationBase.calibrate_image!-Tuple{AbstractMatrix{T} where T}","page":"Calibration","title":"MartaCT.CalibrationBase.calibrate_image!","text":"calibrate_image!(image::AbstractMatrix; min_pos, max_pos, interval=0..1, window=nothing)\n\nPerform calibration of image with reference values.\n\n\n\n\n\n","category":"method"},{"location":"api/Calibration/#MartaCT.CalibrationBase.calibrate_image-Tuple{AbstractMatrix{T} where T, AbstractImageParams}","page":"Calibration","title":"MartaCT.CalibrationBase.calibrate_image","text":"calibrate_image(image::AbstractMatrix, imp::AbstractImageParams; interval=nothing, window=nothing)\n\nPerform calibration of image using image parameters as reference.\n\n\n\n\n\n","category":"method"},{"location":"api/Calibration/#MartaCT.CalibrationBase.calibrate_image-Tuple{AbstractMatrix{T} where T}","page":"Calibration","title":"MartaCT.CalibrationBase.calibrate_image","text":"calibrate_image(image::AbstractMatrix; min_pos, max_pos, interval=0..1, window=nothing)\n\nPerform calibration of image with reference values.\n\n\n\n\n\n","category":"method"},{"location":"api/Calibration/#MartaCT.CalibrationBase.calibrate_tomogram!-Tuple{AbstractMatrix{T} where T, AbstractImageParams}","page":"Calibration","title":"MartaCT.CalibrationBase.calibrate_tomogram!","text":"calibrate_tomogram!(image::AbstractMatrix, imp::AbstractImageParams; interval=nothing, window=nothing)\n\nPerform calibration of reconstructed image with image parameters as reference.\n\n\n\n\n\n","category":"method"},{"location":"api/Calibration/#MartaCT.CalibrationBase.calibrate_tomogram!-Tuple{AbstractMatrix{T} where T}","page":"Calibration","title":"MartaCT.CalibrationBase.calibrate_tomogram!","text":"calibrate_tomogram(image::AbstractMatrix; min_pos, max_pos, interval=0..1, window=nothing)\n\nPerform calibration of reconstructed image with reference values.\n\n\n\n\n\n","category":"method"},{"location":"api/Calibration/#MartaCT.CalibrationBase.calibrate_tomogram-Tuple{AbstractMatrix{T} where T, AbstractImageParams}","page":"Calibration","title":"MartaCT.CalibrationBase.calibrate_tomogram","text":"calibrate_tomogram(image::AbstractMatrix, imp::AbstractImageParams; interval=nothing, window=nothing)\n\nPerform calibration of reconstructed image with image parameters as reference.\n\n\n\n\n\n","category":"method"},{"location":"api/Calibration/#MartaCT.CalibrationBase.calibrate_tomogram-Tuple{AbstractMatrix{T} where T}","page":"Calibration","title":"MartaCT.CalibrationBase.calibrate_tomogram","text":"calibrate_tomogram(image::AbstractMatrix; min_pos, max_pos, interval=0..1, window=nothing)\n\nPerform calibration of reconstructed image with reference values.\n\n\n\n\n\n","category":"method"},{"location":"api/CTImages/#CTImages","page":"CTImages","title":"CTImages","text":"","category":"section"},{"location":"api/CTImages/","page":"CTImages","title":"CTImages","text":"Modules = [MartaCT.CTImages]","category":"page"},{"location":"api/CTImages/#MartaCT.CTImages.rescale!-Union{Tuple{AbstractArray{T, N} where N}, Tuple{T}} where T","page":"CTImages","title":"MartaCT.CTImages.rescale!","text":"rescale!(x::AbstractArray;\n    interval=nothing, calibration=nothing, window=nothing)\n\nRescale array x to the interval specified by interval.\n\nIf calibration is not nothing, then rescaling is done with reference to the values given by calibration. In other words, minimum and maximum are assumed to be the values specified by calibration.\n\nSee also: rescale\n\n\n\n\n\n","category":"method"},{"location":"api/CTImages/#MartaCT.CTImages.rescale!-Union{Tuple{T}, Tuple{AbstractArray{T, N} where N, Number, Number}} where T","page":"CTImages","title":"MartaCT.CTImages.rescale!","text":"rescale!(x::AbstractArray, slope::Number, intercept::Number)\n\nIn place linear rescaling of x.\n\nSee also: rescale\n\n\n\n\n\n","category":"method"},{"location":"api/CTImages/#MartaCT.CTImages.rescale-Union{Tuple{AbstractArray{T, N} where N}, Tuple{T}} where T","page":"CTImages","title":"MartaCT.CTImages.rescale","text":"rescale(x::AbstractArray;\n    interval=nothing, calibration=nothing, window=nothing)\n\nRescale array x to the interval specified by interval.\n\nIf calibration is not nothing, then rescaling is done with reference to the values given by calibration. In other words, minimum and maximum are assumed to be the values specified by calibration.\n\nSee also: rescale!\n\n\n\n\n\n","category":"method"},{"location":"api/CTImages/#MartaCT.CTImages.rescale-Union{Tuple{T}, Tuple{AbstractArray{T, N} where N, Number, Number}} where T","page":"CTImages","title":"MartaCT.CTImages.rescale","text":"rescale(x::AbstractArray, slope::Number, intercept::Number; window)\n\nLinear rescaling of x as x * slope + intercept.\n\nSee also: rescale!\n\n\n\n\n\n","category":"method"},{"location":"api/CTImages/#MartaCT.CTImages.rotate-Union{Tuple{Interp}, Tuple{T}, Tuple{AbstractMatrix{T}, Real}} where {T, Interp<:MartaCT.Interpolation.AbstractInterp2DOrNone}","page":"CTImages","title":"MartaCT.CTImages.rotate","text":"rotate(mat::AbstractMatrix, α::Real; <keyword arguments>)\n\nRotate matrix mat about the center of angle α given in degrees. If rows and cols are not given, the rotated matrix has the same dimensions of the original matrix.\n\nArguments\n\nmat: matrix to rotate.\nα: angle in degrees.\nrows=nothing: number of rows of the rotated matrix.\ncols=nothing: number of columns of the rotated matrix.\ninterpolation: interpolation strategy. By default is   BilinearInterpolation.\n\n\n\n\n\n","category":"method"},{"location":"api/Utils/#Utilities-for-MartaCT","page":"Utilities for MartaCT","title":"Utilities for MartaCT","text":"","category":"section"},{"location":"api/Utils/","page":"Utilities for MartaCT","title":"Utilities for MartaCT","text":"Modules = [MartaCT.Utils]","category":"page"},{"location":"#Marta-CT","page":"Home","title":"Marta CT","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CT images reconstruction and analysis","category":"page"},{"location":"","page":"Home","title":"Home","text":"A package for image reconstruction. The main field of application for this package is medical imaging, but it can also serve as a base for Quantum Tomography applications.","category":"page"},{"location":"#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n  \"Installation.md\",\n  \"Tutorial.md\",\n  \"QuantumOptics.md\",\n  \"GeneralIndex.md\",\n]","category":"page"},{"location":"Installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"Installation/","page":"Installation","title":"Installation","text":"MartaCT can be installed with the following commands:","category":"page"},{"location":"Installation/","page":"Installation","title":"Installation","text":"using Pkg\npkg\"add MartaCT\"","category":"page"}]
}
