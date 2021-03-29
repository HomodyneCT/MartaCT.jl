var documenterSearchIndex = {"docs":
[{"location":"#Marta","page":"Marta","title":"Marta","text":"","category":"section"},{"location":"","page":"Marta","title":"Marta","text":"Documentation for Marta.","category":"page"},{"location":"","page":"Marta","title":"Marta","text":"Modules = [Marta.CTImages, Marta.CTIO, Marta.CTScan, Marta.CTPlots]","category":"page"},{"location":"#Marta.CTImages.rescale!-Union{Tuple{AbstractArray{T,N} where N}, Tuple{W}, Tuple{U}, Tuple{T}} where W<:Number where U<:Number where T<:Number","page":"Marta","title":"Marta.CTImages.rescale!","text":"rescale!(x::AbstractArray;\n    interval=nothing, calibration=nothing, window=nothing)\n\nRescale array x to the interval specified by interval.\n\nIf calibration is not nothing, then rescaling is done with reference to the values given by calibration. In other words, minimum and maximum are assumed to be the values specified by calibration.\n\nSee also: rescale\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTImages.rescale!-Union{Tuple{U}, Tuple{T}, Tuple{AbstractArray{T,N} where N,Number,Number}} where U<:Number where T<:Number","page":"Marta","title":"Marta.CTImages.rescale!","text":"rescale!(x::AbstractArray, slope::Number, intercept::Number)\n\nIn place linear rescaling of x.\n\nSee also: rescale\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTImages.rescale-Tuple{AbstractArray}","page":"Marta","title":"Marta.CTImages.rescale","text":"rescale(x::AbstractArray;\n    interval=nothing, calibration=nothing, window=nothing)\n\nRescale array x to the interval specified by interval.\n\nIf calibration is not nothing, then rescaling is done with reference to the values given by calibration. In other words, minimum and maximum are assumed to be the values specified by calibration.\n\nSee also: rescale!\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTImages.rescale-Union{Tuple{U}, Tuple{T}, Tuple{AbstractArray{T,N} where N,Number,Number}} where U<:Number where T<:Number","page":"Marta","title":"Marta.CTImages.rescale","text":"rescale(x::AbstractArray, slope::Number, intercept::Number; window)\n\nLinear rescaling of x as x * slope + intercept.\n\nSee also: rescale!\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTImages.rotate-Union{Tuple{Interp}, Tuple{T}, Tuple{AbstractArray{T,2},Real}} where Interp<:Union{Marta.Interpolation.NoInterpolation, Marta.Interpolation.AbstractInterpolation2D} where T<:Number","page":"Marta","title":"Marta.CTImages.rotate","text":"rotate(mat::AbstractMatrix, α::Real; <keyword arguments>)\n\nRotate matrix mat about the center of angle α given in degrees. If rows and cols are not given, the rotated matrix has the same dimensions of the original matrix.\n\nArguments\n\nmat: matrix to rotate.\nα: angle in degrees.\nrows=nothing: number of rows of the rotated matrix.\ncols=nothing: number of columns of the rotated matrix.\ninterpolation: interpolation strategy. By default is   BilinearInterpolation.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.load_image-Union{Tuple{AbstractString}, Tuple{T}, Tuple{AbstractString,Type{T}}} where T<:Real","page":"Marta","title":"Marta.CTIO.load_image","text":"load_image(f::AbstractString[, T=Float32]; <keyword arguments>) where {T<:Real}\n\nLoad image from file f.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.load_image-Union{Tuple{IO}, Tuple{T}, Tuple{IO,Type{T}}} where T<:Real","page":"Marta","title":"Marta.CTIO.load_image","text":"load_image(io::IO[, T=Float32]; <keyword arguments>) where {T<:Real}\n\nLoad image from stream io.\n\nArguments\n\nio: the stream where to load the image from.\nheader=false: specify if the file contains the header.\nnϕ=0: number of projections (rows in the file for row major ordering).\nnd=0: number of detectors (columns in the file for row major ordering).\nrow_major=false: tells if the data is stored in row major ordering or not.\nheader_type=Int64: the type of the header data.\n\nSee Also: write_image, read_ct_image, write_ct_image\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.load_sinogram-Union{Tuple{AbstractString}, Tuple{T}, Tuple{AbstractString,Type{T}}} where T<:Real","page":"Marta","title":"Marta.CTIO.load_sinogram","text":"load_sinogram(f::AbstractString[, T=Float32]; <keyword arguments>) where {T<:Real}\n\nLoad sinogram from file f.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.load_sinogram-Union{Tuple{IO}, Tuple{T}, Tuple{IO,Type{T}}} where T<:Real","page":"Marta","title":"Marta.CTIO.load_sinogram","text":"load_sinogram(io::IO[, ::Type{T}=Float32]; <keyword arguments>) where {T<:Real}\n\nLoad sinogram from stream io.\n\nThis function is specifically designed to load a sinogram. Do not expect consistent result with other type of images. Use more general function instead, refer to the See Also section.\n\nArguments\n\nio: the stream where to load the image from.\nheader=false: specify if the file contains the header.\nnϕ=0: number of projections (rows in the file for row major ordering).\nnd=0: number of detectors (columns in the file for row major ordering).\nrow_major=false: tells if the data is stored in row major ordering or not.\nheader_type=Int64: the type of the header data.\n\nSee Also: write_sinogram, read_ct_image, write_ct_image\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.load_tomogram-Union{Tuple{AbstractString}, Tuple{T}, Tuple{AbstractString,Type{T}}} where T<:Real","page":"Marta","title":"Marta.CTIO.load_tomogram","text":"load_tomogram(f::AbstractString[, T=Float32]; <keyword arguments>) where {T<:Real}\n\nLoad reconstructed image from file f.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.load_tomogram-Union{Tuple{IO}, Tuple{T}, Tuple{IO,Type{T}}} where T<:Real","page":"Marta","title":"Marta.CTIO.load_tomogram","text":"load_tomogram(io::IO[, T=Float32]; <keyword arguments>) where {T<:Real}\n\nLoad reconstructed image from stream io.\n\nThis function is specifically designed to load a tomogram. Do not expect consistent result with other type of images. Use more general function instead, refer to the See Also section.\n\nArguments\n\nio: the stream where to load the image from.\nheader=false: specify if the file contains the header.\nnϕ=0: number of projections (rows in the file for row major ordering).\nnd=0: number of detectors (columns in the file for row major ordering).\nrow_major=false: tells if the data is stored in row major ordering or not.\nheader_type=Int64: the type of the header data.\n\nSee Also: write_tomogram, read_ct_image, write_ct_image\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.read_ct_image-Union{Tuple{AbstractString}, Tuple{T}, Tuple{AbstractString,Type{T}}} where T<:Real","page":"Marta","title":"Marta.CTIO.read_ct_image","text":"read_ct_image(f::AbstractString[, T=Float32]; <keyword arguments>) where {T<:Real}\n\nRead CT image from file f.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.read_ct_image-Union{Tuple{IO}, Tuple{T}, Tuple{IO,Type{T}}} where T<:Real","page":"Marta","title":"Marta.CTIO.read_ct_image","text":"read_ct_image(io::IO[, T=Float32]; kwargs...) where {T<:Real}\n\nRead CT image from the stream io.\n\nKeyword arguments\n\nheader: read rows and columns from header in the file. [default: false]\nrows: the image rows.\ncols: the image columns.\b\nheader_type: type for the two values of the header. [default: Int64]\nrow_major: assume the image is stored in row major ordering.\n\nSee also: write_ct_image\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.write_ct_image-Union{Tuple{T}, Tuple{AbstractString,AbstractArray{T,2}}} where T<:Real","page":"Marta","title":"Marta.CTIO.write_ct_image","text":"write_ct_image(f::AbstractString, image::AbstractMatrix{T}; <keyword arguments>) where {T<:Real}\n\nWrite CT image to file f.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.write_ct_image-Union{Tuple{T}, Tuple{IO,AbstractArray{T,2}}} where T<:Real","page":"Marta","title":"Marta.CTIO.write_ct_image","text":"write_ct_image(io::IO, image::AbstractMatrix{T}; kwargs...) where {T<:Real}\n\nWrite CT image to the stream io.\n\nKeyword arguments\n\nheader: if true, write header information with rows and columns. [default: false]\nheader_type: type for the two values of the header. [default: Int64]\nrow_major: write image in row major ordering.\n\nSee also: read_ct_image\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.write_image-Union{Tuple{T}, Tuple{AbstractString,AbstractArray{T,2}}} where T<:Real","page":"Marta","title":"Marta.CTIO.write_image","text":"write_image(f::AbstractString, image::AbstractMatrix{T}; <keyword arguments>) where {T<:Real}\n\nWrite input image to file f.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.write_image-Union{Tuple{T}, Tuple{IO,AbstractArray{T,2}}} where T<:Real","page":"Marta","title":"Marta.CTIO.write_image","text":"write_image(io::IO, image::AbstractMatrix{T}; <keyword arguments>) where {T<:Real}\n\nWrite input image to stream io.\n\nArguments\n\nio: the stream where the image should be written to.\nimage: the image data matrix.\nheader: specify if the image header should be written.\nrow_major: specify if the image should be stored in row major order.\nheader_type=Int64: the type of the header data.\n\nSee Also: read_image, write_ct_image, read_ct_image.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.write_sinogram-Union{Tuple{T}, Tuple{AbstractString,AbstractArray{T,2}}} where T<:Real","page":"Marta","title":"Marta.CTIO.write_sinogram","text":"write_sinogram(f::AbstractString, sinog::AbstractMatrix{T}) where {T<:Real}\n\nWrite sinogram to file f.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.write_sinogram-Union{Tuple{T}, Tuple{IO,AbstractArray{T,2}}} where T<:Real","page":"Marta","title":"Marta.CTIO.write_sinogram","text":"write_sinogram(io::IO, sinog::AbstractMatrix{T}; <keyword arguments>) where {T<:Real}\n\nWrite sinogram sinog to stream io.\n\nThis function is specifically designed to deal with a sinogram image.\n\nArguments\n\nio: the stream where to load the image from.\nheader=false: specify if the file contains the header.\nrow_major=false: tells if the data is stored in row major ordering or not.\nheader_type=Int64: the type of the header data.\n\nSee Also: load_sinogram, read_ct_image, write_ct_image.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.write_tomogram-Union{Tuple{T}, Tuple{AbstractString,AbstractArray{T,2}}} where T<:Real","page":"Marta","title":"Marta.CTIO.write_tomogram","text":"write_tomogram(f::AbstractString, tomog::AbstractMatrix{T}; <keyword arguments>) where {T<:Real}\n\nWrite reconstructed image to file f.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.write_tomogram-Union{Tuple{T}, Tuple{IO,AbstractArray{T,2}}} where T<:Real","page":"Marta","title":"Marta.CTIO.write_tomogram","text":"write_tomogram(io::IO, tomog::AbstractMatrix{T}; <keyword arguments>) where {T<:Real})\n\nWrite reconstructed image to stream io.\n\nArguments\n\nio: the stream where to load the image from.\nheader=false: specify if the file contains the header.\nrow_major=false: tells if the data is stored in row major ordering or not.\nheader_type=Int64: the type of the header data.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTScan.AbstractCTScanner","page":"Marta","title":"Marta.CTScan.AbstractCTScanner","text":"abstract type AbstractCTScanner{<:AbstractArray,<:AbstractGeometry} end\n\nAbstract type for tests in this suite.\n\n\n\n\n\n","category":"type"},{"location":"#Marta.CTScan.FBPScanner-Tuple{AbstractCTScanner}","page":"Marta","title":"Marta.CTScan.FBPScanner","text":"FBPScanner(gst::AbstractCTScanner; name::Optional{String}, study_id::Optional{String})\n\nConstruct a FBPScanner object from another AbstractCTScanner object. This allows to reuse data.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTScan.FBPScanner-Tuple{AbstractGeometry,AbstractTestImage}","page":"Marta","title":"Marta.CTScan.FBPScanner","text":"FBPScanner(\n    geometry::AbstractGeometry,\n    image::AbstractTestImage;\n    name::Optional{String},\n    study_id::Optional{String}\n)\n\nConstruct a FBPScanner object.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTScan.FBPScanner-Tuple{AbstractTestImage}","page":"Marta","title":"Marta.CTScan.FBPScanner","text":"FBPScanner(\n    image::AbstractTestImage;\n    name::Optional{String},\n    study_id::Optional{String}\n)\n\nConstruct a FBPScanner object.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTScan.FBPScanner-Union{Tuple{D}, Tuple{G}, Tuple{G,D}} where D<:Marta.CTScan.CTImageData.AbstractCTData where G<:AbstractGeometry","page":"Marta","title":"Marta.CTScan.FBPScanner","text":"FBPScanner(\n    geometry::AbstractGeometry,\n    data::AbstractCTData;\n    name::Optional{String},\n    study_id::Optional{String}\n)\n\nConstruct a FBPScanner object.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.AbstractAlgorithms.project_image-Union{Tuple{AbstractCTScanner}, Tuple{A}, Tuple{AbstractCTScanner,Union{Nothing, A}}} where A<:AbstractProjectionAlgorithm","page":"Marta","title":"Marta.AbstractAlgorithms.project_image","text":"project_image(gst::AbstractCTScanner, alg::Optional{A}; <keyword arguments>)\n    where {A <: AbstractProjectionAlgorithm}\n\nCompute sinogram for the test gst.\n\nKeyword arguments depend on the algorithm employed, please see the relative documentation.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.AbstractAlgorithms.reconstruct_image-Tuple{AbstractCTScanner}","page":"Marta","title":"Marta.AbstractAlgorithms.reconstruct_image","text":"reconstruct_image(gst::AbstractCTScanner; <keyword arguments>)\n\nReconstruct CT image from test gst.\n\nKeyword arguments depend on the specific algorithm used, please see the relative documentation.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.load_image-Tuple{AbstractString,AbstractCTScanner}","page":"Marta","title":"Marta.CTIO.load_image","text":"load_image(f::AbstractString, gst::AbstractCTScanner)\n\nLoad image from file f.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.load_image-Tuple{IO,AbstractCTScanner}","page":"Marta","title":"Marta.CTIO.load_image","text":"load_image(io::IO, gst::AbstractCTScanner)\n\nLoad image from stream io.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.load_sinogram-Tuple{AbstractString,AbstractCTScanner}","page":"Marta","title":"Marta.CTIO.load_sinogram","text":"load_sinogram(f::AbstractString, gst::AbstractCTScanner)\n\nLoad sinogram from file f.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.load_sinogram-Tuple{IO,AbstractCTScanner}","page":"Marta","title":"Marta.CTIO.load_sinogram","text":"load_sinogram(io::IO, gst::AbstractCTScanner)\n\nLoad sinogram from stream io.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.load_tomogram-Tuple{AbstractString,AbstractCTScanner}","page":"Marta","title":"Marta.CTIO.load_tomogram","text":"load_tomogram(f::AbstractString, gst::AbstractCTScanner)\n\nLoad tomogram from file f.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.load_tomogram-Tuple{IO,AbstractCTScanner}","page":"Marta","title":"Marta.CTIO.load_tomogram","text":"load_tomogram(io::IO, gst::AbstractCTScanner)\n\nLoad tomogram from stream io.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.write_image-Tuple{AbstractString,AbstractCTScanner}","page":"Marta","title":"Marta.CTIO.write_image","text":"write_image(f::AbstractString, gst::AbstractCTScanner)\n\nWrite input image to file f.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.write_image-Tuple{IO,AbstractCTScanner}","page":"Marta","title":"Marta.CTIO.write_image","text":"write_image(io::IO, gst::AbstractCTScanner)\n\nWrite input image to stream io.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.write_sinogram-Tuple{AbstractString,AbstractCTScanner}","page":"Marta","title":"Marta.CTIO.write_sinogram","text":"write_sinogram(f::AbstractString, gst::AbstractCTScanner)\n\nWrite sinogram to file f.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.write_sinogram-Tuple{IO,AbstractCTScanner}","page":"Marta","title":"Marta.CTIO.write_sinogram","text":"write_sinogram(io::IO, gst::AbstractCTScanner)\n\nWrite sinogram to stream io.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.write_tomogram-Tuple{AbstractString,AbstractCTScanner}","page":"Marta","title":"Marta.CTIO.write_tomogram","text":"write_tomogram(f::AbstractString, gst::AbstractCTScanner)\n\nWrite tomogram to file f.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTIO.write_tomogram-Tuple{IO,AbstractCTScanner}","page":"Marta","title":"Marta.CTIO.write_tomogram","text":"write_tomogram(io::IO, gst::AbstractCTScanner)\n\nWrite reconstructed image to stream io.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTScan.compute_gray_scale-Tuple{AbstractArray{T,2} where T,ImageParams}","page":"Marta","title":"Marta.CTScan.compute_gray_scale","text":"compute_gray_scale(image::AbstractMatrix, imp::ImageParams; mean = false)\n\nCompute the corresponding gray scale for the test image image.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTScan.compute_gray_scale-Tuple{AbstractCTScanner,AbstractGrayScale}","page":"Marta","title":"Marta.CTScan.compute_gray_scale","text":"compute_gray_scale(gst::AbstractCTScanner, imp::ImageParams; mean = false)\n\nCompute the corresponding gray scale for the reconstructed image.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTScan.create_image!-Tuple{AbstractCTScanner,ImageParams}","page":"Marta","title":"Marta.CTScan.create_image!","text":"create_image!(gst::AbstractCTScanner, par::ImageParams)\n\nCreate gray scale image for the test gst.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTScan.sample_sinogram-Tuple{AbstractCTScanner}","page":"Marta","title":"Marta.CTScan.sample_sinogram","text":"sample_sinogram(gst::AbstractCTScanner; <keyword arguments>)\n\nReturns a new test object with the resampled sinogram.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CalibrationBase.calibrate_image-Union{Tuple{AbstractCTScanner}, Tuple{T}} where T<:Real","page":"Marta","title":"Marta.CalibrationBase.calibrate_image","text":"calibrate_image(gst::AbstractCTScanner; min_pos, max_pos, interval=0..1, window=nothing)\n\nPerform calibration of input image with reference values.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CalibrationBase.calibrate_image-Union{Tuple{U}, Tuple{T}, Tuple{AbstractCTScanner,AbstractImageParams}} where U<:Real where T<:Real","page":"Marta","title":"Marta.CalibrationBase.calibrate_image","text":"calibrate_image(gst::AbstractCTScanner, imp::ImageParams; interval=nothing, window=nothing)\n\nPerform calibration of input image using image parameters as reference.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CalibrationBase.calibrate_tomogram-Union{Tuple{AbstractCTScanner}, Tuple{T}} where T<:Real","page":"Marta","title":"Marta.CalibrationBase.calibrate_tomogram","text":"calibrate_tomogram(gst::AbstractCTScanner; min_pos, max_pos, interval=0..1, window=nothing)\n\nPerform calibration of reconstructed image with reference values.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CalibrationBase.calibrate_tomogram-Union{Tuple{U}, Tuple{T}, Tuple{AbstractCTScanner,AbstractImageParams}} where U<:Real where T<:Real","page":"Marta","title":"Marta.CalibrationBase.calibrate_tomogram","text":"calibrate_tomogram(gst::AbstractCTScanner, imp::ImageParams; interval=nothing, window=nothing)\n\nPerform calibration of reconstructed image using image parameters as reference.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTPlots.plot_gray_scale-Tuple{ImageParams,Vararg{AbstractArray{T,1} where T,N} where N}","page":"Marta","title":"Marta.CTPlots.plot_gray_scale","text":"plot_gray_scale(\n    imp::ImageParams,\n    gray_scale_data...;\n    options = nothing,\n    <keyword arguments>\n)\n\nPlot gray scale data.\n\nAdditional keyword arguments are passed to the Plots library.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTPlots.plot_image-Tuple{AbstractCTScanner}","page":"Marta","title":"Marta.CTPlots.plot_image","text":"plot_image(gst::AbstractCTScanner)\n\nPlot input image created for the test gst.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTPlots.plot_image-Union{Tuple{AbstractArray{T,2}}, Tuple{T}} where T","page":"Marta","title":"Marta.CTPlots.plot_image","text":"plot_image(image::AbstractMatrix{T}; <keyword arguments>) where {T}\n\nPlot CT image image with predefined options.\n\nAdditional keyword arguments are passed to the Plots library functions.\n\nSee also: plot_sinogram, plot_tomogram\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTPlots.plot_sinogram-Tuple{AbstractCTScanner}","page":"Marta","title":"Marta.CTPlots.plot_sinogram","text":"plot_sinogram(gst::AbstractCTScanner)\n\nPlot the sinogram computed for the test gst.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTPlots.plot_sinogram-Union{Tuple{AbstractArray{T,2}}, Tuple{T}} where T","page":"Marta","title":"Marta.CTPlots.plot_sinogram","text":"plot_sinogram(sinog::AbstractMatrix{T}; <keyword arguments>) where {T}\n\nPlot sinogram sinog with predefined options.\n\nAdditional keyword arguments are passed to the underlying plotting driver.\n\nSee also: plot_image, plot_tomogram\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTPlots.plot_tomogram-Tuple{AbstractCTScanner}","page":"Marta","title":"Marta.CTPlots.plot_tomogram","text":"plot_tomogram(gst::AbstractCTScanner)\n\nPlot the reconstructed image for the test gst.\n\n\n\n\n\n","category":"method"},{"location":"#Marta.CTPlots.plot_tomogram-Union{Tuple{AbstractArray{T,2}}, Tuple{T}} where T","page":"Marta","title":"Marta.CTPlots.plot_tomogram","text":"plot_tomogram(tomog::AbstractMatrix{T}; <keyword arguments>) where {T}\n\nPlot reconstructed CT image tomog with predefined options.\n\nAlias for plot_image. See plot_image for full reference. Additional keyword arguments are passed to the Plots library.\n\nSee also: plot_sinogram, plot_image\n\n\n\n\n\n","category":"method"}]
}
