# Use
#    @warnpcfail precompile(args...)
# if you want to be warned when a precompile directive fails
macro warnpcfail(ex::Expr)
    modl = __module__
    file = __source__.file === nothing ? "?" : String(__source__.file)
    line = __source__.line
    quote
        $(esc(ex)) || @warn """precompile directive
     $($(Expr(:quote, ex)))
 failed. Please report an issue in $($modl) (after checking for duplicates) or remove this directive.""" _file=$file _line=$line
    end
end


function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id), Tuple{Nothing, Nothing}},Type{FBPScanner{G, D} where {G, D}},ParallelBeamGeometry{Float32, DefaultTomograph},CTImageMat{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id), Tuple{Nothing, Nothing}},Type{FBPScanner{G, D} where {G, D}},ParallelBeamGeometry{Float64, DefaultTomograph},CTImageMat{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,), Tuple{Int64}},Type{FanBeamGeometry},CircleImage{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,), Tuple{Int64}},Type{FanBeamGeometry},CircleImage{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,), Tuple{Int64}},Type{ParallelBeamGeometry},CircleImage{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,), Tuple{Int64}},Type{ParallelBeamGeometry},CircleImage{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width, :height, :nϕ), Tuple{Int64, Int64, Int64}},Type{ParallelBeamGeometry},Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,), Tuple{Int64}},Type{CircleImage},Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,), Tuple{Int64}},Type{CircleImage},Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,), Tuple{Int64}},Type{GrayScaleLine},Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,), Tuple{Int64}},Type{GrayScaleLine},Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,), Tuple{Int64}},Type{GrayScalePyramid},Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,), Tuple{Int64}},Type{GrayScalePyramid},Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,), Tuple{Int64}},Type{WhiteRect},Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,), Tuple{Int64}},Type{WhiteRect},Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(circle_image)),NamedTuple{(:radius, :calibration_value, :background, :rows), Tuple{Int64, Float32, Float32, Int64}},typeof(circle_image),Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(circle_image)),NamedTuple{(:radius, :calibration_value, :background, :rows), Tuple{Int64, Float64, Float64, Int64}},typeof(circle_image),Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(gray_scale_image)),NamedTuple{(:swidth, :sheight, :gray_scale, :background), Tuple{Int64, Int64, IntervalSets.ClosedInterval{Float32}, Float32}},typeof(gray_scale_image),Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(gray_scale_image)),NamedTuple{(:swidth, :sheight, :gray_scale, :background), Tuple{Int64, Int64, IntervalSets.ClosedInterval{Float64}, Float64}},typeof(gray_scale_image),Type{Float64}})
    Base.precompile(Tuple{Marta.CTImages.var"#168#threadsfor_fun#13"{Float32, StepRangeLen{Float32, Float64, Float64}, StepRangeLen{Float32, Float64, Float64}, Matrix{Float32}, Vector{Tuple{Int64, Int64}}, Marta.Interpolation.var"#1#2"{Matrix{Float32}, Marta.Interpolation.BilinearInterpolation}, Float32, Float32, Base.OneTo{Int64}}})
    Base.precompile(Tuple{Marta.CTImages.var"#168#threadsfor_fun#13"{Float32, StepRangeLen{Float32, Float64, Float64}, StepRangeLen{Float32, Float64, Float64}, Matrix{Float32}, Vector{Tuple{Int64, Int64}}, Marta.Interpolation.var"#1#2"{Matrix{Float32}, Marta.Interpolation.NoInterpolation}, Float32, Float32, Base.OneTo{Int64}}})
    Base.precompile(Tuple{Marta.CTImages.var"#168#threadsfor_fun#13"{Float64, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}, Matrix{Float64}, Vector{Tuple{Int64, Int64}}, Marta.Interpolation.var"#1#2"{Matrix{Float64}, Marta.Interpolation.BilinearInterpolation}, Float64, Float64, Base.OneTo{Int64}}})
    Base.precompile(Tuple{Marta.CTImages.var"#168#threadsfor_fun#13"{Float64, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}, Matrix{Float64}, Vector{Tuple{Int64, Int64}}, Marta.Interpolation.var"#1#2"{Matrix{Float64}, Marta.Interpolation.NoInterpolation}, Float64, Float64, Base.OneTo{Int64}}})
    Base.precompile(Tuple{Marta.FanBeam.var"#198#threadsfor_fun#3"{Float32, Marta.FanBeam.var"#compute_value#2"{Float32, Float32, Marta.Interpolation.var"#1#2"{CTSinogMat{Float32}, Marta.Interpolation.BilinearInterpolation}, Float32, Float32, Float32, Int64, Int64}, CTSinogMat{Float32}, StepRangeLen{Float32, Float64, Float64}, StepRangeLen{Float32, Float64, Float64}, Float32, Int64, UnitRange{Int64}}})
    Base.precompile(Tuple{Marta.FanBeam.var"#198#threadsfor_fun#3"{Float64, Marta.FanBeam.var"#compute_value#2"{Float64, Float64, Marta.Interpolation.var"#1#2"{CTSinogMat{Float64}, Marta.Interpolation.BilinearInterpolation}, Float64, Float64, Float64, Int64, Int64}, CTSinogMat{Float64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}, Float64, Int64, UnitRange{Int64}}})
    Base.precompile(Tuple{Marta.FanBeam.var"#218#threadsfor_fun#10"{Float32, Marta.FanBeam.var"#compute_value#9"{Float32, Float32, Marta.Interpolation.var"#1#2"{CTSinogMat{Float32}, Marta.Interpolation.BilinearInterpolation}, Float32, Float32, Float32, Int64, Int64}, CTSinogMat{Float32}, StepRangeLen{Float32, Float64, Float64}, Float32, StepRangeLen{Float32, Float64, Float64}, Float32, Int64, UnitRange{Int64}}})
    Base.precompile(Tuple{Marta.FanBeam.var"#218#threadsfor_fun#10"{Float64, Marta.FanBeam.var"#compute_value#9"{Float64, Float64, Marta.Interpolation.var"#1#2"{CTSinogMat{Float64}, Marta.Interpolation.BilinearInterpolation}, Float64, Float64, Float64, Int64, Int64}, CTSinogMat{Float64}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}, Float64, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}, Float64, Int64, UnitRange{Int64}}})
    Base.precompile(Tuple{Type{FBPScanner{G, D} where {G, D}},CircleImage{Float32}})
    Base.precompile(Tuple{Type{FBPScanner{G, D} where {G, D}},CircleImage{Float64}})
    Base.precompile(Tuple{Type{FBPScanner{G, D} where {G, D}},GrayScaleLine{Float32}})
    Base.precompile(Tuple{Type{FBPScanner{G, D} where {G, D}},GrayScaleLine{Float64}})
    Base.precompile(Tuple{Type{FBPScanner{G, D} where {G, D}},GrayScalePyramid{Float32}})
    Base.precompile(Tuple{Type{FBPScanner{G, D} where {G, D}},GrayScalePyramid{Float64}})
    Base.precompile(Tuple{Type{FBPScanner{G, D} where {G, D}},WhiteRect{Float32}})
    Base.precompile(Tuple{Type{FBPScanner{G, D} where {G, D}},WhiteRect{Float64}})
    Base.precompile(Tuple{typeof(Marta.Utils._wrap_angle),Float32,Int64})
    Base.precompile(Tuple{typeof(Marta.Utils._wrap_angle),Float64,Int64})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBPInfo{ParallelBeamGeometry{Float32, DefaultTomograph}, FBP{Marta.RadonAlgorithm.Filters.RamLak}}, Marta.CTScan.CTImageData.CTData{Float32, Matrix{Float32}}},CircleImage{Float32}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBPInfo{ParallelBeamGeometry{Float32, DefaultTomograph}, FBP{Marta.RadonAlgorithm.Filters.RamLak}}, Marta.CTScan.CTImageData.CTData{Float32, Matrix{Float32}}},CircleParams{Float32}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBPInfo{ParallelBeamGeometry{Float32, DefaultTomograph}, FBP{Marta.RadonAlgorithm.Filters.RamLak}}, Marta.CTScan.CTImageData.CTData{Float32, Matrix{Float32}}},ImageParams{Float32}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBPInfo{ParallelBeamGeometry{Float64, DefaultTomograph}, FBP{Marta.RadonAlgorithm.Filters.RamLak}}, Marta.CTScan.CTImageData.CTData{Float64, Matrix{Float64}}},CircleImage{Float64}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBPInfo{ParallelBeamGeometry{Float64, DefaultTomograph}, FBP{Marta.RadonAlgorithm.Filters.RamLak}}, Marta.CTScan.CTImageData.CTData{Float64, Matrix{Float64}}},CircleParams{Float64}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBPInfo{ParallelBeamGeometry{Float64, DefaultTomograph}, FBP{Marta.RadonAlgorithm.Filters.RamLak}}, Marta.CTScan.CTImageData.CTData{Float64, Matrix{Float64}}},GrayScaleLine{Float64}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBPInfo{ParallelBeamGeometry{Float64, DefaultTomograph}, FBP{Marta.RadonAlgorithm.Filters.RamLak}}, Marta.CTScan.CTImageData.CTData{Float64, Matrix{Float64}}},ImageParams{Float64}})
    Base.precompile(Tuple{typeof(combine_images),ImageParams{Float32},Matrix{Float32},Matrix{Float32}})
    Base.precompile(Tuple{typeof(combine_images),ImageParams{Float64},Matrix{Float64},Matrix{Float64}})
end
