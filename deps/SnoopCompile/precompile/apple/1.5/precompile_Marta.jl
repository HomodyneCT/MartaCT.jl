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
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id),Tuple{Nothing,Nothing}},Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},FanBeamGeometry{Float32,DefaultTomograph},CTImage{Float32,Array{Float32,2}}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id),Tuple{Nothing,Nothing}},Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},FanBeamGeometry{Float64,DefaultTomograph},CTImage{Float64,Array{Float64,2}}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id),Tuple{Nothing,Nothing}},Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},ParallelBeamGeometry{Float32,DefaultTomograph},CTImage{Float32,Array{Float32,2}}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id),Tuple{Nothing,Nothing}},Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},ParallelBeamGeometry{Float64,DefaultTomograph},CTImage{Float64,Array{Float64,2}}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id),Tuple{Symbol,String}},Type{CTScanner},FBPInfo{ParallelBeamGeometry{Float32,DefaultTomograph},FBP{Marta.RadonAlgorithm.Filters.RamLak}},Marta.CTScan.CTImageData.CTData{Float32,Array{Float32,2}}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,),Tuple{Int64}},Type{FanBeamGeometry},CircleImage{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,),Tuple{Int64}},Type{FanBeamGeometry},CircleImage{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,),Tuple{Int64}},Type{FanBeamGeometry},WhiteRect{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,),Tuple{Int64}},Type{ParallelBeamGeometry},CircleImage{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,),Tuple{Int64}},Type{ParallelBeamGeometry},CircleImage{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{CircleImage},Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{CircleImage},Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{GrayScaleLine},Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{GrayScaleLine},Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{GrayScalePyramid},Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{GrayScalePyramid},Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{WhiteRect},Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{WhiteRect},Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(circle_image)),NamedTuple{(:radius, :calibration_value, :background, :rows),Tuple{Int64,Float32,Float32,Int64}},typeof(circle_image),Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(circle_image)),NamedTuple{(:radius, :calibration_value, :background, :rows),Tuple{Int64,Float64,Float64,Int64}},typeof(circle_image),Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(gray_scale_image)),NamedTuple{(:swidth, :sheight, :gray_scale, :background),Tuple{Int64,Int64,IntervalSets.Interval{:closed,:closed,Float32},Float32}},typeof(gray_scale_image),Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(gray_scale_image)),NamedTuple{(:swidth, :sheight, :gray_scale, :background),Tuple{Int64,Int64,IntervalSets.Interval{:closed,:closed,Float64},Float64}},typeof(gray_scale_image),Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(polar2cart)),NamedTuple{(:rows, :cols, :background, :ν),Tuple{Int64,Int64,Float32,Int64}},typeof(polar2cart),Array{Float32,2}})
    Base.precompile(Tuple{Core.kwftype(typeof(polar2cart)),NamedTuple{(:rows, :cols, :background, :ν),Tuple{Int64,Int64,Float64,Int64}},typeof(polar2cart),Array{Float64,2}})
    Base.precompile(Tuple{Core.kwftype(typeof(polar2cart)),NamedTuple{(:rows, :cols, :background, :ν, :interpolation),Tuple{Int64,Int64,Float32,Float64,Marta.Interpolation.NoInterpolation}},typeof(polar2cart),Array{Float32,2}})
    Base.precompile(Tuple{Core.kwftype(typeof(polar2cart)),NamedTuple{(:rows, :cols, :background, :ν, :interpolation),Tuple{Int64,Int64,Float64,Float64,Marta.Interpolation.NoInterpolation}},typeof(polar2cart),Array{Float64,2}})
    Base.precompile(Tuple{Marta.CTImages.var"#147#threadsfor_fun#14"{Float32,StepRangeLen{Float32,Float64,Float64},StepRangeLen{Float32,Float64,Float64},Int64,Int64,Float32,Float32,Marta.Interpolation.var"#1#2"{Array{Float32,2},Marta.Interpolation.BilinearInterpolation},Array{Tuple{Int64,Int64},1},Array{Float32,2},Marta.CTImages.var"#_compute_radius#13",Base.OneTo{Int64}}})
    Base.precompile(Tuple{Marta.CTImages.var"#147#threadsfor_fun#14"{Float32,StepRangeLen{Float32,Float64,Float64},StepRangeLen{Float32,Float64,Float64},Int64,Int64,Float32,Float32,Marta.Interpolation.var"#1#2"{Array{Float32,2},Marta.Interpolation.NoInterpolation},Array{Tuple{Int64,Int64},1},Array{Float32,2},Marta.CTImages.var"#_compute_radius#13",Base.OneTo{Int64}}})
    Base.precompile(Tuple{Marta.CTImages.var"#147#threadsfor_fun#14"{Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},Int64,Int64,Float64,Float64,Marta.Interpolation.var"#1#2"{Array{Float64,2},Marta.Interpolation.BilinearInterpolation},Array{Tuple{Int64,Int64},1},Array{Float64,2},Marta.CTImages.var"#_compute_radius#13",Base.OneTo{Int64}}})
    Base.precompile(Tuple{Marta.CTImages.var"#147#threadsfor_fun#14"{Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},Int64,Int64,Float64,Float64,Marta.Interpolation.var"#1#2"{Array{Float64,2},Marta.Interpolation.NoInterpolation},Array{Tuple{Int64,Int64},1},Array{Float64,2},Marta.CTImages.var"#_compute_radius#13",Base.OneTo{Int64}}})
    Base.precompile(Tuple{Marta.FanBeam.var"#170#threadsfor_fun#3"{Float32,Int64,Float32,Float32,Float32,Float32,Float32,CTSinogram{Float32,Array{Float32,2}},Marta.FanBeam.var"#compute_value#2"{Float32,Int64,Int64,Float32,Float32,Float32,Marta.Interpolation.var"#1#2"{CTSinogram{Float32,Array{Float32,2}},Marta.Interpolation.BilinearInterpolation},Float32},UnitRange{Int64}}})
    Base.precompile(Tuple{Marta.FanBeam.var"#170#threadsfor_fun#3"{Float64,Int64,Float64,Float64,Float64,Float64,Float64,CTSinogram{Float64,Array{Float64,2}},Marta.FanBeam.var"#compute_value#2"{Float64,Int64,Int64,Float64,Float64,Float64,Marta.Interpolation.var"#1#2"{CTSinogram{Float64,Array{Float64,2}},Marta.Interpolation.BilinearInterpolation},Float64},UnitRange{Int64}}})
    Base.precompile(Tuple{Marta.FanBeam.var"#189#threadsfor_fun#10"{Float32,Int64,Float32,Float32,Float32,Float32,Float32,CTSinogram{Float32,Array{Float32,2}},Marta.FanBeam.var"#compute_value#9"{Float32,Int64,Int64,Float32,Float32,Float32,Marta.Interpolation.var"#1#2"{CTSinogram{Float32,Array{Float32,2}},Marta.Interpolation.BilinearInterpolation},Float32},UnitRange{Int64}}})
    Base.precompile(Tuple{Marta.FanBeam.var"#189#threadsfor_fun#10"{Float64,Int64,Float64,Float64,Float64,Float64,Float64,CTSinogram{Float64,Array{Float64,2}},Marta.FanBeam.var"#compute_value#9"{Float64,Int64,Int64,Float64,Float64,Float64,Marta.Interpolation.var"#1#2"{CTSinogram{Float64,Array{Float64,2}},Marta.Interpolation.BilinearInterpolation},Float64},UnitRange{Int64}}})
    Base.precompile(Tuple{Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},CircleImage{Float32}})
    Base.precompile(Tuple{Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},CircleImage{Float64}})
    Base.precompile(Tuple{Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},GrayScalePyramid{Float32}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBPInfo{ParallelBeamGeometry{Float32,DefaultTomograph},FBP{Marta.RadonAlgorithm.Filters.RamLak}},Marta.CTScan.CTImageData.CTData{Float32,Array{Float32,2}}},CircleImage{Float32}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBPInfo{ParallelBeamGeometry{Float32,DefaultTomograph},FBP{Marta.RadonAlgorithm.Filters.RamLak}},Marta.CTScan.CTImageData.CTData{Float32,Array{Float32,2}}},CircleParams{Float32}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBPInfo{ParallelBeamGeometry{Float32,DefaultTomograph},FBP{Marta.RadonAlgorithm.Filters.RamLak}},Marta.CTScan.CTImageData.CTData{Float32,Array{Float32,2}}},ImageParams{Float32}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBPInfo{ParallelBeamGeometry{Float64,DefaultTomograph},FBP{Marta.RadonAlgorithm.Filters.RamLak}},Marta.CTScan.CTImageData.CTData{Float64,Array{Float64,2}}},CircleImage{Float64}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBPInfo{ParallelBeamGeometry{Float64,DefaultTomograph},FBP{Marta.RadonAlgorithm.Filters.RamLak}},Marta.CTScan.CTImageData.CTData{Float64,Array{Float64,2}}},CircleParams{Float64}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBPInfo{ParallelBeamGeometry{Float64,DefaultTomograph},FBP{Marta.RadonAlgorithm.Filters.RamLak}},Marta.CTScan.CTImageData.CTData{Float64,Array{Float64,2}}},ImageParams{Float64}})
    Base.precompile(Tuple{typeof(project_image),Some{CTImage{Float32,Array{Float32,2}}},RadonInfo{FanBeamGeometry{Float32,DefaultTomograph},Radon}})
    Base.precompile(Tuple{typeof(project_image),Some{CTImage{Float32,Array{Float32,2}}},RadonInfo{ParallelBeamGeometry{Float32,DefaultTomograph},Radon}})
    Base.precompile(Tuple{typeof(project_image),Some{CTImage{Float64,Array{Float64,2}}},RadonInfo{FanBeamGeometry{Float64,DefaultTomograph},Radon}})
    Base.precompile(Tuple{typeof(project_image),Some{CTImage{Float64,Array{Float64,2}}},RadonInfo{ParallelBeamGeometry{Float64,DefaultTomograph},Radon}})
end
