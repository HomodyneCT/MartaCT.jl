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
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id),Tuple{Nothing,Nothing}},Type{CTScanner{FBP,M} where M},FanBeamGeometry{Float32,DefaultTomograph},CTImage{Array{Float32,2}}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id),Tuple{Nothing,Nothing}},Type{CTScanner{FBP,M} where M},FanBeamGeometry{Float64,DefaultTomograph},CTImage{Array{Float64,2}}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id),Tuple{Nothing,Nothing}},Type{CTScanner{FBP,M} where M},ParallelBeamGeometry{Float32,DefaultTomograph},CTImage{Array{Float32,2}}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id),Tuple{Nothing,Nothing}},Type{CTScanner{FBP,M} where M},ParallelBeamGeometry{Float64,DefaultTomograph},CTImage{Array{Float64,2}}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id),Tuple{Symbol,String}},Type{CTScanner},FBP{ParallelBeamGeometry{Float32,DefaultTomograph},Marta.Filters.RamLak},Marta.CTData.GrayScaleData{Array{Float32,2}}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,),Tuple{Int64}},Type{FanBeamGeometry},GrayScaleLine{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,),Tuple{Int64}},Type{FanBeamGeometry},GrayScaleLine{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,),Tuple{Int64}},Type{ParallelBeamGeometry},GrayScaleLine{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{CircleImage},Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{CircleImage},Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{GrayScaleLine},Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{GrayScaleLine},Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{GrayScalePyramid},Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{GrayScalePyramid},Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{WhiteRect},Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{WhiteRect},Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(calibrate_image)),NamedTuple{(:interval, :min_pos, :max_pos, :window),Tuple{IntervalSets.Interval{:closed,:closed,Float32},Tuple{Int64,Int64},Tuple{Int64,Int64},Nothing}},typeof(calibrate_image),Array{Float32,2}})
    Base.precompile(Tuple{Core.kwftype(typeof(calibrate_image)),NamedTuple{(:interval, :min_pos, :max_pos, :window),Tuple{IntervalSets.Interval{:closed,:closed,Float64},Tuple{Int64,Int64},Tuple{Int64,Int64},Nothing}},typeof(calibrate_image),Array{Float64,2}})
    Base.precompile(Tuple{Core.kwftype(typeof(circle_image)),NamedTuple{(:radius, :calibration_value, :background, :rows),Tuple{Int64,Float32,Float32,Int64}},typeof(circle_image),Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(circle_image)),NamedTuple{(:radius, :calibration_value, :background, :rows),Tuple{Int64,Float64,Float64,Int64}},typeof(circle_image),Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(gray_scale_image)),NamedTuple{(:swidth, :sheight, :gray_scale, :background),Tuple{Int64,Int64,IntervalSets.Interval{:closed,:closed,Float32},Float32}},typeof(gray_scale_image),Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(gray_scale_image)),NamedTuple{(:swidth, :sheight, :gray_scale, :background),Tuple{Int64,Int64,IntervalSets.Interval{:closed,:closed,Float64},Float64}},typeof(gray_scale_image),Type{Float64}})
    Base.precompile(Tuple{Marta.CTImages.var"#136#threadsfor_fun#41"{Float32,Int64,Int64,Float32,Float32,Marta.Interpolation.var"#1#2"{Array{Float32,2},Marta.Interpolation.BilinearInterpolation},Marta.CTImages.var"#39#40",StepRangeLen{Float32,Float64,Float64},StepRangeLen{Float32,Float64,Float64},Array{Tuple{Int64,Int64},1},Array{Float32,2},Base.OneTo{Int64}}})
    Base.precompile(Tuple{Marta.CTImages.var"#136#threadsfor_fun#41"{Float32,Int64,Int64,Float32,Float32,Marta.Interpolation.var"#1#2"{Array{Float32,2},Marta.Interpolation.NoInterpolation},Marta.CTImages.var"#39#40",StepRangeLen{Float32,Float64,Float64},StepRangeLen{Float32,Float64,Float64},Array{Tuple{Int64,Int64},1},Array{Float32,2},Base.OneTo{Int64}}})
    Base.precompile(Tuple{Marta.CTImages.var"#136#threadsfor_fun#41"{Float64,Int64,Int64,Float64,Float64,Marta.Interpolation.var"#1#2"{Array{Float64,2},Marta.Interpolation.BilinearInterpolation},Marta.CTImages.var"#39#40",StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},Array{Tuple{Int64,Int64},1},Array{Float64,2},Base.OneTo{Int64}}})
    Base.precompile(Tuple{Marta.CTImages.var"#136#threadsfor_fun#41"{Float64,Int64,Int64,Float64,Float64,Marta.Interpolation.var"#1#2"{Array{Float64,2},Marta.Interpolation.NoInterpolation},Marta.CTImages.var"#39#40",StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},Array{Tuple{Int64,Int64},1},Array{Float64,2},Base.OneTo{Int64}}})
    Base.precompile(Tuple{Marta.FanBeam.var"#161#threadsfor_fun#3"{Float32,Int64,Array{Float32,2},Float32,Float32,Float32,Float32,Float32,Marta.FanBeam.var"#compute_value#2"{Float32,Int64,Int64,Marta.Interpolation.var"#1#2"{Array{Float32,2},Marta.Interpolation.BilinearInterpolation},Float32,Float32,Float32},UnitRange{Int64}}})
    Base.precompile(Tuple{Marta.FanBeam.var"#161#threadsfor_fun#3"{Float64,Int64,Array{Float64,2},Float64,Float64,Float64,Float64,Float64,Marta.FanBeam.var"#compute_value#2"{Float64,Int64,Int64,Marta.Interpolation.var"#1#2"{Array{Float64,2},Marta.Interpolation.BilinearInterpolation},Float64,Float64,Float64},UnitRange{Int64}}})
    Base.precompile(Tuple{Marta.FanBeam.var"#180#threadsfor_fun#13"{Float32,Int64,Array{Float32,2},Float32,Float32,Float32,Float32,Float32,Marta.FanBeam.var"#compute_value#12"{Float32,Int64,Int64,Marta.Interpolation.var"#1#2"{Array{Float32,2},Marta.Interpolation.BilinearInterpolation},Float32,Float32,Float32},UnitRange{Int64}}})
    Base.precompile(Tuple{Marta.FanBeam.var"#180#threadsfor_fun#13"{Float64,Int64,Array{Float64,2},Float64,Float64,Float64,Float64,Float64,Marta.FanBeam.var"#compute_value#12"{Float64,Int64,Int64,Marta.Interpolation.var"#1#2"{Array{Float64,2},Marta.Interpolation.BilinearInterpolation},Float64,Float64,Float64},UnitRange{Int64}}})
    Base.precompile(Tuple{Type{CTScanner{FBP,M} where M},CircleImage{Float64}})
    Base.precompile(Tuple{Type{CTScanner{FBP,M} where M},GrayScaleLine{Float32}})
    Base.precompile(Tuple{Type{CTScanner{FBP,M} where M},GrayScaleLine{Float64}})
    Base.precompile(Tuple{Type{CTScanner{FBP,M} where M},GrayScalePyramid{Float32}})
    Base.precompile(Tuple{Type{CTScanner{FBP,M} where M},GrayScalePyramid{Float64}})
    Base.precompile(Tuple{Type{CTScanner{FBP,M} where M},WhiteRect{Float32}})
    Base.precompile(Tuple{typeof(Marta.FanBeam.para2fan),CTSinogram{Array{Float32,2}},FanBeamGeometry{Float32,DefaultTomograph}})
    Base.precompile(Tuple{typeof(Marta.FanBeam.para2fan),CTSinogram{Array{Float64,2}},FanBeamGeometry{Float64,DefaultTomograph}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBP{ParallelBeamGeometry{Float32,DefaultTomograph},Marta.Filters.RamLak},Array{Float32,2}},CircleParams{Float32}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBP{ParallelBeamGeometry{Float32,DefaultTomograph},Marta.Filters.RamLak},Array{Float32,2}},GrayScaleLine{Float32}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBP{ParallelBeamGeometry{Float32,DefaultTomograph},Marta.Filters.RamLak},Array{Float32,2}},ImageParams{Float32}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBP{ParallelBeamGeometry{Float64,DefaultTomograph},Marta.Filters.RamLak},Array{Float64,2}},CircleParams{Float64}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBP{ParallelBeamGeometry{Float64,DefaultTomograph},Marta.Filters.RamLak},Array{Float64,2}},ImageParams{Float64}})
    Base.precompile(Tuple{typeof(project_image),Some{CTImage{Array{Float32,2}}},Radon{FanBeamGeometry{Float32,DefaultTomograph}}})
    Base.precompile(Tuple{typeof(project_image),Some{CTImage{Array{Float32,2}}},Radon{ParallelBeamGeometry{Float32,DefaultTomograph}}})
    Base.precompile(Tuple{typeof(project_image),Some{CTImage{Array{Float64,2}}},Radon{FanBeamGeometry{Float64,DefaultTomograph}}})
    Base.precompile(Tuple{typeof(project_image),Some{CTImage{Array{Float64,2}}},Radon{ParallelBeamGeometry{Float64,DefaultTomograph}}})
    isdefined(Marta.AbstractAlgorithms, Symbol("#2#3")) && Base.precompile(Tuple{getfield(Marta.AbstractAlgorithms, Symbol("#2#3")),Array{Float32,2}})
    isdefined(Marta.AbstractAlgorithms, Symbol("#2#3")) && Base.precompile(Tuple{getfield(Marta.AbstractAlgorithms, Symbol("#2#3")),Array{Float64,2}})
end
