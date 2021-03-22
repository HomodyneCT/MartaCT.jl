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
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id),Tuple{Symbol,String}},Type{CTScanner},FBPInfo{ParallelBeamGeometry{Float64,DefaultTomograph},FBP{Marta.RadonAlgorithm.Filters.RamLak}},Marta.CTScan.CTImageData.CTData{Float64,Array{Float64,2}}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,),Tuple{Int64}},Type{FanBeamGeometry},CircleImage{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,),Tuple{Int64}},Type{FanBeamGeometry},CircleImage{Float64}})
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
    Base.precompile(Tuple{Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},CircleImage{Float32}})
    Base.precompile(Tuple{Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},CircleImage{Float64}})
    Base.precompile(Tuple{Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},FanBeamGeometry{Float64,DefaultTomograph},CircleImage{Float64}})
    Base.precompile(Tuple{Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},GrayScalePyramid{Float32}})
    Base.precompile(Tuple{Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},GrayScalePyramid{Float64}})
    Base.precompile(Tuple{Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},WhiteRect{Float32}})
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
