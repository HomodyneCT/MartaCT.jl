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


const __bodyfunction__ = Dict{Method,Any}()

# Find keyword "body functions" (the function that contains the body
# as written by the developer, called after all missing keyword-arguments
# have been assigned values), in a manner that doesn't depend on
# gensymmed names.
# `mnokw` is the method that gets called when you invoke it without
# supplying any keywords.
function __lookup_kwbody__(mnokw::Method)
    function getsym(arg)
        isa(arg, Symbol) && return arg
        @assert isa(arg, GlobalRef)
        return arg.name
    end

    f = get(__bodyfunction__, mnokw, nothing)
    if f === nothing
        fmod = mnokw.module
        # The lowered code for `mnokw` should look like
        #   %1 = mkw(kwvalues..., #self#, args...)
        #        return %1
        # where `mkw` is the name of the "active" keyword body-function.
        ast = Base.uncompressed_ast(mnokw)
        if isa(ast, Core.CodeInfo) && length(ast.code) >= 2
            callexpr = ast.code[end-1]
            if isa(callexpr, Expr) && callexpr.head == :call
                fsym = callexpr.args[1]
                if isa(fsym, Symbol)
                    f = getfield(fmod, fsym)
                elseif isa(fsym, GlobalRef)
                    if fsym.mod === Core && fsym.name === :_apply
                        f = getfield(mnokw.module, getsym(callexpr.args[2]))
                    elseif fsym.mod === Core && fsym.name === :_apply_iterate
                        f = getfield(mnokw.module, getsym(callexpr.args[3]))
                    else
                        f = getfield(fsym.mod, fsym.name)
                    end
                else
                    f = missing
                end
            else
                f = missing
            end
        else
            f = missing
        end
        __bodyfunction__[mnokw] = f
    end
    return f
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
    Base.precompile(Tuple{Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},CircleImage{Float64}})
    Base.precompile(Tuple{Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},GrayScaleLine{Float64}})
    Base.precompile(Tuple{Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},GrayScalePyramid{Float32}})
    Base.precompile(Tuple{Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},GrayScalePyramid{Float64}})
    Base.precompile(Tuple{Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},WhiteRect{Float32}})
    Base.precompile(Tuple{Type{CTScanner{FBPInfo{G,A} where A<:AbstractFBP,D} where D where G},WhiteRect{Float64}})
    Base.precompile(Tuple{typeof(Marta.Utils._wrap_angle),Float32,Int64})
    Base.precompile(Tuple{typeof(Marta.Utils._wrap_angle),Float64,Int64})
    Base.precompile(Tuple{typeof(combine_images),ImageParams{Float32},Array{Float32,2},Array{Float32,2}})
    Base.precompile(Tuple{typeof(combine_images),ImageParams{Float64},Array{Float64,2},Array{Float64,2}})
    Base.precompile(Tuple{typeof(project_image),Some{CTImage{Float32,Array{Float32,2}}},RadonInfo{FanBeamGeometry{Float32,DefaultTomograph},Radon}})
    Base.precompile(Tuple{typeof(project_image),Some{CTImage{Float32,Array{Float32,2}}},RadonInfo{ParallelBeamGeometry{Float32,DefaultTomograph},Radon}})
    Base.precompile(Tuple{typeof(project_image),Some{CTImage{Float64,Array{Float64,2}}},RadonInfo{FanBeamGeometry{Float64,DefaultTomograph},Radon}})
    Base.precompile(Tuple{typeof(project_image),Some{CTImage{Float64,Array{Float64,2}}},RadonInfo{ParallelBeamGeometry{Float64,DefaultTomograph},Radon}})
    Base.precompile(Tuple{typeof(similar),Marta.CTScan.CTImageData.CTData{Float32,Array{Float32,2}},CTSinogram{Float32,Array{Float32,2}}})
    Base.precompile(Tuple{typeof(similar),Marta.CTScan.CTImageData.CTData{Float64,Array{Float64,2}},CTSinogram{Float64,Array{Float64,2}}})
    let fbody = try __lookup_kwbody__(which(polar2cart, (Array{Float32,2},StepRangeLen{Float32,Float64,Float64},StepRangeLen{Float32,Float64,Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Int64,Bool,Float32,Bool,Marta.Interpolation.NoInterpolation,typeof(polar2cart),Array{Float32,2},StepRangeLen{Float32,Float64,Float64},StepRangeLen{Float32,Float64,Float64},))
        end
    end
    let fbody = try __lookup_kwbody__(which(polar2cart, (Array{Float32,2},StepRangeLen{Float32,Float64,Float64},StepRangeLen{Float32,Float64,Float64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Int64,Bool,Float32,Bool,Nothing,typeof(polar2cart),Array{Float32,2},StepRangeLen{Float32,Float64,Float64},StepRangeLen{Float32,Float64,Float64},))
        end
    end
    let fbody = try __lookup_kwbody__(which(polar2cart, (Array{Float64,2},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Int64,Bool,Float64,Bool,Marta.Interpolation.NoInterpolation,typeof(polar2cart),Array{Float64,2},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},))
        end
    end
    let fbody = try __lookup_kwbody__(which(polar2cart, (Array{Float64,2},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Int64,Bool,Float64,Bool,Nothing,typeof(polar2cart),Array{Float64,2},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},))
        end
    end
end
