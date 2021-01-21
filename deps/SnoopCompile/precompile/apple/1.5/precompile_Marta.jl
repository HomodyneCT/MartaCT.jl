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
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id),Tuple{Nothing,Nothing}},Type{CTScanner{FBP,M} where M},FanBeamGeometry{Float32,DefaultTomograph},CTImage{Array{Float32,2}}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id),Tuple{Nothing,Nothing}},Type{CTScanner{FBP,M} where M},FanBeamGeometry{Float64,DefaultTomograph},CTImage{Array{Float64,2}}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id),Tuple{Nothing,Nothing}},Type{CTScanner{FBP,M} where M},ParallelBeamGeometry{Float32,DefaultTomograph},CTImage{Array{Float32,2}}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id),Tuple{Nothing,Nothing}},Type{CTScanner{FBP,M} where M},ParallelBeamGeometry{Float64,DefaultTomograph},CTImage{Array{Float64,2}}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:name, :study_id),Tuple{Symbol,String}},Type{CTScanner},FBP{ParallelBeamGeometry{Float32,DefaultTomograph},Marta.Filters.RamLak},Marta.CTData.GrayScaleData{Array{Float32,2}}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,),Tuple{Int64}},Type{FanBeamGeometry},GrayScaleLine{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,),Tuple{Int64}},Type{FanBeamGeometry},GrayScaleLine{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:nϕ,),Tuple{Int64}},Type{ParallelBeamGeometry},GrayScaleLine{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{GrayScaleLine},Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{GrayScaleLine},Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{GrayScalePyramid},Type{Float32}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{GrayScalePyramid},Type{Float64}})
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:width,),Tuple{Int64}},Type{WhiteRect},Type{Float32}})
    Base.precompile(Tuple{Marta.FanBeam.var"#103#threadsfor_fun#3"{Float32,Int64,Array{Float32,2},Float32,Float32,Float32,Float32,Float32,Marta.FanBeam.var"#compute_value#2"{Float32,Int64,Int64,Marta.Interpolation.BilinearInterpolation{Array{Float32,2}},Float32,Float32,Float32},UnitRange{Int64}}})
    Base.precompile(Tuple{Marta.FanBeam.var"#103#threadsfor_fun#3"{Float64,Int64,Array{Float64,2},Float64,Float64,Float64,Float64,Float64,Marta.FanBeam.var"#compute_value#2"{Float64,Int64,Int64,Marta.Interpolation.BilinearInterpolation{Array{Float64,2}},Float64,Float64,Float64},UnitRange{Int64}}})
    Base.precompile(Tuple{Marta.FanBeam.var"#122#threadsfor_fun#13"{Float32,Int64,Array{Float32,2},Float32,Float32,Float32,Float32,Float32,Marta.FanBeam.var"#compute_value#12"{Float32,Int64,Int64,Marta.Interpolation.BilinearInterpolation{Array{Float32,2}},Float32,Float32,Float32},UnitRange{Int64}}})
    Base.precompile(Tuple{Marta.FanBeam.var"#122#threadsfor_fun#13"{Float64,Int64,Array{Float64,2},Float64,Float64,Float64,Float64,Float64,Marta.FanBeam.var"#compute_value#12"{Float64,Int64,Int64,Marta.Interpolation.BilinearInterpolation{Array{Float64,2}},Float64,Float64,Float64},UnitRange{Int64}}})
    Base.precompile(Tuple{Marta.RadonAlgorithm.var"#149#threadsfor_fun#5"{Array{Float64,1},Array{Tuple{Float64,Float64},1},Array{Tuple{Int64,Int64},1},Array{Tuple{Float64,Float64},1},LinearIndices{2,Tuple{Base.OneTo{Int64},Base.OneTo{Int64}}}}})
    Base.precompile(Tuple{Type{CTScanner{FBP,M} where M},GrayScaleLine{Float32}})
    Base.precompile(Tuple{Type{CTScanner{FBP,M} where M},GrayScaleLine{Float64}})
    Base.precompile(Tuple{Type{CTScanner{FBP,M} where M},GrayScalePyramid{Float32}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBP{ParallelBeamGeometry{Float32,DefaultTomograph},Marta.Filters.RamLak},Array{Float32,2}},GrayScaleLine{Float32}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBP{ParallelBeamGeometry{Float32,DefaultTomograph},Marta.Filters.RamLak},Array{Float32,2}},ImageParams{Float32}})
    Base.precompile(Tuple{typeof(calibrate_tomogram),CTScanner{FBP{ParallelBeamGeometry{Float64,DefaultTomograph},Marta.Filters.RamLak},Array{Float64,2}},ImageParams{Float64}})
    Base.precompile(Tuple{typeof(project_image),Some{CTImage{Array{Float32,2}}},Radon{FanBeamGeometry{Float32,DefaultTomograph}}})
    Base.precompile(Tuple{typeof(project_image),Some{CTImage{Array{Float32,2}}},Radon{ParallelBeamGeometry{Float32,DefaultTomograph}}})
    Base.precompile(Tuple{typeof(project_image),Some{CTImage{Array{Float64,2}}},Radon{FanBeamGeometry{Float64,DefaultTomograph}}})
    Base.precompile(Tuple{typeof(project_image),Some{CTImage{Array{Float64,2}}},Radon{ParallelBeamGeometry{Float64,DefaultTomograph}}})
    isdefined(Marta.AbstractAlgorithms, Symbol("#14#16")) && Base.precompile(Tuple{getfield(Marta.AbstractAlgorithms, Symbol("#14#16")),CTSinogram{Array{Float32,2}}})
    isdefined(Marta.AbstractAlgorithms, Symbol("#14#16")) && Base.precompile(Tuple{getfield(Marta.AbstractAlgorithms, Symbol("#14#16")),CTSinogram{Array{Float64,2}}})
    let fbody = try __lookup_kwbody__(which(calibrate_tomogram, (CTScanner{FBP{ParallelBeamGeometry{Float32,DefaultTomograph},Marta.Filters.RamLak},Array{Float32,2}},GrayScaleLine{Float32},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Iterators.Pairs{Union{},Union{},Tuple{},NamedTuple{(),Tuple{}}},typeof(calibrate_tomogram),CTScanner{FBP{ParallelBeamGeometry{Float32,DefaultTomograph},Marta.Filters.RamLak},Array{Float32,2}},GrayScaleLine{Float32},))
        end
    end
end
