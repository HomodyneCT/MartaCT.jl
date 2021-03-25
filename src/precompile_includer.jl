should_precompile =   true


# Don't edit the following! Instead change the script for `snoop_bot`.
ismultios = true
ismultiversion = true
# precompile_enclosure
@static if !should_precompile
    # nothing
elseif !ismultios && !ismultiversion
    @static if isfile(joinpath(@__DIR__, "../deps/SnoopCompile/precompile/precompile_Marta.jl"))
        include("../deps/SnoopCompile/precompile/precompile_Marta.jl")
        _precompile_()
    end
else
    @static if Sys.islinux() 
    @static if v"1.5.0-DEV" <= VERSION <= v"1.5.9" 
    @static if isfile(joinpath(@__DIR__, "../deps/SnoopCompile/precompile/linux/1.5/precompile_Marta.jl"))
        include("../deps/SnoopCompile/precompile/linux/1.5/precompile_Marta.jl")
        _precompile_()
    end
else 
    @static if isfile(joinpath(@__DIR__, "../deps/SnoopCompile/precompile/linux/1.5/precompile_Marta.jl"))
        include("../deps/SnoopCompile/precompile/linux/1.5/precompile_Marta.jl")
        _precompile_()
    end
    end

elseif Sys.isapple() 
    @static if v"1.5.0-DEV" <= VERSION <= v"1.5.9" 
    @static if isfile(joinpath(@__DIR__, "../deps/SnoopCompile/precompile/apple/1.5/precompile_Marta.jl"))
        include("../deps/SnoopCompile/precompile/apple/1.5/precompile_Marta.jl")
        _precompile_()
    end
else 
    @static if isfile(joinpath(@__DIR__, "../deps/SnoopCompile/precompile/apple/1.5/precompile_Marta.jl"))
        include("../deps/SnoopCompile/precompile/apple/1.5/precompile_Marta.jl")
        _precompile_()
    end
    end

elseif Sys.iswindows() 
    @static if v"1.5.0-DEV" <= VERSION <= v"1.5.9" 
    @static if isfile(joinpath(@__DIR__, "../deps/SnoopCompile/precompile/windows/1.5/precompile_Marta.jl"))
        include("../deps/SnoopCompile/precompile/windows/1.5/precompile_Marta.jl")
        _precompile_()
    end
else 
    @static if isfile(joinpath(@__DIR__, "../deps/SnoopCompile/precompile/windows/1.5/precompile_Marta.jl"))
        include("../deps/SnoopCompile/precompile/windows/1.5/precompile_Marta.jl")
        _precompile_()
    end
    end

else 
    @static if v"1.5.0-DEV" <= VERSION <= v"1.5.9" 
    @static if isfile(joinpath(@__DIR__, "../deps/SnoopCompile/precompile/linux/1.5/precompile_Marta.jl"))
        include("../deps/SnoopCompile/precompile/linux/1.5/precompile_Marta.jl")
        _precompile_()
    end
else 
    @static if isfile(joinpath(@__DIR__, "../deps/SnoopCompile/precompile/linux/1.5/precompile_Marta.jl"))
        include("../deps/SnoopCompile/precompile/linux/1.5/precompile_Marta.jl")
        _precompile_()
    end
    end

    end

end # precompile_enclosure
