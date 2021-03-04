using IterTools: repeatedly

@deprecate resample_sinogram sample_sinogram_external false


"""generate_photons(n::Integer, nx::Integer, nϕ::Integer)

Generate ``⟨n⟩`` Poisson distributed photons per projection
angle over an array of `nx` detectors.

Return a `nd × nϕ` matrix with the generated photons.
```
"""
@inline function generate_photons(n::Integer, nx::Integer, nϕ::Integer)
    d = Distributions.Poisson(n)
    emit = Distributions.DiscreteUniform(1,nx)
    rngs_gen = [Random.MersenneTwister() for _ ∈ 1:Threads.nthreads()]
    rngs_measure = [Random.MersenneTwister() for _ ∈ 1:Threads.nthreads()]
    photons = zeros(Int, nx, nϕ)
    p = Progress(nϕ, 0.2, "Generating photons: ")
    Threads.@threads for i ∈ 1:nϕ
        id = Threads.threadid()
        @inbounds g = rngs_gen[id]
        @inbounds m = rngs_measure[id]
        num = rand(g, d)
        foreach(1:num) do _
            x = rand(m, emit)
            @inbounds photons[x, i] += 1
        end
        next!(p)
    end
    photons
end


"""sample_sinogram(sinog::AbstractMatrix{T}; <keyword arguments>) where {T <: Real}

Simulate a low dose CT scan. This samples `sinog` with ``⟨n⟩``
random photons per projection angle.

# Arguments
- `sinog`: sinogram data.
- `nphotons::Integer=10000`: mean number of photons.
- `ϵ::Real=1`: detectors quantum efficiency.
- `take_log::Bool=true`: whether to take the logarithm of the
  resampled intensities to obtain the corresponding
  sinogram.
"""
function sample_sinogram(
    sinog::AbstractMatrix{T};
    nphotons::Integer = 10000,
    ϵ::Real = 1,
    take_log::Bool = true,
) where {T<:Real}
    nd, nϕ = size(sinog)
    photons = generate_photons(nphotons, nd, nϕ)
    rngs_absorp = [Random.MersenneTwister() for _ ∈ 1:Threads.nthreads()]
    rngs_clicks = [Random.MersenneTwister() for _ ∈ 1:Threads.nthreads()]
    measure = Distributions.Uniform()
    detector = Distributions.Uniform()
    low_sample = similar(sinog, nd, nϕ)
    p = Progress(length(low_sample), 0.2, "Tracking photons: ")
    Threads.@threads for k ∈ eachindex(low_sample)
        id = Threads.threadid()
        @inbounds a = rngs_absorp[id]
        @inbounds c = rngs_clicks[id]
        @inbounds u::T = exp(-sinog[k])
        @inbounds n = photons[k]
        @inbounds low_sample[k] = count(1:n) do _
            rand(a, measure) <= u && rand(c, detector) <= ϵ
        end
        next!(p)
    end
    max_photons = maximum(low_sample)
    if take_log
        smax = maximum(sinog)
        @. low_sample = -log(low_sample / max_photons)
    else
        @. low_sample = low_sample / max_photons
    end
    CTSinogram(low_sample)
end

sample_sinogram(nphotons::Integer; kwargs...) =
    x -> sample_sinogram(x; nphotons, kwargs...)
sample_sinogram(; kwargs...) = x -> sample_sinogram(x; kwargs...)

function Random.rand(sinog::CTSinogram; kwargs...)
    sinog ↣ sample_sinogram(; kwargs...)
end

function Random.rand(sinog::CTSinogram, nphotons::Integer; kwargs...)
    rand(sinog; nphotons, kwargs...)
end

"""
    sample_sinogram_external(sinog::AbstractMatrix{T}; <keyword arguments>) where {T <: Real}

Simulate a low dose CT scan. This Resample `sinog` with `n` photons.

# Arguments
- `sinog`: sinogram data.
- `sinog_path='_tmp_sinog.dat'`: path to a file where to write the sinogram to.
- `resampled_path='_tmp_resampled.dat'`: path to a file where to read the resampled sinogram from.
- `n=10000`: mean number of photons.
- `ϵ=1`: detectors quantum efficiency.
- `take_log=true`: whether to take the logarithm of the resampled intensities to obtain the corresponding sinogram.
- `verbosity=0`: set verbosity level.
- `progress=false`: show progress bar.
- `options=[]`: Additional options to be passed to the program.
"""
function sample_sinogram_external(
    sinog::AbstractMatrix{T};
    sinog_path = "_tmp_sinog.dat",
    resampled_path = "_tmp_resampled.dat",
    log_path = "_tmp_log.txt",
    n = 10000,
    ϵ = 1,
    take_log = true,
    verbosity = nothing,
    progress = false,
    options = [],
) where {T<:Real}
    tmp_path = mk_temp_dir()
    sinog_path = joinpath(tmp_path, sinog_path)
    resampled_path = joinpath(tmp_path, resampled_path)
    log_path = joinpath(tmp_path, log_path)
    vlevel = verbosity_level(verbosity)
    write_ct_image(sinog_path, sinog)
    nd, nϕ = size(sinog)
    cmd_args = []
    if !isnothing(vlevel)
        append!(cmd_args, ["$vlevel"])
    end
    if !progress
        append!(cmd_args, ["--xpbar"])
    end
    append!(cmd_args, options)
    cmd = @cmd "ctrun resample -ILT1 -D \$nd -P \$nϕ -n \$n -e \$ϵ \$cmd_args"
    run(pipeline(
        cmd,
        stdin = sinog_path,
        stdout = resampled_path,
        stderr = log_path,
    ))
    read_ct_image(resampled_path; rows = nd, cols = nϕ)
end
