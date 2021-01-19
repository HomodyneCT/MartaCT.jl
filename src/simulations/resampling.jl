"""
    resample_sinogram(sinog::AbstractMatrix{T}; <keyword arguments>) where {T <: Real}

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

See also: [`compute_radon`](@ref), [`compute_iradon`](@ref), [`compute_qct`](@ref)
"""
function resample_sinogram(
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


resample_sinogram(; kwargs...) =
    x -> CTSinogram <| resample_sinogram(x; kwargs...)
