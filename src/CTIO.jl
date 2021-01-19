module CTIO

Base.Experimental.@optlevel 0

export load_data, sinog_extract, sinog_load
export load_yaml, write_yaml, struct2dict, yaml_repr
export write_ct_image, read_ct_image
export load_image, write_image
export load_sinogram, write_sinogram
export load_tomogram, write_tomogram

import YAML
using Mmap: mmap
using IterTools: imap

include("TypeDict.jl")
using .TypeDict: standardize_type


struct DataFileHelper
    header::Int
    size::Int
    width::Int
    height::Int
    tot_len::Int
    offset::Int
    dtype::DataType
    stype::DataType
end


function DataFileHelper(kwargs...)
    arg_dict = Dict(kwargs)
    header = get(arg_dict, :header, 0)
    width = get(arg_dict, :width, 0)
    tot_len = get(arg_dict, :length, 0)
    @assert (width != 0 || tot_len != 0) "Specify at least either width or length"
    width != 0 && tot_len == 0 && (tot_len = width)
    height = get(arg_dict, :height, 0)
    offset = get(arg_dict, :offset, 0)
    dtype = standardize_type(get(arg_dict, :dtype, UInt16))
    stype = standardize_type(get(arg_dict, :stype, Float32))
    size = tot_len * height
    DataFileHelper(header, size, width, height, tot_len, offset, dtype, stype)
end


function load_data(
    io::IO;
    header::Int = 0,
    dtype::Union{String,DataType} = UInt16,
)
    if header != 0
        seek(io, header)
    end
    mmap(io, Vector{standardize_type(dtype)})
end


function load_data(f::AbstractString; kwargs...)
    open(f) do s
        load_data(s; kwargs...)
    end
end


function sinog_extract(data::Vector{<:Real}, dhp::DataFileHelper)
    if dhp.size == 0
        eff_size = length(data)
    else
        eff_size = dhp.size
    end
    width = dhp.width
    height = dhp.height
    if height == 0
        height = eff_size ÷ dhp.tot_len
        eff_size = width * height
    end
    conv_data = convert(Vector{dhp.stype}, data[1:eff_size])
    reshape(conv_data, height, width)
end


sinog_extract(data; kwargs...) = sinog_extract(data, DataFileHelper(kwargs...))


function sinog_load(s::IO; kwargs...)
    dhp = DataFileHelper(kwargs...)
    raw_data = load_data(s, header = dhp.header, dtype = dhp.dtype)
    sinog_extract(raw_data, dhp)
end


function sinog_load(f::AbstractString; kwargs...)
    open(f) do s
        sinog_load(s; kwargs...)
    end
end


function load_yaml(s::IO)
    yaml_doc = YAML.load(s; dicttype = Dict{Symbol,Any})
    (; pairs(yaml_doc)...)
end


function load_yaml(f::AbstractString)
    open(f) do s
        load_yaml(s)
    end
end


yaml_repr(obj) = obj
yaml_repr(ntup::NamedTuple) = struct2dict(ntup)


function struct2dict(obj)
    Dict(map(propertynames(obj)) do p
        field = getproperty(obj, p)
        p => yaml_repr(field)
    end)
end


function write_yaml(io::IO, obj)
    YAML.write(io, yaml_repr(obj))
end


function write_yaml(io::IO; kwargs...)
    write_yaml(io, Dict(imap(yaml_repr, kwargs)))
end


function write_yaml(f::String, obj)
    open(f, "w") do s
        write_yaml(s, obj)
    end
end


function write_yaml(f::String; kwargs...)
    open(f, "w") do s
        write_yaml(s; kwargs...)
    end
end


"""
    read_ct_image(io::IO[, T=Float32]; kwargs...) where {T<:Real}

Read CT image from the stream `io`.

# Keyword arguments

- `header`: read rows and columns from header in the file. [default: false]
- `rows`: the image rows.
- `cols`: the image columns.
- `header_type`: type for the two values of the header. [default: `Int64`]
- `row_major`: assume the image is stored in row major ordering.

See also: [`write_ct_image`](@ref)
"""
function read_ct_image(
    io::IO,
    ::Type{T} = Float32;
    header = false,
    rows = 0,
    cols = 0,
    header_type = Int64,
    row_major = false,
) where {T<:Real}
    if header
        rows = read(io, header_type)
        cols = read(io, header_type)
    elseif cols == 0
        cols = rows
    end
    if row_major
        data = Matrix{T}(undef, cols, rows)
        return read!(io, data)' |> collect
    end
    data = Matrix{T}(undef, rows, cols)
    read!(io, data)
end



"""
    read_ct_image(f::AbstractString[, T=Float32]; <keyword arguments>) where {T<:Real}

Read CT image from file `f`.
"""
function read_ct_image(
    f::AbstractString,
    ::Type{T} = Float32;
    kwargs...,
) where {T<:Real}
    open(f) do s
        read_ct_image(s, T; kwargs...)
    end
end



"""
    write_ct_image(io::IO, image::AbstractMatrix{T}; kwargs...) where {T<:Real}

Write CT image to the stream `io`.

# Keyword arguments

- `header`: if `true`, write header information with rows and columns. [default: false]
- `header_type`: type for the two values of the header. [default: `Int64`]
- `row_major`: write image in row major ordering.

See also: [`read_ct_image`](@ref)
"""
function write_ct_image(
    io::IO,
    image::AbstractMatrix{T};
    header = false,
    header_type = Int64,
    row_major = false,
) where {T<:Real}
    if header
        write(io, header_type(rows))
        write(io, header_type(cols))
    end
    if row_major
        write(io, image')
        return
    end
    write(io, image)
end



"""
    write_ct_image(f::AbstractString, image::AbstractMatrix{T}; <keyword arguments>) where {T<:Real}

Write CT image to file `f`.
"""
function write_ct_image(
    f::AbstractString,
    image::AbstractMatrix{T};
    kwargs...,
) where {T<:Real}
    open(f, "w") do s
        write_ct_image(s, image; kwargs...)
    end
end



"""
    load_image(io::IO[, T=Float32]; <keyword arguments>) where {T<:Real}

Load image from stream `io`.

# Arguments
- `io`: the stream where to load the image from.
- `header=false`: specify if the file contains the header.
- `nϕ=0`: number of projections (rows in the file for row major ordering).
- `nd=0`: number of detectors (columns in the file for row major ordering).
- `row_major=false`: tells if the data is stored in row major ordering or not.
- `header_type=Int64`: the type of the header data.

See Also: [`write_image`](@ref), [`read_ct_image`](@ref), [`write_ct_image`](@ref)
"""
function load_image(
    io::IO,
    ::Type{T} = Float32;
    header = false,
    rows = 0,
    cols = 0,
    row_major = false,
    header_type = Int64,
) where {T<:Real}
    read_ct_image(
        io,
        T;
        header = header,
        rows = rows,
        cols = cols,
        row_major = false,
        header_type = header_type,
    )
end



"""
    load_image(f::AbstractString[, T=Float32]; <keyword arguments>) where {T<:Real}

Load image from file `f`.
"""
function load_image(
    f::AbstractString,
    ::Type{T} = Float32;
    kwargs...,
) where {T<:Real}
    open(f) do s
        load_image(s, T; kwargs...)
    end
end



"""
    write_image(io::IO, image::AbstractMatrix{T}; <keyword arguments>) where {T<:Real}

Write input image to stream `io`.

# Arguments
- `io`: the stream where the image should be written to.
- `image`: the image data matrix.
- `header`: specify if the image header should be written.
- `row_major`: specify if the image should be stored in row major order.
- `header_type=Int64`: the type of the header data.

See Also: [`read_image`](@ref), [`write_ct_image`](@ref), [`read_ct_image`](@ref).
"""
function write_image(
    io::IO,
    image::AbstractMatrix{T};
    header = false,
    row_major = false,
    header_type = Int64,
) where {T<:Real}
    write_ct_image(
        io,
        image;
        header = header,
        row_major = row_major,
        header_type = header_type,
    )
end


"""
    write_image(f::AbstractString, image::AbstractMatrix{T}; <keyword arguments>) where {T<:Real}

Write input image to file `f`.
"""
function write_image(
    f::AbstractString,
    image::AbstractMatrix{T};
    kwargs...,
) where {T<:Real}
    open(f, "w") do s
        write_image(s, image; kwargs...)
    end
end


"""
    load_sinogram(io::IO[, ::Type{T}=Float32]; <keyword arguments>) where {T<:Real}

Load sinogram from stream `io`.

This function is specifically designed to load a sinogram. Do not
expect consistent result with other type of images. Use more general
function instead, refer to the `See Also` section.

# Arguments
- `io`: the stream where to load the image from.
- `header=false`: specify if the file contains the header.
- `nϕ=0`: number of projections (rows in the file for row major ordering).
- `nd=0`: number of detectors (columns in the file for row major ordering).
- `row_major=false`: tells if the data is stored in row major ordering or not.
- `header_type=Int64`: the type of the header data.

See Also: [`write_sinogram`](@ref), [`read_ct_image`](@ref), [`write_ct_image`](@ref)
"""
function load_sinogram(
    io::IO,
    ::Type{T} = Float32;
    header = false,
    nϕ = 0,
    nd = 0,
    row_major = false,
    header_type = Int64,
) where {T<:Real}
    read_ct_image(
        io,
        T;
        header = header,
        rows = nd,
        cols = nϕ,
        row_major = false,
        header_type = header_type,
    )
end



"""
    load_sinogram(f::AbstractString[, T=Float32]; <keyword arguments>) where {T<:Real}

Load sinogram from file `f`.
"""
function load_sinogram(
    f::AbstractString,
    ::Type{T} = Float32;
    kwargs...,
) where {T<:Real}
    open(f) do s
        load_sinogram(s, T; kwargs...)
    end
end



"""
    write_sinogram(io::IO, sinog::AbstractMatrix{T}; <keyword arguments>) where {T<:Real}

Write sinogram `sinog` to stream `io`.

This function is specifically designed to deal with a sinogram image.

# Arguments
- `io`: the stream where to load the image from.
- `header=false`: specify if the file contains the header.
- `row_major=false`: tells if the data is stored in row major ordering or not.
- `header_type=Int64`: the type of the header data.

See Also: [`load_sinogram`](@ref), [`read_ct_image`](@ref), [`write_ct_image`](@ref).
"""
function write_sinogram(
    io::IO,
    sinog::AbstractMatrix{T};
    header = false,
    row_major = false,
    header_type = Int64,
) where {T<:Real}
    write_ct_image(
        io,
        sinog;
        header = header,
        row_major = row_major,
        header_type = header_type,
    )
end



"""
    write_sinogram(f::AbstractString, sinog::AbstractMatrix{T}) where {T<:Real}

Write sinogram to file `f`.
"""
function write_sinogram(
    f::AbstractString,
    sinog::AbstractMatrix{T};
    kwargs...,
) where {T<:Real}
    open(f, "w") do s
        write_sinogram(s, sinog; kwargs...)
    end
end



"""
    load_tomogram(io::IO[, T=Float32]; <keyword arguments>) where {T<:Real}

Load reconstructed image from stream `io`.

This function is specifically designed to load a tomogram. Do not
expect consistent result with other type of images. Use more general
function instead, refer to the `See Also` section.

# Arguments
- `io`: the stream where to load the image from.
- `header=false`: specify if the file contains the header.
- `nϕ=0`: number of projections (rows in the file for row major ordering).
- `nd=0`: number of detectors (columns in the file for row major ordering).
- `row_major=false`: tells if the data is stored in row major ordering or not.
- `header_type=Int64`: the type of the header data.

See Also: [`write_tomogram`](@ref), [`read_ct_image`](@ref), [`write_ct_image`](@ref)
"""
function load_tomogram(
    io::IO,
    ::Type{T} = Float32;
    header = false,
    rows = 0,
    cols = 0,
    row_major = false,
    header_type = Int64,
) where {T<:Real}
    load_image(
        io,
        T;
        header = header,
        rows = rows,
        cols = cols,
        row_major = row_major,
        header_type = header_type,
    )
end



"""
    load_tomogram(f::AbstractString[, T=Float32]; <keyword arguments>) where {T<:Real}

Load reconstructed image from file `f`.
"""
function load_tomogram(
    f::AbstractString,
    ::Type{T} = Float32;
    kwargs...,
) where {T<:Real}
    open(f) do s
        load_tomogram(s, T; kwargs...)
    end
end



"""
    write_tomogram(io::IO, tomog::AbstractMatrix{T}; <keyword arguments>) where {T<:Real})

Write reconstructed image to stream `io`.

# Arguments
- `io`: the stream where to load the image from.
- `header=false`: specify if the file contains the header.
- `row_major=false`: tells if the data is stored in row major ordering or not.
- `header_type=Int64`: the type of the header data.
"""
function write_tomogram(
    io::IO,
    tomog::AbstractMatrix{T};
    header = false,
    row_major = false,
    header_type = Int64,
) where {T<:Real}
    write_image(
        io,
        tomog;
        header = header,
        row_major = row_major,
        header_type = header_type,
    )
end


"""
    write_tomogram(f::AbstractString, tomog::AbstractMatrix{T}; <keyword arguments>) where {T<:Real}

Write reconstructed image to file `f`.
"""
function write_tomogram(
    f::AbstractString,
    tomog::AbstractMatrix{T};
    kwargs...,
) where {T<:Real}
    open(f, "w") do s
        write_tomogram(s, tomog; kwargs...)
    end
end

end  # module
