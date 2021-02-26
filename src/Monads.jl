module Monads

export Monad, Maybe, Optional
export mreturn, mjoin, mbind, mmap
export maybe, ↣

import Base: something, isnothing, convert

function mbind end
function mreturn end
function mjoin end

abstract type Monad end

const Optional{T} = Union{Nothing,T}
const Maybe{T} = Optional{Some{T}}
const MonadOrMaybe = Union{Monad,Maybe}
const FuncOrType = Union{Function,Type}


mreturn(::Type{M}, x) where {M<:MonadOrMaybe} = M(x)
mjoin(m::MonadOrMaybe) = mbind(identity, m)
mmap(f::FuncOrType, m::M) where {M<:MonadOrMaybe} =
    mbind(m) do x
        mreturn(M, f(x))
    end
↣(m::MonadOrMaybe, f::FuncOrType) = mbind(f, m)


# struct Maybe{T} <: Monad
#     value::Optional{Some{T}}
#     Maybe{T}(::Nothing) where {T} = new(nothing)
#     Maybe{T}(s::Some{T}) where {T} = new(s)
#     Maybe{T}(x::T) where {T} = Maybe{T}(Some(x))
# end

Maybe{T}(x::T) where {T} = Some(x)
Maybe{T}(::Nothing) where {T} = nothing
#Maybe{T}(::Maybe{Nothing}) where {T} = Maybe{T}(nothing)
#Maybe(::Nothing) = Maybe{Nothing}(nothing)
Maybe(::Nothing) = nothing
#Maybe(x::T) where {T} = Maybe{T}(Some(x))
Maybe(x::T) where {T} = Maybe{T}(x)
Maybe(m::Maybe) = m

# We cannot export this method, otherwise bad things happen to other's code.
#convert(::Type{<:Maybe{T}}, s::S) where {T,S} = Maybe(convert(T,s))

#something(m::Maybe) = something(m.value)
#isnothing(m::Maybe) = isnothing(m.value)

mbind(f::FuncOrType, m::Maybe) =
    isnothing(m) ? nothing : f(something(m))

maybe(b, m::Optional) = isnothing(m) ? b : m
#maybe(b, m::Maybe) = maybe(b, something(m))
#maybe(b) = m::Union{Maybe,Optional} -> maybe(b, m)
maybe(b) = m::Optional -> maybe(b, m)
# maybe(f::Function, b, m::Union{Maybe,Optional}) =
#     maybe(b, mbind(f, m))
#maybe(f::Function, b, m::Optional) = maybe(b, mbind(f, Maybe(m)))
maybe(f::Function, b, m::Optional) = f(maybe(b, m)) # Not sure which one is correct, this works for FBP
#maybe(f::Function, b) = maybe(b) ∘ f

end # module
