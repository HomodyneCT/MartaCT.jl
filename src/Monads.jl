module Monads

export Monad, Maybe, Optional
export mreturn, mjoin, mbind, mmap
export maybe, ↣

using SimpleTraits
using ..Applicative: Callable
import Base: something, isnothing, convert

@traitdef Monad{T}

function mbind end
function mreturn end
function mjoin end

@traitfn @inline mreturn(::Type{M}, x) where {M; Monad{M}} = M(x)
@traitfn @inline mjoin(m::M) where {M; Monad{M}} = mbind(identity, m)
@traitfn @inline mbind(f::F, m::M) where {F,M; Callable{F},Monad{M}} =
    f(mjoin(m))
@traitfn @inline function mmap(f::F, m::M) where {F,M; Callable{F},Monad{M}}
    mbind(m) do x
        mreturn(M, f(x))
    end
end
@traitfn @inline ↣(m::M, f::F) where {M,F; Monad{M},Callable{F}} = mbind(f, m)

const Optional{T} = Union{Nothing,T}
const Maybe{T} = Optional{Some{T}}

@traitimpl Monad{Maybe}

Maybe{T}(x::T) where {T} = Some(x)
Maybe{T}(::Nothing) where {T} = nothing
Maybe(::Nothing) = nothing
Maybe(x::T) where {T} = Maybe{T}(x)
Maybe(m::Maybe) = m

@traitfn @inline mbind(f::F, m::Maybe) where {F; Callable{F}} =
    isnothing(m) ? nothing : f(something(m))

@inline maybe(b, m::Optional) = isnothing(m) ? b : m
@inline maybe(b) = m::Optional -> maybe(b, m)
@traitfn @inline maybe(f::F, b, m::Optional) where {F; Callable{F}} =
    maybe(b, Maybe(m) ↣ f)

end # module
