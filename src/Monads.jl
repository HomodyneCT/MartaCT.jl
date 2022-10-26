module Monads



export Monad, Optional
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

const Optional{T} = Union{T,Nothing}

@inline maybe(b, m::Optional) = isnothing(m) ? b : m
@inline maybe(b) = m::Optional -> maybe(b, m)
@traitfn @inline maybe(f::F, b, m::Optional) where {F; Callable{F}} =
    maybe(b, Maybe(m) ↣ f)

end # module
