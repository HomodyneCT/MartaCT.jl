module Applicative

export Callable
export apply, <|
export fimap, fmap, iflip, trpos, make_indices

using SimpleTraits
using IterTools: imap

@traitdef Callable{F}

@traitimpl Callable{Function}
@traitimpl Callable{Type}

apply(f, x) = f(x)
<|(f, x) = f(x)

fimap(f) = a -> imap(f, a)
fmap(f) = a -> map(f, a)

iflip(x, y) = (y, x)
iflip((x, y)) = (y, x)
iflip(i::CartesianIndex{2}) = Tuple(i) |> iflip |> CartesianIndex

trpos(t) = p -> CartesianIndex(Tuple(p) .+ t)

function make_indices(rows, cols, itr)
    (
        fimap(trpos((rows ÷ 2 + 1, cols ÷ 2 + 1))) ∘ fimap(iflip) ∘
        fimap(x -> CartesianIndex(round.(Int, Tuple(x))))
    )(
        itr,
    )
end

end # module
