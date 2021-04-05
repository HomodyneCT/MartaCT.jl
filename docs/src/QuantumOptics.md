# Quantum Optics example

```@setup qoptics
ENV["GKS_WSTYPE"] = "svg"
using Marta, Plots, QuantumOptics, IntervalSets, LinearAlgebra, Plots.Measures
gr()
Plots.reset_defaults()
default(size=(500,300), rightmargin=0.5cm)
```

Marta can also be used for Quantum Tomography as the basic
_filtered backprojection_ algorithm can be used to
reconstruct a quantum state (or more precisely the Wigner
function) from the measurements of the position.

First we need to import the packages for this example. We
use the QuantumOptics Julia package for the definition of
the quantum objects.

```julia
using Marta, Plots, QuantumOptics, IntervalSets, LinearAlgebra
```

At this point, let's construct an optical cat state:

```@example qoptics
N = 64 # density matrix dimension
α = 5
bs = FockBasis(N-1)
ν = inv(√(2*(1+exp(-2abs2(α)))))
ψ = ν * (coherentstate(bs, α) + coherentstate(bs, -α))
ρ = ψ ⊗ ψ'
@show tr(ρ)
nothing # hide
```

```@setup qoptics
bar(diag(ρ.data)|>real; c=:red, leg=:none)
savefig("rho-diag.svg"); nothing # hide
```

![](rho-diag.svg)

We need to get the Wigner function representation in order
to compute the marginal distributions of the position.

```@example qoptics
ζ = 10
xs = linspace(-ζ..ζ, 200)
W = wigner(ρ, xs, xs) |> permutedims |> CTTomogram
δW = 4ζ^2 / length(W)
@show sum(W) * δW
heatmap(xs, xs, W)
savefig("Wcat.svg"); nothing # hide
```

![](Wcat.svg)

The marginal distributions can be computed as a Radon
transform of the Wigner function:

```@example qoptics
ϕs = linspace(ORI(0..2π), 800)
marg = radon(W, xs, ϕs, RadonSquare())
@show sum(marg[:,1])
heatmap(ϕs, xs, marg)
savefig("marg.svg"); nothing # hide
```

![](marg.svg)

Now we can employ the standard FBP algorithm to recover the
Wigner distribution:

```@example qoptics
marg′ = marg * length
Wrec = iradon(marg, xs, xs, FBPFFTSquare())
δW = 2ζ^3 / π * (length(xs)-1) / length(ϕs) / length(Wrec)
@show sum(Wrec) * δW
heatmap(xs, xs, Wrec)
savefig("Wrec.svg"); nothing # hide
```

![](Wrec.svg)

We can check that the normalization is preserved:

```@example qoptics
δWrec = 4ζ^2 / (length(xs)-1)^2
@show sum(Wrec) * δWrec
nothing # hide
```
