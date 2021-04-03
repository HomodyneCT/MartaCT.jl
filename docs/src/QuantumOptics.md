# Quantum Optics example

```@setup qoptics
ENV["GKS_WSTYPE"] = "svg"
using Marta, Plots, QuantumOptics, IntervalSets, LinearAlgebra, Plots.Measures
gr()
Plots.reset_defaults()
default(size=(500,300), rightmargin=1cm)
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
nothing # hide
```

We need to get the Wigner function representation in order
to compute the marginal distributions of the position.

```@example qoptics
ζ = 10
xs = linspace(-ζ..ζ, 200)
W = wigner(ρ, xs, xs)
heatmap(xs, xs, W')
savefig("Wcat.svg"); nothing # hide
```

![](Wcat.svg)

The marginal distributions can be computed as a Radon
transform of the Wigner function:

```@example qoptics
ϕs = linspace(ORI(0..2π), 800)
marg = radon(W, xs, ϕs, RadonSquare())
heatmap(ϕs, xs, marg)
savefig("marg.svg"); nothing # hide
```

![](marg.svg)

Now we can employ the standard FBP algorithm to recover the
Wigner distribution:

```@example qoptics
Wrec = iradon(W, xs, xs, FBPFFTSquare())
heatmap(xs, xs, Wrec)
savefig("Wrec.svg"); nothing # hide
```

![](Wrec.svg)
