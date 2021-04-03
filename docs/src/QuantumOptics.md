# Quantum Optics example

Marta can also be used for Quantum Tomography as the basic
_filtered backprojection_ algorithm can be used to
reconstruct a quantum state (or more precisely the Wigner
function) from the measurements of the position.

First we need to import the packages for this example. We
use the QuantumOptics Julia package for the definition of
the quantum objects.

```@example qoptics
using Marta, Plots, QuantumOptics, IntervalSets
gc() # hide
Plots.reset_defaults() # hide
default(size=(500,300)); nothing # hide
```
