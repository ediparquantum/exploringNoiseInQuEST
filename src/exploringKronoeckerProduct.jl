## Pkg.add("Distributions") # If you don't already have it installed.
using Distributions
A = rand(1:10, 3,3)
B = rand(1:10, 3,3)
C = rand(1:10, 3,3)
μ = rand(1:100)
μA = μ .* A
μB = μ .* B

kron(μA,B) == kron(A,μB) == μ*kron(A,B)

kron(A+B,C) == (kron(A,C) + kron(B,C))


