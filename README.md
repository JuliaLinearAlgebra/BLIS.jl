BLIS.jl
=======

This repository to provides:

- Wrapper for typed interface of [BLIS](https://github.com/flame/blis).
- As BLIS itself is using actively templates and macros, 
  this package also aims to maximize usage of Julia's
  metaprogramming features.
  

## Installation
```julia
]add https://github.com/xrq-phys/BLIS.jl
```

## Usage

### LinearAlgebra Frontend

**WIP** The first thing BLIS.jl provides is a linear algebra frontend which resembles BLAS functions very much akin to `LinearAlgebra.BLAS.gemm!`, but defined for both column-major and generic strided matrices (namely all `strides(A)[1] >= 1`). An example for direct computation of generic strided case could be:

```julia
using BLIS
using BLIS.BLASInterface
A = rand(ComplexF64, 4, 4);
B = rand(ComplexF64, 4, 4);
C = ones(ComplexF64, 4, 4);
# Create matrix views. No memory copied!
wA = view(A, 1:2:4, 1:2:4);
wB = view(B, 1:2:4, 1:2:4);
wC = view(C, 1:2:4, 1:2:4);
# This call then computes 2×2×2 product:
#  wA*wB -> wC
# directly without reallocating wA, wB or wC to other arrays.
BLASInterface.gemm!('N', 'N', 1.0+0im, wA, wB, 1.0+0im, wC)
```

## Roadmaps

- Provide support also for BLIS' object interface.
- Implement `LinearAlgebra.BLAS` frontend against the object interface to enable mixed precision support.
- Define all `LinearAlgebra.BLAS`-compatible frontends.
- Define interface methods like `Base.*` and `Base.mul!` for basic
  Julia types like `Array` and `StridedArray`.
- Provide option to set BLIS as BLAS provider.
- Incorporate [HPAC/Linnea](https://github.com/HPAC/linnea) in this
  or another repository.
- Introduce BLIS' testsuite and fallback to `LinearAlgebra.BLAS` for
  routines that failed the tests.
