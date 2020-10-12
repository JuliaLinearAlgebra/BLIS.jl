BLIS.jl
=======

This repository to provides:

- Wrapper for typed interface of [BLIS](https://github.com/flame/blis).
- As BLIS itself is using actively templates and macros, 
  this package also aims to maximize usage of Julia's
  metaprogramming features.
  

## Installation
```
]add https://github.com/xrq-phys/BLIS.jl
```

## Usage

### LinearAlgebra Frontend

The first thing BLIS.jl provides is a linear algebra frontend which resembles BLAS functions very much akin to `LinearAlgebra.BLAS.gemm!`, but defined for both column-major and generic strided matrices (namely all `strides(A)[1] >= 1`). An example for direct computation of generic strided case could be:

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

Mixed precision is also directly supported by this interface.

As BLAS interface defined by Julia's `LinearAlgebra.BLAS` is a few routines beyond standard BLAS, [here is a list for status of support](src/interface_linalg/ABOUT.md).

### BLIS Typed Backend

[Typed-API](https://github.com/flame/blis/blob/master/docs/BLISTypedAPI.md) backend of BLIS is directly accesible via the `BLIS.TypedBackend` submodule. The following example demonstrates usage of BLIS's [`bli_dmkherm`](https://github.com/flame/blis/blob/master/docs/BLISTypedAPI.md#mkherm) method:

```julia
using BLIS
# Import BLIS' enumerate types
# like BLIS_UPPER and BLIS_TRANSPOSE.
using BLIS.Types

A = rand(4, 4)
# NOTE: of course for rsa and csa parameters
# one can directly pass strides(A) and expand.
BLIS.TypedBackend.bli_dmkherm!(BLIS_UPPER, 4, A, strides(A)...)
A # Check that A is projected to be Hermitian.
```

### BLIS Object Backend

[Object-API](https://github.com/flame/blis/blob/master/docs/BLISObjectAPI.md) backend is exposed as `BLIS.ObjectBackend` submodule. Mixed-precision BLAS functions are only available via Object-API methods:

```julia
using BLIS
using BLIS.ObjectBackend

# α and β should always be created from numbers.
oα = BliObj(1.0);
oβ = BliObj(1.0);
# Matrix objects can be created from Arrays or views.
oA = BliObj(rand(Float32, 4, 4));
oB = BliObj(rand(Float32, 4, 4));
C = zeros(Float64, 8, 8);
wC = view(C, 1:2:8, 1:2:8);
oC = BliObj(wC);
# Mixed-precision GEMM into scattered target storage:
#      A * B -> C[1:2:7, 1:2:7],
#  where A and B have single precision and multiply into
#  a double precision target matrix.
BLIS.ObjectBackend.bli_gemm!(oα, oA, oB, oβ, oC)
# View result at target matrix C.
# Here's an example output.
C
[
 1.59116  0.0  1.34804   0.0  1.68991  0.0  1.58798  0.0;
 0.0      0.0  0.0       0.0  0.0      0.0  0.0      0.0;
 1.27603  0.0  1.31893   0.0  1.40484  0.0  1.06297  0.0;
 0.0      0.0  0.0       0.0  0.0      0.0  0.0      0.0;
 2.12725  0.0  1.72495   0.0  2.37218  0.0  1.56937  0.0;
 0.0      0.0  0.0       0.0  0.0      0.0  0.0      0.0;
 1.0284   0.0  0.858046  0.0  1.05108  0.0  1.25478  0.0;
 0.0      0.0  0.0       0.0  0.0      0.0  0.0      0.0;
]
```

## Roadmaps

- Define interface methods like `Base.*` and `Base.mul!` for basic
  Julia types like `Array` and `StridedArray`.
- Provide option to set BLIS as BLAS provider.
- Incorporate [HPAC/Linnea](https://github.com/HPAC/linnea) in this
  or another repository.
- Introduce BLIS' testsuite and fallback to `LinearAlgebra.BLAS` for
  routines that failed the tests.
