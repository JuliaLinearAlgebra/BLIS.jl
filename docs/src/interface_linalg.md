# LinearAlgebra Frontend

The first thing BLIS.jl provides is a linear algebra frontend which resembles BLAS functions very much akin to `LinearAlgebra.BLAS.gemm!`, but defined for both column-major and generic strided matrices (namely all `strides(A)[1] >= 1`). An example for direct computation of generic strided case could be:

```julia
using BLIS

A = rand(ComplexF64, 4, 4);
B = rand(ComplexF64, 4, 4);
C = ones(ComplexF64, 4, 4);
# Create matrix views. No memory copied!
wA = view(A, 1:2:4, 1:2:4);
wB = view(B, 1:2:4, 1:2:4);
wC = view(C, 1:2:4, 1:2:4);
# wA * wB' can be directly computed via:
sC = wA * wB';
# without reallocating wA or wB to other arrays as
# mul! is overriden to call BLIS within this module.

# One can also use the same gemm! as provided in
# LinearAlgebra.BLAS to compute this 2×2×2 product:
#  wA*wB -> wC
BLIS.BLASInterface.gemm!('N', 'N', 1.0+0im, wA, wB, 1.0+0im, wC)
# Note that storage target can also be generic-strided.
```

Mixed precision is also directly supported by this interface.

## Testing LinearAlgebra Frontend

As the `LinearAlgebra` frontend overwrites Julia's Base methods, it's recommended to run tests before using. Currently level-3 tests are implemented within the package and can be invoked by:

```
using Pkg
Pkg.test("BLIS")
# or simply `]test BLIS`
```

The testsuite is method-resolved and comes with a failure workaround. If failure occurs with a message like:

```
[ Info: `gemm` test failed. Consider adding it to ~/.blis_jlbla_blacklist.
BLAS level-3 LinearAlgebra interface: Test Failed at ...
```

One can make the following blacklist to block the very failing method and use the rest:

```bash
echo "gemm" > ~/.blis_jlbla_blacklist
# Trigger JIT reload.
touch $(find ~/.julia/packages/ -name "BLIS.jl" | tail -n 1)
```



## `BLIS.BLASInterface`: Level-1

```@docs
BLIS.LinearAlgebra.BLAS.scal!
BLIS.LinearAlgebra.BLAS.blascopy!
BLIS.LinearAlgebra.BLAS.axpy!
BLIS.LinearAlgebra.BLAS.axpby!
```

## `BLIS.BLASInterface`: Level-2

```@docs
BLIS.LinearAlgebra.BLAS.gemv!
BLIS.LinearAlgebra.BLAS.hemv!
BLIS.LinearAlgebra.BLAS.symv!
BLIS.LinearAlgebra.BLAS.trmv!
BLIS.LinearAlgebra.BLAS.trsv!
BLIS.LinearAlgebra.BLAS.ger!
BLIS.LinearAlgebra.BLAS.her!
BLIS.LinearAlgebra.BLAS.syr!
```

## `BLIS.BLASInterface`: Level-3

```@docs
BLIS.LinearAlgebra.BLAS.gemm!
BLIS.LinearAlgebra.BLAS.hemm!
BLIS.LinearAlgebra.BLAS.symm!
BLIS.LinearAlgebra.BLAS.her2k!
BLIS.LinearAlgebra.BLAS.syr2k!
BLIS.LinearAlgebra.BLAS.herk!
BLIS.LinearAlgebra.BLAS.syrk!
BLIS.LinearAlgebra.BLAS.trmm!
BLIS.LinearAlgebra.BLAS.trsm!
```

