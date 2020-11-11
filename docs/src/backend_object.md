# BLIS Object Backend

LinearAlgebra-like frontend communicates with the BLIS library through its
[Object-API](https://github.com/flame/blis/blob/master/docs/BLISObjectAPI.md), which is exposed as `BLIS.ObjectBackend` submodule.
Mixed-precision BLAS functions are only available via Object-API methods:

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

## `BLIS.ObjectBackend`: Object Creation

```@docs
BLIS.ObjectBackend.BliObj
BLIS.ObjectBackend.BliObjBase
```

## `BLIS.ObjectBackend`: Level-1 Vector

```@docs
BLIS.ObjectBackend.bli_addv!
BLIS.ObjectBackend.bli_copyv!
BLIS.ObjectBackend.bli_subv!
BLIS.ObjectBackend.bli_swapv!
BLIS.ObjectBackend.bli_axpyv!
BLIS.ObjectBackend.bli_scal2v!
BLIS.ObjectBackend.bli_scalv!
BLIS.ObjectBackend.bli_setv!
BLIS.ObjectBackend.bli_setrv!
BLIS.ObjectBackend.bli_setiv!
BLIS.ObjectBackend.bli_amaxv!
BLIS.ObjectBackend.bli_axpbyv!
BLIS.ObjectBackend.bli_dotv!
BLIS.ObjectBackend.bli_dotxv!
BLIS.ObjectBackend.bli_invertv!
BLIS.ObjectBackend.bli_xpbyv!
```

## `BLIS.ObjectBackend`: Level-1 Diagonal

```@docs
BLIS.ObjectBackend.bli_addd!
BLIS.ObjectBackend.bli_copyd!
BLIS.ObjectBackend.bli_subd!
BLIS.ObjectBackend.bli_axpyd!
BLIS.ObjectBackend.bli_scal2d!
BLIS.ObjectBackend.bli_scald!
BLIS.ObjectBackend.bli_setd!
BLIS.ObjectBackend.bli_setid!
BLIS.ObjectBackend.bli_shiftd!
BLIS.ObjectBackend.bli_invertd!
BLIS.ObjectBackend.bli_xpbyd!
```

## `BLIS.ObjectBackend`: Level-1 Matrix

```@docs
BLIS.ObjectBackend.bli_addm!
BLIS.ObjectBackend.bli_copym!
BLIS.ObjectBackend.bli_subm!
BLIS.ObjectBackend.bli_axpym!
BLIS.ObjectBackend.bli_scal2m!
BLIS.ObjectBackend.bli_scalm!
BLIS.ObjectBackend.bli_setm!
BLIS.ObjectBackend.bli_xpbym!
BLIS.ObjectBackend.bli_xpbym_md!
```

## `BLIS.ObjectBackend`: Level-1 Fused-vector

```@docs
BLIS.ObjectBackend.bli_axpy2v!
BLIS.ObjectBackend.bli_axpyf!
BLIS.ObjectBackend.bli_dotxf!
```

## `BLIS.ObjectBackend`: Level-2

```@docs
BLIS.ObjectBackend.bli_gemv!
BLIS.ObjectBackend.bli_hemv!
BLIS.ObjectBackend.bli_symv!
BLIS.ObjectBackend.bli_ger!
BLIS.ObjectBackend.bli_her2!
BLIS.ObjectBackend.bli_syr2!
BLIS.ObjectBackend.bli_her!
BLIS.ObjectBackend.bli_syr!
BLIS.ObjectBackend.bli_trsv!
BLIS.ObjectBackend.bli_trmv!
```

## `BLIS.ObjectBackend`: Level-3

```@docs
BLIS.ObjectBackend.bli_gemm!
BLIS.ObjectBackend.bli_hemm!
BLIS.ObjectBackend.bli_symm!
BLIS.ObjectBackend.bli_her2k!
BLIS.ObjectBackend.bli_syr2k!
BLIS.ObjectBackend.bli_herk!
BLIS.ObjectBackend.bli_syrk!
BLIS.ObjectBackend.bli_trmm!
BLIS.ObjectBackend.bli_trmm3!
BLIS.ObjectBackend.bli_trsm!
```

## `BLIS.ObjectBackend`: Utility-level

```@docs
BLIS.ObjectBackend.bli_mkherm!
BLIS.ObjectBackend.bli_mksymm!
BLIS.ObjectBackend.bli_mktrim!
BLIS.ObjectBackend.bli_norm1v!
BLIS.ObjectBackend.bli_normfv!
BLIS.ObjectBackend.bli_normiv!
BLIS.ObjectBackend.bli_norm1m!
BLIS.ObjectBackend.bli_normfm!
BLIS.ObjectBackend.bli_normim!
BLIS.ObjectBackend.bli_randv!
BLIS.ObjectBackend.bli_randnv!
BLIS.ObjectBackend.bli_randm!
BLIS.ObjectBackend.bli_randnm!
BLIS.ObjectBackend.bli_asumv!
BLIS.ObjectBackend.bli_sumsqv!
```

