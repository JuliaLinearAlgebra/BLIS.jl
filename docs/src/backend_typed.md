# BLIS Typed Backend

[Typed-API](https://github.com/flame/blis/blob/master/docs/BLISTypedAPI.md) backend of BLIS is also exposed via the `BLIS.TypedBackend` submodule. The following example demonstrates usage of BLIS's [`bli_dmkherm`](https://github.com/flame/blis/blob/master/docs/BLISTypedAPI.md#mkherm) method:

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



## `BLIS.TypedBackend`: Level-1 Vector

```@docs
BLIS.TypedBackend.bli_saddv!
BLIS.TypedBackend.bli_scopyv!
BLIS.TypedBackend.bli_ssubv!
BLIS.TypedBackend.bli_saxpyv!
BLIS.TypedBackend.bli_sscal2v!
BLIS.TypedBackend.bli_sscalv!
BLIS.TypedBackend.bli_ssetv!
BLIS.TypedBackend.bli_samaxv!
BLIS.TypedBackend.bli_saxpbyv!
BLIS.TypedBackend.bli_sdotv!
BLIS.TypedBackend.bli_sdotxv!
BLIS.TypedBackend.bli_sinvertv!
BLIS.TypedBackend.bli_sswapv!
BLIS.TypedBackend.bli_sxpbyv!

BLIS.TypedBackend.bli_daddv!
BLIS.TypedBackend.bli_dcopyv!
BLIS.TypedBackend.bli_dsubv!
BLIS.TypedBackend.bli_daxpyv!
BLIS.TypedBackend.bli_dscal2v!
BLIS.TypedBackend.bli_dscalv!
BLIS.TypedBackend.bli_dsetv!
BLIS.TypedBackend.bli_damaxv!
BLIS.TypedBackend.bli_daxpbyv!
BLIS.TypedBackend.bli_ddotv!
BLIS.TypedBackend.bli_ddotxv!
BLIS.TypedBackend.bli_dinvertv!
BLIS.TypedBackend.bli_dswapv!
BLIS.TypedBackend.bli_dxpbyv!

BLIS.TypedBackend.bli_caddv!
BLIS.TypedBackend.bli_ccopyv!
BLIS.TypedBackend.bli_csubv!
BLIS.TypedBackend.bli_caxpyv!
BLIS.TypedBackend.bli_cscal2v!
BLIS.TypedBackend.bli_cscalv!
BLIS.TypedBackend.bli_csetv!
BLIS.TypedBackend.bli_camaxv!
BLIS.TypedBackend.bli_caxpbyv!
BLIS.TypedBackend.bli_cdotv!
BLIS.TypedBackend.bli_cdotxv!
BLIS.TypedBackend.bli_cinvertv!
BLIS.TypedBackend.bli_cswapv!
BLIS.TypedBackend.bli_cxpbyv!

BLIS.TypedBackend.bli_zaddv!
BLIS.TypedBackend.bli_zcopyv!
BLIS.TypedBackend.bli_zsubv!
BLIS.TypedBackend.bli_zaxpyv!
BLIS.TypedBackend.bli_zscal2v!
BLIS.TypedBackend.bli_zscalv!
BLIS.TypedBackend.bli_zsetv!
BLIS.TypedBackend.bli_zamaxv!
BLIS.TypedBackend.bli_zaxpbyv!
BLIS.TypedBackend.bli_zdotv!
BLIS.TypedBackend.bli_zdotxv!
BLIS.TypedBackend.bli_zinvertv!
BLIS.TypedBackend.bli_zswapv!
BLIS.TypedBackend.bli_zxpbyv!
```

## `BLIS.TypedBackend`: Level-1 Diagonal

```@docs
BLIS.TypedBackend.bli_saddd!
BLIS.TypedBackend.bli_scopyd!
BLIS.TypedBackend.bli_ssubd!
BLIS.TypedBackend.bli_saxpyd!
BLIS.TypedBackend.bli_sscal2d!
BLIS.TypedBackend.bli_sscald!
BLIS.TypedBackend.bli_ssetd!
BLIS.TypedBackend.bli_sinvertd!
BLIS.TypedBackend.bli_ssetid!
BLIS.TypedBackend.bli_sshiftd!
BLIS.TypedBackend.bli_sxpbyd!

BLIS.TypedBackend.bli_daddd!
BLIS.TypedBackend.bli_dcopyd!
BLIS.TypedBackend.bli_dsubd!
BLIS.TypedBackend.bli_daxpyd!
BLIS.TypedBackend.bli_dscal2d!
BLIS.TypedBackend.bli_dscald!
BLIS.TypedBackend.bli_dsetd!
BLIS.TypedBackend.bli_dinvertd!
BLIS.TypedBackend.bli_dsetid!
BLIS.TypedBackend.bli_dshiftd!
BLIS.TypedBackend.bli_dxpbyd!

BLIS.TypedBackend.bli_caddd!
BLIS.TypedBackend.bli_ccopyd!
BLIS.TypedBackend.bli_csubd!
BLIS.TypedBackend.bli_caxpyd!
BLIS.TypedBackend.bli_cscal2d!
BLIS.TypedBackend.bli_cscald!
BLIS.TypedBackend.bli_csetd!
BLIS.TypedBackend.bli_cinvertd!
BLIS.TypedBackend.bli_csetid!
BLIS.TypedBackend.bli_cshiftd!
BLIS.TypedBackend.bli_cxpbyd!

BLIS.TypedBackend.bli_zaddd!
BLIS.TypedBackend.bli_zcopyd!
BLIS.TypedBackend.bli_zsubd!
BLIS.TypedBackend.bli_zaxpyd!
BLIS.TypedBackend.bli_zscal2d!
BLIS.TypedBackend.bli_zscald!
BLIS.TypedBackend.bli_zsetd!
BLIS.TypedBackend.bli_zinvertd!
BLIS.TypedBackend.bli_zsetid!
BLIS.TypedBackend.bli_zshiftd!
BLIS.TypedBackend.bli_zxpbyd!
```

## `BLIS.TypedBackend`: Level-1 Matrix

```@docs
BLIS.TypedBackend.bli_saddm!
BLIS.TypedBackend.bli_scopym!
BLIS.TypedBackend.bli_ssubm!
BLIS.TypedBackend.bli_saxpym!
BLIS.TypedBackend.bli_sscal2m!
BLIS.TypedBackend.bli_sscalm!
BLIS.TypedBackend.bli_ssetm!

BLIS.TypedBackend.bli_daddm!
BLIS.TypedBackend.bli_dcopym!
BLIS.TypedBackend.bli_dsubm!
BLIS.TypedBackend.bli_daxpym!
BLIS.TypedBackend.bli_dscal2m!
BLIS.TypedBackend.bli_dscalm!
BLIS.TypedBackend.bli_dsetm!

BLIS.TypedBackend.bli_caddm!
BLIS.TypedBackend.bli_ccopym!
BLIS.TypedBackend.bli_csubm!
BLIS.TypedBackend.bli_caxpym!
BLIS.TypedBackend.bli_cscal2m!
BLIS.TypedBackend.bli_cscalm!
BLIS.TypedBackend.bli_csetm!

BLIS.TypedBackend.bli_zaddm!
BLIS.TypedBackend.bli_zcopym!
BLIS.TypedBackend.bli_zsubm!
BLIS.TypedBackend.bli_zaxpym!
BLIS.TypedBackend.bli_zscal2m!
BLIS.TypedBackend.bli_zscalm!
BLIS.TypedBackend.bli_zsetm!
```

## `BLIS.TypedBackend`: Level-1 Fused-vector

```@docs
BLIS.TypedBackend.bli_saxpy2v!
BLIS.TypedBackend.bli_sdotaxpyv!
BLIS.TypedBackend.bli_saxpyf!
BLIS.TypedBackend.bli_sdotxf!
BLIS.TypedBackend.bli_sdotxaxpyf!

BLIS.TypedBackend.bli_daxpy2v!
BLIS.TypedBackend.bli_ddotaxpyv!
BLIS.TypedBackend.bli_daxpyf!
BLIS.TypedBackend.bli_ddotxf!
BLIS.TypedBackend.bli_ddotxaxpyf!

BLIS.TypedBackend.bli_caxpy2v!
BLIS.TypedBackend.bli_cdotaxpyv!
BLIS.TypedBackend.bli_caxpyf!
BLIS.TypedBackend.bli_cdotxf!
BLIS.TypedBackend.bli_cdotxaxpyf!

BLIS.TypedBackend.bli_zaxpy2v!
BLIS.TypedBackend.bli_zdotaxpyv!
BLIS.TypedBackend.bli_zaxpyf!
BLIS.TypedBackend.bli_zdotxf!
BLIS.TypedBackend.bli_zdotxaxpyf!
```

## `BLIS.TypedBackend`: Level-2

```@docs
BLIS.TypedBackend.bli_sher2!
BLIS.TypedBackend.bli_ssyr2!
BLIS.TypedBackend.bli_shemv!
BLIS.TypedBackend.bli_ssymv!
BLIS.TypedBackend.bli_strmv!
BLIS.TypedBackend.bli_strsv!
BLIS.TypedBackend.bli_sgemv!
BLIS.TypedBackend.bli_sger!
BLIS.TypedBackend.bli_sher!
BLIS.TypedBackend.bli_ssyr!

BLIS.TypedBackend.bli_dher2!
BLIS.TypedBackend.bli_dsyr2!
BLIS.TypedBackend.bli_dhemv!
BLIS.TypedBackend.bli_dsymv!
BLIS.TypedBackend.bli_dtrmv!
BLIS.TypedBackend.bli_dtrsv!
BLIS.TypedBackend.bli_dgemv!
BLIS.TypedBackend.bli_dger!
BLIS.TypedBackend.bli_dher!
BLIS.TypedBackend.bli_dsyr!

BLIS.TypedBackend.bli_cher2!
BLIS.TypedBackend.bli_csyr2!
BLIS.TypedBackend.bli_chemv!
BLIS.TypedBackend.bli_csymv!
BLIS.TypedBackend.bli_ctrmv!
BLIS.TypedBackend.bli_ctrsv!
BLIS.TypedBackend.bli_cgemv!
BLIS.TypedBackend.bli_cger!
BLIS.TypedBackend.bli_cher!
BLIS.TypedBackend.bli_csyr!

BLIS.TypedBackend.bli_zher2!
BLIS.TypedBackend.bli_zsyr2!
BLIS.TypedBackend.bli_zhemv!
BLIS.TypedBackend.bli_zsymv!
BLIS.TypedBackend.bli_ztrmv!
BLIS.TypedBackend.bli_ztrsv!
BLIS.TypedBackend.bli_zgemv!
BLIS.TypedBackend.bli_zger!
BLIS.TypedBackend.bli_zher!
BLIS.TypedBackend.bli_zsyr!
```

## `BLIS.TypedBackend`: Level-3

```@docs
BLIS.TypedBackend.bli_shemm!
BLIS.TypedBackend.bli_ssymm!
BLIS.TypedBackend.bli_strmm!
BLIS.TypedBackend.bli_strsm!
BLIS.TypedBackend.bli_sgemm!
BLIS.TypedBackend.bli_sherk!
BLIS.TypedBackend.bli_sher2k!
BLIS.TypedBackend.bli_ssyrk!
BLIS.TypedBackend.bli_ssyr2k!
BLIS.TypedBackend.bli_strmm3!

BLIS.TypedBackend.bli_dhemm!
BLIS.TypedBackend.bli_dsymm!
BLIS.TypedBackend.bli_dtrmm!
BLIS.TypedBackend.bli_dtrsm!
BLIS.TypedBackend.bli_dgemm!
BLIS.TypedBackend.bli_dherk!
BLIS.TypedBackend.bli_dher2k!
BLIS.TypedBackend.bli_dsyrk!
BLIS.TypedBackend.bli_dsyr2k!
BLIS.TypedBackend.bli_dtrmm3!

BLIS.TypedBackend.bli_chemm!
BLIS.TypedBackend.bli_csymm!
BLIS.TypedBackend.bli_ctrmm!
BLIS.TypedBackend.bli_ctrsm!
BLIS.TypedBackend.bli_cgemm!
BLIS.TypedBackend.bli_cherk!
BLIS.TypedBackend.bli_cher2k!
BLIS.TypedBackend.bli_csyrk!
BLIS.TypedBackend.bli_csyr2k!
BLIS.TypedBackend.bli_ctrmm3!

BLIS.TypedBackend.bli_zhemm!
BLIS.TypedBackend.bli_zsymm!
BLIS.TypedBackend.bli_ztrmm!
BLIS.TypedBackend.bli_ztrsm!
BLIS.TypedBackend.bli_zgemm!
BLIS.TypedBackend.bli_zherk!
BLIS.TypedBackend.bli_zher2k!
BLIS.TypedBackend.bli_zsyrk!
BLIS.TypedBackend.bli_zsyr2k!
BLIS.TypedBackend.bli_ztrmm3!
```

## `BLIS.TypedBackend`: Utility-level

```@docs
BLIS.TypedBackend.bli_snorm1v!
BLIS.TypedBackend.bli_snormfv!
BLIS.TypedBackend.bli_snormiv!
BLIS.TypedBackend.bli_snorm1m!
BLIS.TypedBackend.bli_snormfm!
BLIS.TypedBackend.bli_snormim!
BLIS.TypedBackend.bli_smkherm!
BLIS.TypedBackend.bli_smksymm!
BLIS.TypedBackend.bli_smktrim!
BLIS.TypedBackend.bli_srandv!
BLIS.TypedBackend.bli_srandnv!
BLIS.TypedBackend.bli_srandm!
BLIS.TypedBackend.bli_srandnm!
BLIS.TypedBackend.bli_sasumv!
BLIS.TypedBackend.bli_sprintv!
BLIS.TypedBackend.bli_sprintm!
BLIS.TypedBackend.bli_ssumsqv!

BLIS.TypedBackend.bli_dnorm1v!
BLIS.TypedBackend.bli_dnormfv!
BLIS.TypedBackend.bli_dnormiv!
BLIS.TypedBackend.bli_dnorm1m!
BLIS.TypedBackend.bli_dnormfm!
BLIS.TypedBackend.bli_dnormim!
BLIS.TypedBackend.bli_dmkherm!
BLIS.TypedBackend.bli_dmksymm!
BLIS.TypedBackend.bli_dmktrim!
BLIS.TypedBackend.bli_drandv!
BLIS.TypedBackend.bli_drandnv!
BLIS.TypedBackend.bli_drandm!
BLIS.TypedBackend.bli_drandnm!
BLIS.TypedBackend.bli_dasumv!
BLIS.TypedBackend.bli_dprintv!
BLIS.TypedBackend.bli_dprintm!
BLIS.TypedBackend.bli_dsumsqv!

BLIS.TypedBackend.bli_cnorm1v!
BLIS.TypedBackend.bli_cnormfv!
BLIS.TypedBackend.bli_cnormiv!
BLIS.TypedBackend.bli_cnorm1m!
BLIS.TypedBackend.bli_cnormfm!
BLIS.TypedBackend.bli_cnormim!
BLIS.TypedBackend.bli_cmkherm!
BLIS.TypedBackend.bli_cmksymm!
BLIS.TypedBackend.bli_cmktrim!
BLIS.TypedBackend.bli_crandv!
BLIS.TypedBackend.bli_crandnv!
BLIS.TypedBackend.bli_crandm!
BLIS.TypedBackend.bli_crandnm!
BLIS.TypedBackend.bli_casumv!
BLIS.TypedBackend.bli_cprintv!
BLIS.TypedBackend.bli_cprintm!
BLIS.TypedBackend.bli_csumsqv!

BLIS.TypedBackend.bli_znorm1v!
BLIS.TypedBackend.bli_znormfv!
BLIS.TypedBackend.bli_znormiv!
BLIS.TypedBackend.bli_znorm1m!
BLIS.TypedBackend.bli_znormfm!
BLIS.TypedBackend.bli_znormim!
BLIS.TypedBackend.bli_zmkherm!
BLIS.TypedBackend.bli_zmksymm!
BLIS.TypedBackend.bli_zmktrim!
BLIS.TypedBackend.bli_zrandv!
BLIS.TypedBackend.bli_zrandnv!
BLIS.TypedBackend.bli_zrandm!
BLIS.TypedBackend.bli_zrandnm!
BLIS.TypedBackend.bli_zasumv!
BLIS.TypedBackend.bli_zprintv!
BLIS.TypedBackend.bli_zprintm!
BLIS.TypedBackend.bli_zsumsqv!
```

