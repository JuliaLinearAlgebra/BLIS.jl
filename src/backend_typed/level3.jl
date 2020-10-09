# Wrapper for level-3 BLIS routines.
#

# Level3 API of common forms are put together.
#
macro blis_group_level3_hemm(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          sidea,   BliSide,
                          uploa,   BliUpLo,
                          conja,   BliConj,
                          transb,  BliTrans,
                          m,       BliDim,
                          n,       BliDim,
                          α,       Ptr{xType},
                          a,       Ptr{xType}, rsa, BliInc, csa, BliInc,
                          b,       Ptr{xType}, rsb, BliInc, csb, BliInc,
                          β,       Ptr{xType},
                          c,       Ptr{xType}, rsc, BliInc, csc, BliInc)
    end
end

@blis_group_level3_hemm hemm
@blis_group_level3_hemm symm

macro blis_group_level3_tr(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          sidea,   BliSide,
                          uploa,   BliUpLo,
                          transa,  BliTrans,
                          diaga,   BliDiag,
                          m,       BliDim,
                          n,       BliDim,
                          α,       Ptr{xType},
                          a,       Ptr{xType}, rsa, BliInc, csa, BliInc,
                          b,       Ptr{xType}, rsb, BliInc, csb, BliInc)
    end
end

@blis_group_level3_tr trmm
@blis_group_level3_tr trsm


# Forms with only 1 member are defined directly.
#
@blis_ccall_group(gemm,
                  Cvoid,
                  transa,  BliTrans,
                  transb,  BliTrans,
                  m,       BliDim,
                  n,       BliDim,
                  k,       BliDim,
                  α,       Ptr{xType},
                  a,       Ptr{xType}, rsa, BliInc, csa, BliInc,
                  b,       Ptr{xType}, rsb, BliInc, csb, BliInc,
                  β,       Ptr{xType},
                  c,       Ptr{xType}, rsc, BliInc, csc, BliInc)
@blis_ccall_group(herk,
                  Cvoid,
                  uploc,   BliUpLo,
                  transa,  BliTrans,
                  m,       BliDim,
                  k,       BliDim,
                  α,       Ptr{rType},
                  a,       Ptr{xType}, rsa, BliInc, csa, BliInc,
                  β,       Ptr{rType},
                  c,       Ptr{xType}, rsc, BliInc, csc, BliInc)
@blis_ccall_group(her2k,
                  Cvoid,
                  uploc,   BliUpLo,
                  transab, BliTrans,
                  m,       BliDim,
                  k,       BliDim,
                  α,       Ptr{xType},
                  a,       Ptr{xType}, rsa, BliInc, csa, BliInc,
                  b,       Ptr{xType}, rsb, BliInc, csb, BliInc,
                  β,       Ptr{rType},
                  c,       Ptr{xType}, rsc, BliInc, csc, BliInc)
@blis_ccall_group(syrk,
                  Cvoid,
                  uploc,   BliUpLo,
                  transa,  BliTrans,
                  m,       BliDim,
                  k,       BliDim,
                  α,       Ptr{xType},
                  a,       Ptr{xType}, rsa, BliInc, csa, BliInc,
                  β,       Ptr{xType},
                  c,       Ptr{xType}, rsc, BliInc, csc, BliInc)
@blis_ccall_group(syr2k,
                  Cvoid,
                  uploc,   BliUpLo,
                  transab, BliTrans,
                  m,       BliDim,
                  k,       BliDim,
                  α,       Ptr{xType},
                  a,       Ptr{xType}, rsa, BliInc, csa, BliInc,
                  b,       Ptr{xType}, rsb, BliInc, csb, BliInc,
                  β,       Ptr{xType},
                  c,       Ptr{xType}, rsc, BliInc, csc, BliInc)
@blis_ccall_group(trmm3,
                  Cvoid,
                  sidea,   BliSide,
                  uploa,   BliUpLo,
                  transa,  BliTrans,
                  diaga,   BliDiag,
                  transb,  BliTrans,
                  m,       BliDim,
                  n,       BliDim,
                  α,       Ptr{xType},
                  a,       Ptr{xType}, rsa, BliInc, csa, BliInc,
                  b,       Ptr{xType}, rsb, BliInc, csb, BliInc,
                  β,       Ptr{xType},
                  c,       Ptr{xType}, rsc, BliInc, csc, BliInc)

