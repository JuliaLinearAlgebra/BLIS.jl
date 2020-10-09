# Wrapper for level-2 BLIS routines.
#

# Level2 API of common forms are put together.
#
macro blis_group_level2_her2(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          uploa,  BliUpLo,
                          conjx,  BliConj,
                          conjy,  BliConj,
                          m,      BliDim,
                          α,      Ptr{xType},
                          x,      Ptr{xType}, incx, BliInc,
                          y,      Ptr{xType}, incy, BliInc,
                          a,      Ptr{xType}, rsa,  BliInc, csa, BliInc)
    end
end

@blis_group_level2_her2 her2
@blis_group_level2_her2 syr2

macro blis_group_level2_hemv(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          uploa,  BliUpLo,
                          conja,  BliConj,
                          conjx,  BliConj,
                          m,      BliDim,
                          α,      Ptr{xType},
                          a,      Ptr{xType}, rsa,  BliInc, csa, BliInc,
                          x,      Ptr{xType}, incx, BliInc,
                          β,      Ptr{xType},
                          y,      Ptr{xType}, incy, BliInc)
    end
end

@blis_group_level2_hemv hemv
@blis_group_level2_hemv symv

macro blis_group_level2_tr(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          uploa,  BliUpLo,
                          transa, BliTrans,
                          diaga,  BliDiag,
                          m,      BliDim,
                          α,      Ptr{xType},
                          a,      Ptr{xType}, rsa,  BliInc, csa, BliInc,
                          x,      Ptr{xType}, incx, BliInc)
    end
end

@blis_group_level2_tr trmv
@blis_group_level2_tr trsv


# Forms with only 1 member are defined directly.
#
@blis_ccall_group(gemv,
                  Cvoid,
                  transa, BliTrans,
                  conjx,  BliConj,
                  m,      BliDim,
                  n,      BliDim,
                  α,      Ptr{xType},
                  a,      Ptr{xType}, rsa,  BliInc, csa, BliInc,
                  x,      Ptr{xType}, incx, BliInc,
                  β,      Ptr{xType},
                  y,      Ptr{xType}, incy, BliInc)
@blis_ccall_group(ger,
                  Cvoid,
                  conjx,  BliConj,
                  conjy,  BliConj,
                  m,      BliDim,
                  n,      BliDim,
                  α,      Ptr{xType},
                  x,      Ptr{xType}, incx, BliInc,
                  y,      Ptr{xType}, incy, BliInc,
                  a,      Ptr{xType}, rsa,  BliInc, csa, BliInc)
@blis_ccall_group(her,
                  Cvoid,
                  uploa,  BliUpLo,
                  conjx,  BliConj,
                  m,      BliDim,
                  α,      Ptr{rType},
                  x,      Ptr{xType}, incx, BliInc,
                  a,      Ptr{xType}, rsa,  BliInc, csa, BliInc)
@blis_ccall_group(syr,
                  Cvoid,
                  uploa,  BliUpLo,
                  conjx,  BliConj,
                  m,      BliDim,
                  α,      Ptr{xType},
                  x,      Ptr{xType}, incx, BliInc,
                  a,      Ptr{xType}, rsa,  BliInc, csa, BliInc)

