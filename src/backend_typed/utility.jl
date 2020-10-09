# Wrapper for some BLIS utility routines.
#

macro blis_group_lvutil_normv(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          n,        BliDim,
                          x,        Ptr{xType}, incx, BliInc,
                          norm,     Ptr{rType})
    end
end

@blis_group_lvutil_normv norm1v
@blis_group_lvutil_normv normfv
@blis_group_lvutil_normv normiv

macro blis_group_lvutil_normm(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          diagoffa, BliDoff,
                          diaga,    BliDoff,
                          uploa,    BliUpLo,
                          m,        BliDim,
                          n,        BliDim,
                          a,        Ptr{xType}, rsa,  BliInc, csa, BliInc,
                          norm,     Ptr{rType})
    end
end

@blis_group_lvutil_normm norm1m
@blis_group_lvutil_normm normfm
@blis_group_lvutil_normm normim

macro blis_group_lvutil_mkspecial(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          uploa,    BliUpLo,
                          m,        BliDim,
                          a,        Ptr{xType}, rsa,  BliInc, csa, BliInc)
    end
end

@blis_group_lvutil_mkspecial mkherm
@blis_group_lvutil_mkspecial mksymm
@blis_group_lvutil_mkspecial mktrim

macro blis_group_lvutil_randv(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          n,        BliDim,
                          x,        Ptr{xType}, incx, BliInc)
    end
end

@blis_group_lvutil_randv randv
@blis_group_lvutil_randv randnv

macro blis_group_lvutil_randm(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          diagoffa, BliDoff,
                          uploa,    BliUpLo,
                          m,        BliDim,
                          n,        BliDim,
                          a,        Ptr{xType}, rsa,  BliInc, csa, BliInc)
    end
end

@blis_group_lvutil_randm randm
@blis_group_lvutil_randm randnm

@blis_ccall_group(asumv,
                  Cvoid,
                  n,        BliDim,
                  x,        Ptr{xType}, incx, BliInc,
                  asum,     Ptr{rType})
@blis_ccall_group(printv,
                  Cvoid,
                  s1,       Ptr{Cchar},
                  m,        BliDim,
                  x,        Ptr{xType}, incx, BliInc,
                  format,   Ptr{Cchar},
                  s2,       Ptr{Cchar})
@blis_ccall_group(printm,
                  Cvoid,
                  s1,       Ptr{Cchar},
                  m,        BliDim,
                  n,        BliDim,
                  a,        Ptr{xType}, rsa,  BliInc, csa, BliInc,
                  format,   Ptr{Cchar},
                  s2,       Ptr{Cchar})
@blis_ccall_group(sumsqv,
                  Cvoid,
                  n,        BliDim,
                  x,        Ptr{xType}, incx, BliInc,
                  scale,    Ptr{rType},
                  sumsq,    Ptr{rType})

