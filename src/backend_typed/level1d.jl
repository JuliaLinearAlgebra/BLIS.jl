# Wrapper for level-1 BLIS matrix-diagonal routines.
#

# Level1d API of common forms are put together.
#
macro blis_group_level1d_form1(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          diagoffa, BliDoff,
                          diaga,    BliDiag,
                          transa,   BliTrans,
                          m,        BliDim,
                          n,        BliDim,
                          a,        Ptr{xType}, rsa, BliInc, csa, BliInc,
                          b,        Ptr{xType}, rsb, BliInc, csb, BliInc)
    end
end

@blis_group_level1d_form1 addd
@blis_group_level1d_form1 copyd
@blis_group_level1d_form1 subd

macro blis_group_level1d_form2(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          diagoffa, BliDoff,
                          diaga,    BliDiag,
                          transa,   BliTrans,
                          m,        BliDim,
                          n,        BliDim,
                          α,        Ptr{xType},
                          a,        Ptr{xType}, rsa, BliInc, csa, BliInc,
                          b,        Ptr{xType}, rsb, BliInc, csb, BliInc)
    end
end

@blis_group_level1d_form2 axpyd
@blis_group_level1d_form2 scal2d

macro blis_group_level1d_form3(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          conjα,    BliConj,
                          diagoffa, BliDoff,
                          m,        BliDim,
                          n,        BliDim,
                          α,        Ptr{xType},
                          a,        Ptr{xType}, rsa, BliInc, csa, BliInc)
    end
end

@blis_group_level1d_form3 scald
@blis_group_level1d_form3 setd


# Forms with only 1 member are defined directly.
#
@blis_ccall_group(invertd,
                  Cvoid,
                  diagoffa, BliDoff,
                  m,        BliDim,
                  n,        BliDim,
                  a,        Ptr{xType}, rsa, BliInc, csa, BliInc)
@blis_ccall_group(setid,
                  Cvoid,
                  diagoffa, BliDoff,
                  m,        BliDim,
                  n,        BliDim,
                  α,        Ptr{rType},
                  a,        Ptr{xType}, rsa, BliInc, csa, BliInc)
@blis_ccall_group(shiftd, Cvoid,
                  diagoffa, BliDoff,
                  m,        BliDim,
                  n,        BliDim,
                  α,        Ptr{rType},
                  a,        Ptr{xType}, rsa, BliInc, csa, BliInc)
@blis_ccall_group(xpbyd, Cvoid,
                  diagoffa, BliDoff,
                  diaga,    BliDiag,
                  transa,   BliTrans,
                  m,        BliDim,
                  n,        BliDim,
                  a,        Ptr{xType}, rsa, BliInc, csa, BliInc,
                  β,        Ptr{xType},
                  b,        Ptr{xType}, rsb, BliInc, csb, BliInc)

