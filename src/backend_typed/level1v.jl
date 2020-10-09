# Wrapper for level-1 BLIS vector routines.
#

# Level1v API of common forms are put together.
#
macro blis_group_level1v_form1(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          conjx, BliConj,
                          n,     BliDim,
                          x,     Ptr{xType}, incx, BliInc,
                          y,     Ptr{xType}, incy, BliInc)
    end
end

@blis_group_level1v_form1 addv
@blis_group_level1v_form1 copyv
@blis_group_level1v_form1 subv

macro blis_group_level1v_form2(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          conjx, BliConj,
                          n,     BliDim,
                          α,     Ptr{xType},
                          x,     Ptr{xType}, incx, BliInc,
                          y,     Ptr{xType}, incy, BliInc)
    end
end

@blis_group_level1v_form2 axpyv
@blis_group_level1v_form2 scal2v

macro blis_group_level1v_form3(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          conjα, BliConj,
                          n,     BliDim,
                          α,     Ptr{xType},
                          x,     Ptr{xType}, incx, BliInc)
    end
end

@blis_group_level1v_form3 scalv
@blis_group_level1v_form3 setv


# Forms with only 1 member are defined directly.
#
@blis_ccall_group(amaxv,
                  Cvoid, 
                  n,     BliDim,
                  x,     Ptr{xType}, incx, BliInc,
                  index, Ptr{Clong})
@blis_ccall_group(axpbyv,
                  Cvoid, 
                  conjx, BliConj,
                  n,     BliDim,
                  α,     Ptr{xType},
                  x,     Ptr{xType}, incx, BliInc,
                  beta,  Ptr{xType},
                  y,     Ptr{xType}, incy, BliInc)
@blis_ccall_group(dotv,
                  Cvoid, 
                  conjx, BliConj,
                  conjy, BliConj,
                  n,     BliDim,
                  x,     Ptr{xType}, incx, BliInc,
                  y,     Ptr{xType}, incy, BliInc,
                  ρ,     Ptr{xType})
@blis_ccall_group(dotxv,
                  Cvoid, 
                  conjx, BliConj,
                  conjy, BliConj,
                  n,     BliDim,
                  α,     Ptr{xType},
                  x,     Ptr{xType}, incx, BliInc,
                  y,     Ptr{xType}, incy, BliInc,
                  β,     Ptr{xType},
                  ρ,     Ptr{xType})
@blis_ccall_group(invertv,
                  Cvoid,
                  n,     BliDim,
                  x,     Ptr{xType}, incx, BliInc)
@blis_ccall_group(swapv,
                  Cvoid, 
                  n,     BliDim,
                  x,     Ptr{xType}, incx, BliInc,
                  y,     Ptr{xType}, incy, BliInc)
@blis_ccall_group(xpbyv,
                  Cvoid, 
                  conjx, BliConj,
                  n,     BliDim,
                  x,     Ptr{xType}, incx, BliInc,
                  β,     Ptr{xType},
                  y,     Ptr{xType}, incy, BliInc)

