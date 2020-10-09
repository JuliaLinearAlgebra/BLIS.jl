# Wrapper for level-1 BLIS matrix routines.
#

# Level1v API of common forms are put together.
#
macro blis_group_level1m_form1(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid,
                          diagoffa, BliDoff,
                          diaga,    BliDiag,
                          uploa,    BliUpLo,
                          transa,   BliTrans,
                          m,        BliDim,
                          n,        BliDim,
                          a,        Ptr{xType}, rsa, BliInc, csa, BliInc,
                          b,        Ptr{xType}, rsb, BliInc, csb, BliInc)
    end
end

@blis_group_level1m_form1 addm
@blis_group_level1m_form1 copym
@blis_group_level1m_form1 subm

macro blis_group_level1m_form2(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid, 
                          diagoffa, BliDoff,
                          diaga,    BliDiag,
                          uploa,    BliUpLo,
                          transa,   BliTrans,
                          m,        BliDim,
                          n,        BliDim,
                          α,        Ptr{xType},
                          a,        Ptr{xType}, rsa, BliInc, csa, BliInc,
                          b,        Ptr{xType}, rsb, BliInc, csb, BliInc)
    end
end

@blis_group_level1m_form2 axpym
@blis_group_level1m_form2 scal2m

macro blis_group_level1m_form3(funcname)
    return quote
        @blis_ccall_group($funcname,
                          Cvoid, 
                          conjα,    BliConj,
                          diagoffa, BliDoff,
                          diaga,    BliDiag,
                          uploa,    BliUpLo,
                          m,        BliDim,
                          n,        BliDim,
                          α,        Ptr{xType},
                          a,        Ptr{xType}, rsa, BliInc, csa, BliInc)
    end
end

@blis_group_level1m_form3 scalm
@blis_group_level1m_form3 setm

