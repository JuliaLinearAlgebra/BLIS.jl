# Object-based wrapper for level-3 BLIS routines.
#

macro blis_object_api_level3_gemm(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           α,     xObj,
                           a,     xObj,
                           b,     xObj,
                           β,     xObj,
                           c,     xObj)
    end
end

@blis_object_api_level3_gemm gemm
@blis_object_api_level3_gemm gemmt
@blis_object_api_level3_gemm her2k
@blis_object_api_level3_gemm syr2k

macro blis_object_api_level3_hemm(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           sidea, BliSide,
                           α,     xObj,
                           a,     xObj,
                           b,     xObj,
                           β,     xObj,
                           c,     xObj)
    end
end

@blis_object_api_level3_hemm hemm
@blis_object_api_level3_hemm symm
@blis_object_api_level3_hemm trmm3

macro blis_object_api_level3_herk(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           α,     xObj,
                           a,     xObj,
                           β,     xObj,
                           c,     xObj)
    end
end

@blis_object_api_level3_herk herk
@blis_object_api_level3_herk syrk

macro blis_object_api_level3_trmm(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           sidea, BliSide,
                           α,     xObj,
                           a,     xObj,
                           b,     xObj)
    end
end

@blis_object_api_level3_trmm trmm
@blis_object_api_level3_trmm trsm

