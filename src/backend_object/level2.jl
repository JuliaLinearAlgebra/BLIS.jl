# Object-based wrapper for level-2 BLIS routines.
#

macro blis_object_api_level2_mv(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           α, xObj,
                           a, xObj,
                           x, xObj,
                           β, xObj,
                           y, xObj)
    end
end

@blis_object_api_level2_mv gemv
@blis_object_api_level2_mv hemv
@blis_object_api_level2_mv symv

macro blis_object_api_level2_r2(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           α, xObj,
                           x, xObj,
                           y, xObj,
                           a, xObj)
    end
end

@blis_object_api_level2_r2 ger
@blis_object_api_level2_r2 her2
@blis_object_api_level2_r2 syr2

macro blis_object_api_level2_r(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           α, xObj,
                           x, xObj,
                           a, xObj)
    end
end

@blis_object_api_level2_r her
@blis_object_api_level2_r syr

macro blis_object_api_level2_tr(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           α, xObj,
                           a, xObj,
                           x, xObj)
    end
end

@blis_object_api_level2_tr trmv
@blis_object_api_level2_tr trsv

