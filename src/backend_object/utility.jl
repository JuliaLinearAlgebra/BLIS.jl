# Object-based wrapper for BLIS utility routines.
#

macro blis_object_api_lvutil_mkspetial(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           a,      xObj)
    end
end

@blis_object_api_lvutil_mkspetial mkherm
@blis_object_api_lvutil_mkspetial mksymm
@blis_object_api_lvutil_mkspetial mktrim

macro blis_object_api_lvutil_norm(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           x,      xObj,
                           norm,   xObj)
    end
end

@blis_object_api_lvutil_norm norm1v
@blis_object_api_lvutil_norm normfv
@blis_object_api_lvutil_norm normiv
@blis_object_api_lvutil_norm norm1m
@blis_object_api_lvutil_norm normfm
@blis_object_api_lvutil_norm normim

macro blis_object_api_lvutil_rand(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           x,      xObj)
    end
end

@blis_object_api_lvutil_rand randv
@blis_object_api_lvutil_rand randnv
@blis_object_api_lvutil_rand randm
@blis_object_api_lvutil_rand randnm

macro blis_object_api_lvutil_print(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           s1,     Ptr{Cchar},
                           x,      xObj,
                           format, Ptr{Cchar},
                           s2,     Ptr{Cchar})
    end
end

@blis_ccall_object(asumv,
                   Cvoid,
                   x,     xObj,
                   asum,  xObj)
@blis_ccall_object(sumsqv,
                   Cvoid,
                   x,     xObj,
                   scale, xObj,
                   sumsq, xObj)

