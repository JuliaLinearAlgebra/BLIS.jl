# Object-based wrapper for level-1 BLIS matrix-diagonal routines.
#

macro blis_object_api_level1d_form1(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           a, xObj,
                           b, xObj)
    end
end

@blis_object_api_level1d_form1 addd
@blis_object_api_level1d_form1 copyd
@blis_object_api_level1d_form1 subd

macro blis_object_api_level1d_form2(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           α, xObj,
                           a, xObj,
                           b, xObj)
    end
end

@blis_object_api_level1d_form2 axpyd
@blis_object_api_level1d_form2 scal2d

macro blis_object_api_level1d_form3(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           α, xObj,
                           a, xObj)
    end
end

@blis_object_api_level1d_form3 scald
@blis_object_api_level1d_form3 setd
@blis_object_api_level1d_form3 setid
@blis_object_api_level1d_form3 shiftd

@blis_ccall_object(invertd,
                   Cvoid,
                   a, xObj)
@blis_ccall_object(xpbyd,
                   Cvoid,
                   a, xObj,
                   β, xObj,
                   b, xObj)

