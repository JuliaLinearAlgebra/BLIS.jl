# Object-based wrapper for level-1 BLIS vector routines.
#

# Some methods here might share the same type lists.
# Still, they're made separate to give help message a clear meaning.
#
macro blis_object_api_level1v_form1(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           x, xObj,
                           y, xObj)
    end
end

@blis_object_api_level1v_form1 addv
@blis_object_api_level1v_form1 copyv
@blis_object_api_level1v_form1 subv
@blis_object_api_level1v_form1 swapv

macro blis_object_api_level1v_form2(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           α, xObj,
                           x, xObj,
                           y, xObj)
    end
end

@blis_object_api_level1v_form2 axpyv
@blis_object_api_level1v_form2 scal2v

macro blis_object_api_level1v_form3(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           α, xObj,
                           x, xObj)
    end
end

@blis_object_api_level1v_form3 scalv
@blis_object_api_level1v_form3 setv
@blis_object_api_level1v_form3 setrv
@blis_object_api_level1v_form3 setiv

@blis_ccall_object(amaxv,
                   Cvoid,
                   x,     xObj,
                   index, xObj)
@blis_ccall_object(axpbyv,
                   Cvoid,
                   α, xObj,
                   x, xObj,
                   β, xObj,
                   y, xObj)
@blis_ccall_object(dotv,
                   Cvoid,
                   x, xObj,
                   y, xObj,
                   ρ, xObj)
@blis_ccall_object(dotxv,
                   Cvoid,
                   α, xObj,
                   x, xObj,
                   y, xObj,
                   β, xObj,
                   ρ, xObj)
@blis_ccall_object(invertv,
                   Cvoid,
                   x, xObj)
@blis_ccall_object(xpbyv,
                   Cvoid,
                   x, xObj,
                   β, xObj,
                   y, xObj)

