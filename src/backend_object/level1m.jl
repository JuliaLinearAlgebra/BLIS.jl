# Object-based wrapper for level-1 BLIS matrix routines.
#

macro blis_object_api_level1m_form1(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           a, xObj,
                           b, xObj)
    end
end

@blis_object_api_level1m_form1 addm
@blis_object_api_level1m_form1 copym
@blis_object_api_level1m_form1 subm

macro blis_object_api_level1m_form2(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           α, xObj,
                           a, xObj,
                           b, xObj)
    end
end

@blis_object_api_level1m_form2 axpym
@blis_object_api_level1m_form2 scal2m

macro blis_object_api_level1m_form3(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           α, xObj,
                           a, xObj)
    end
end

@blis_object_api_level1m_form3 scalm
@blis_object_api_level1m_form3 setm

macro blis_object_api_level1m_form4(funcname)
    return quote
        @blis_ccall_object($funcname,
                           Cvoid,
                           a, xObj,
                           β, xObj,
                           b, xObj)
    end
end

@blis_object_api_level1m_form4 xpbym
@blis_object_api_level1m_form4 xpbym_md

