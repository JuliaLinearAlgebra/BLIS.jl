# Macro for defining object-based backend methods.
#

# Generate single call.
macro blis_ccall_object(funcname, ret, named_param_pairs...)
    # C and Julia function names.
    cfunc  = string("bli_", funcname)
    jlfunc = Symbol("bli_", funcname, "!")
    ccfunc = Expr(:string, cfunc)

    # Separate variable names and their types.
    names  = named_param_pairs[1:2:end]
    params = named_param_pairs[2:2:end]

    # For types, replace xObj with appropriate C and Jl call types.
    # ccparams = ( iType -> eval(:( ( xObj -> $iType )(Ptr{BliObjBase}) )) ).(params)
    # jlparams = ( iType -> eval(:( ( xObj -> $iType )(    BliObj     ) )) ).(params)
    # None-Base types cannot be "verified" with the above `eval`
    #  (will yield "BLIS.ObjeckBackend.BliObj" as variable name hence not found).
    # Do this string (symbol) replacing.
    ccparams = (iSym -> if iSym == :xObj
                :( Ptr{BliObjBase} ) else iSym end).(params)
    jlparams = (iSym -> if iSym == :xObj
                :BliObj else iSym end).(params)

    # Convert C typenames into a static symbol of tuple.
    ccparam_syms = Expr(:tuple, ccparams...)

    # Convert Jl typenames into typed-asserted parameter definition.
    jldefs = ( (iArg, iType) -> :( $iArg::$iType ) ).(names, jlparams)

    # Actual arguments with BliObj type should be
    #  converted to BliObjBase when passed to ccall.
    jlargs = ( (iArg, iType) -> begin 
                   if iType == :BliObj
                       :( Ref($iArg.obj) )
                   else
                       iArg
                   end 
               end ).(names, jlparams)

    # Docstring.
    funcdoc = Expr(:string, """
        Autogenerated BLIS Object-API function.
        For usage refer the
        [BLIS docs for `$cfunc`.](https://github.com/flame/blis/blob/0.7.0/docs/BLISObjectAPI.md#$funcname)
        """)

    return quote
        $(esc(jlfunc))( $(jldefs...) ) = begin 
            ccall(dlsym(libblis, $(esc(ccfunc))), $ret, $ccparam_syms, $(jlargs...))
        end

        @doc $funcdoc $jlfunc
    end
end
