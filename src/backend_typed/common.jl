# Common backend methods.
#

# Generate single call, replace [xr]Type with typenames and assign letter.
macro blis_ccall_typed(funcname, ret, ctypename, rtypename, named_param_pairs...)
    # C and Julia function names.
    cfunc  = string("bli_", cblas_typechar[ eval(ctypename) ], funcname)
    jlfunc = Symbol("bli_", cblas_typechar[ eval(ctypename) ], funcname, "!")
    ccfunc = Expr(:string, cfunc)

    # Separate variable names and their types.
    names  = named_param_pairs[1:2:end]
    params = named_param_pairs[2:2:end]

    # For typed parameters and return value, replace [xr]Type with real typenames.
    typed_params = ( iType -> eval(:( ( (xType, rType) -> $iType )($ctypename, $rtypename) )) ).(params)
    typed_return =            eval(:( ( (xType, rType) -> $ret   )($ctypename, $rtypename) ))

    # Convert type specifiers to static symbol.
    typed_param_syms = Expr(:tuple, typed_params...)
    typed_return_sym = Symbol(typed_return)

    # Docstring.
    funcdoc = Expr(:string, """
        Autogenerated BLIS Typed-API function.
        For usage refer the
        [BLIS docs for `$cfunc`.](https://github.com/flame/blis/blob/0.7.0/docs/BLISTypedAPI.md#$funcname)
        """)

    # Define Julia function as ccall.
    return quote
        $(esc(jlfunc))( $(names...) ) = begin 
            ccall(dlsym(libblis, $(esc(ccfunc))), $typed_return_sym, $typed_param_syms, $(names...))
        end

        @doc $funcdoc $jlfunc
    end
end

# Ganerate a famity of calls.
macro blis_ccall_group(funcname, ret, named_param_pairs...)
    return quote 
        @blis_ccall_typed($funcname, $ret, Cfloat    , Cfloat , $(named_param_pairs...))
        @blis_ccall_typed($funcname, $ret, Cdouble   , Cdouble, $(named_param_pairs...))
        @blis_ccall_typed($funcname, $ret, ComplexF32, Cfloat , $(named_param_pairs...))
        @blis_ccall_typed($funcname, $ret, ComplexF64, Cdouble, $(named_param_pairs...))
    end
end

