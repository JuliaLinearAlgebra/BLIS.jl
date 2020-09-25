# Common backend methods.
#

cblas_typechar = Dict( Cfloat     => 's',
                       Cdouble    => 'd',
                       ComplexF32 => 'c',
                       ComplexF64 => 'z' )

# Generate single call, replace [xr]Type with typenames and assign letter.
macro blis_ccall_typed(funcname, ret, params, ctypename, rtypename)
    # C and Julia function names.
    cfunc  = string("bli_", cblas_typechar[ eval(ctypename) ], funcname)
    jlfunc = Symbol("bli_", cblas_typechar[ eval(ctypename) ], funcname, "!")
    ccfunc = Expr(:string, cfunc)

    # Typed parameters and return value.
    typed_params = Expr(:tuple, eval(:( ( (xType, rType) -> $params )($ctypename, $rtypename) ))...)
    typed_ret    = :( $(        eval(:( ( (xType, rType) -> $ret    )($ctypename, $rtypename) ))) )

    # Generate parameter list as Julia argument list.
    jlarg = ( :( $(Symbol(string("var", i))) ) for (i, dtype)=enumerate(eval(typed_params)) )

    # Define Julia function as ccall.
    return quote
             $(esc(jlfunc))( $(jlarg...) ) = begin 
                 ccall(dlsym(libblis, $(esc(ccfunc))), $typed_ret, $typed_params, $(jlarg...))
             end
         end
end

# Ganerate a famity of calls.
macro blis_ccall_group(funcname, ret, params)
    return quote 
        @blis_ccall_typed($funcname, $ret, $params, Cfloat    , Cfloat )
        @blis_ccall_typed($funcname, $ret, $params, Cdouble   , Cdouble)
        @blis_ccall_typed($funcname, $ret, $params, ComplexF32, Cfloat )
        @blis_ccall_typed($funcname, $ret, $params, ComplexF64, Cdouble)
    end
end

