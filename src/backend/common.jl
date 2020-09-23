# Common backend methods.

cblas_typechar = Dict( Cfloat     => 's',
                       Cdouble    => 'd',
                       ComplexF32 => 'c',
                       ComplexF64 => 'z' )

# Common interface for committing backend calls.
backend_commit(func::Symbol, ret::Type, params::NTuple{NPar, Type}, 
               kwargs...) where {NPar} = ccall((func, blis), ret, params, kwargs...)

# Macros for quick definition of call rules.
# Generate single call
macro blis_ccall_direct(funcname, ret, params)
    return :( Symbol("bli_", $funcname) => CCallInfo(($ret, $params)) )
end

# Generate single call, replace xType with typename and assign letter.
macro blis_ccall_typed(funcname, ret, params, typename)
    return :( begin
        typed_params = ( xType -> $params )($typename)
        typed_ret    = ( xType -> $ret    )($typename)
        Symbol("bli_", cblas_typechar[$typename], $funcname) => CCallInfo((typed_ret, typed_params))
    end )
end

macro blis_ccall_group(funcname, ret, params)
    return :( ( @blis_ccall_typed($funcname, $ret, $params, Cfloat    ),
                @blis_ccall_typed($funcname, $ret, $params, Cdouble   ),
                @blis_ccall_typed($funcname, $ret, $params, ComplexF32),
                @blis_ccall_typed($funcname, $ret, $params, ComplexF64) ) )
end

