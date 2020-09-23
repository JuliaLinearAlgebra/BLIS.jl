# Wrapper for level-1 BLIS routines.
#

level1_interface = begin 
    Dict{Symbol,
         CCallInfo}(@blis_ccall_group("addv", Cvoid, (Cuint, Clong, Ptr{xType}, Clong, Ptr{xType}, Clong))..., 
                    @blis_ccall_group("amaxv", Cvoid, (Clong, Ptr{xType}, Clong, Ptr{Clong}))...,
                    @blis_ccall_group("axpyv", Cvoid, (Cuint, Clong, Ptr{xType}, Ptr{xType}, Clong, Ptr{xType}, Clong))...
                   )
end

level1_commit(func::Symbol, kwargs...) = begin 
    ret, params = level1_interface[func]
    backend_commit(func, ret, params, kwargs...)
end

