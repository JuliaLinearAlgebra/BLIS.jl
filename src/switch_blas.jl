# LBT-forward BLAS calls
#

global fallback_lib = []

lbt_enable_blis(; clear=false) = begin
    if length(fallback_lib) == 0
        libs = LinearAlgebra.BLAS.lbt_get_config().loaded_libs
        if length(libs) > 0
            # Save currently loaded library path for resetting.
            fallback_lib = [fallback_lib..., ]
        end

        LinearAlgebra.BLAS.lbt_forward(libblis_path; clear=clear)
    else
        @warn "BLIS.lbt_enable_blis: BLIS already loaded. Not doing anything."
    end
end

lbt_disable_blis() = begin
    if length(fallback_lib) > 0
        LinearAlgebra.BLAS.lbt_forward(fallback_lib[1]; clear=true)
    else
        @warn "BLIS.lbt_enable_blis: Fallback library not found. Not doing anything."
    end
end
