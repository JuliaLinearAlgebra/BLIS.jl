module BLIS

using Libdl
using blis_jll: blis
using LinearAlgebra

global libblis = C_NULL

__init__() = begin
    if length(get(ENV, "BLISDIR", "")) > 0
        # BLIS installation overriden by environmental variables.
        @info "Using custom defined BLIS installation instead of blis_jll."
        global libblis = dlopen(string(get(ENV, "BLISDIR", ""), "/lib/libblis"))
    else
        # Use BinaryBuilder provided BLIS library.
        global libblis = dlopen(blis)
    end
end

# Data types.
module Types
include("types.jl")
end

# Typed-API backend.
module TypedBackend
import ..blis
import ..libblis
using ..Libdl
using ..Types
using ..Types: cblas_typechar
using ..Types: BliCompatibleType
using ..Types: BliDim, BliInc, BliDoff, BliSiz
using ..Types: BliTrans, BliConj, BliUpLo, BliSide, BliDiag, BliInvDiag, BliStruc
include("backend_typed/common.jl")
include("backend_typed/level1v.jl")
include("backend_typed/level1d.jl")
include("backend_typed/level1m.jl")
include("backend_typed/level1f.jl")
include("backend_typed/level2.jl")
include("backend_typed/level3.jl")
include("backend_typed/utility.jl")
end

# Object-API backend.
module ObjectBackend
import ..blis
import ..libblis
using ..Libdl
using ..Types: BliDim, BliInc, BliDoff, BliSiz
using ..Types: BliObjBits, BliAtomic, BliNum
using ..Types: BliSide
include("backend_object/object.jl")
include("backend_object/common.jl")
include("backend_object/level1v.jl")
include("backend_object/level1d.jl")
include("backend_object/level1m.jl")
include("backend_object/level1f.jl")
include("backend_object/level2.jl")
include("backend_object/level3.jl")
include("backend_object/utility.jl")
end

# LinearAlgebra BLAS interface.
module BLASInterface
import ..TypedBackend
import ..LinearAlgebra.BLAS: gemm!, gemm
import ..LinearAlgebra.BLAS: hemm!, hemm
import ..LinearAlgebra.BLAS: symm!, symm
import ..LinearAlgebra.BLAS: trmm!, trmm
import ..LinearAlgebra.BLAS: her2k!, her2k
import ..LinearAlgebra.BLAS: syr2k!, syr2k
using ..Types
using ..Types: cblas_typechar
using ..Types: char_to_trans, char_to_conj, char_to_side, char_to_uplo, char_to_diag
include("interface_linalg/gemm.jl")
include("interface_linalg/hemm.jl")
include("interface_linalg/trmm.jl")
include("interface_linalg/her2k.jl")
end

end

