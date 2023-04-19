module BLIS

using Base
using Libdl
using blis_jll
using LinearAlgebra

global libblis = C_NULL

__init__() = begin
    if length(get(ENV, "BLISDIR", "")) > 0
        # BLIS installation overriden by environmental variables.
        @info "Using custom defined BLIS installation instead of blis_jll."
        global libblis = dlopen(joinpath(get(ENV, "BLISDIR", ""), "lib/libblis"))
    else
        blis_path = blis_jll.blis_path
        # Use BinaryBuilder provided BLIS library.
        @info "blis_jll yields BLIS installation: $blis_path."
        global libblis = dlopen(blis_path)
    end
end

# Data types.
module Types
include("types.jl")
end

# Typed-API backend.
module TypedBackend
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
import ..libblis
using ..Libdl
using ..Base: RefValue, unsafe_convert
using ..Types: ctype_to_bli_num
using ..Types: BliCompatibleType
using ..Types: BliDim, BliInc, BliDoff, BliSiz
using ..Types: BliObjBits, BliAtomic, BliNum
using ..Types: BliTrans, BliConj, BliUpLo, BliSide, BliDiag, BliInvDiag, BliStruc
using ..Types: BLIS_CONJTRANS_BITS, BLIS_DATATYPE_BITS, BLIS_DOMAIN_BIT, BLIS_PRECISION_BIT
using ..Types: BLIS_CONJTRANS_BITS, BLIS_TRANS_BIT, BLIS_CONJ_BIT, BLIS_UPLO_BITS
using ..Types: BLIS_UPPER_BIT, BLIS_DIAG_BIT, BLIS_LOWER_BIT, BLIS_UNIT_DIAG_BIT
using ..Types: BLIS_INVERT_DIAG_BIT, BLIS_STRUC_BITS
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
import ..LinearAlgebra.BLAS: gemm!, hemm!, symm!, trmm!, trsm!, herk!, syrk!, her2k!, syr2k!
import ..LinearAlgebra.BLAS: gemv!, hemv!, symv!, trmv!, trsv!, ger!, her!, syr!
import ..LinearAlgebra.BLAS: axpby!, axpy!, scal!, blascopy!
import ..LinearAlgebra: lapack_size, gemm_wrapper!
import ..LinearAlgebra: Adjoint, StridedVecOrMat
export gemm!, hemm!, symm!, trmm!, trsm!, herk!, syrk!, her2k!, syr2k!
export gemv!, hemv!, symv!, trmv!, trsv!, ger!, her!, syr!
export axpby!, axpy!, scal!, blascopy!
export gemm_wrapper!
using ..TypedBackend
using ..ObjectBackend
using ..Types
using ..Types: cblas_typechar
using ..Types: char_to_trans, char_to_conj, char_to_side, char_to_uplo, char_to_diag
using ..Types: BliTrans, BliConj, BliUpLo, BliSide, BliDiag, BliInvDiag, BliStruc
using ..Types: BliCompatibleType
using ..Types: BLIS_TRANS_BIT
using DelimitedFiles
# Load blacklist.
blacklist = []
if isfile(joinpath(homedir(), ".blis_jlbla_blacklist"))
    blacklist = readdlm(joinpath(homedir(), ".blis_jlbla_blacklist"))
end
include("interface_linalg/level3.jl")
include("interface_linalg/level2.jl")
include("interface_linalg/level1.jl")
end

end

