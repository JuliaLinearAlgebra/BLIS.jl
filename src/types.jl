# BLIS datatypes for uses within Julia and for calling backends.
#

# BLAS-like type characters.
cblas_typechar = Dict( Float32    => 's',
                       Float64    => 'd',
                       ComplexF32 => 'c',
                       ComplexF64 => 'z' )

BliCompatibleType = Union{keys(cblas_typechar)...}

# Int type should be consistent with system word size.
if Sys.WORD_SIZE == 64
    BliInt  = Int64
    BliUInt = UInt64
elseif Sys.WORD_SIZE == 32
    BliInt  = Int32
    BliUInt = UInt32
end

# Integer parameters.
BliDim  = BliInt
BliInc  = BliInt
BliDoff = BliInt
BliSiz  = BliUInt

# Integer params for object API.
BliObjBits = UInt32
BliAtomic  = ComplexF64

# Enumerates, defined as structures.
#--------------------------------------------------
# Enumerate types generated from list_blis_enums.c 
struct BliTrans enum::BliInt end
const BLIS_NO_TRANSPOSE      = BliTrans(0)
const BLIS_TRANSPOSE         = BliTrans(8)
const BLIS_CONJ_NO_TRANSPOSE = BliTrans(16)
const BLIS_CONJ_TRANSPOSE    = BliTrans(24)

struct BliConj enum::BliInt end
const BLIS_NO_CONJUGATE      = BliConj(0)
const BLIS_CONJUGATE         = BliConj(16)

struct BliUpLo enum::BliInt end
const BLIS_ZEROS             = BliUpLo(0)
const BLIS_LOWER             = BliUpLo(192)
const BLIS_UPPER             = BliUpLo(96)
const BLIS_DENSE             = BliUpLo(224)

struct BliSide enum::BliInt end
const BLIS_LEFT              = BliSide(0)
const BLIS_RIGHT             = BliSide(1)

struct BliDiag enum::BliInt end
const BLIS_NONUNIT_DIAG      = BliDiag(0)
const BLIS_UNIT_DIAG         = BliDiag(256)

struct BliInvDiag enum::BliInt end
const BLIS_NO_INVERT_DIAG    = BliInvDiag(0)
const BLIS_INVERT_DIAG       = BliInvDiag(512)

struct BliStruc enum::BliInt end
const BLIS_GENERAL           = BliStruc(0)
const BLIS_HERMITIAN         = BliStruc(134217728)
const BLIS_SYMMETRIC         = BliStruc(268435456)
const BLIS_TRIANGULAR        = BliStruc(402653184)

struct BliNum enum::BliInt end
const BLIS_FLOAT             = BliNum(0)
const BLIS_DOUBLE            = BliNum(2)
const BLIS_SCOMPLEX          = BliNum(1)
const BLIS_DCOMPLEX          = BliNum(3)
const BLIS_INT               = BliNum(4)
const BLIS_CONSTANT          = BliNum(5)
const BLIS_DT_LO             = BliNum(0)
const BLIS_DT_HI             = BliNum(3)

const BLIS_CONJTRANS_BITS    = BliObjBits(24)
const BLIS_DATATYPE_BITS     = BliObjBits(7)
const BLIS_DOMAIN_BIT        = BliObjBits(1)
const BLIS_PRECISION_BIT     = BliObjBits(2)
const BLIS_CONJTRANS_BITS    = BliObjBits(24)
const BLIS_TRANS_BIT         = BliObjBits(8)
const BLIS_CONJ_BIT          = BliObjBits(16)
const BLIS_UPLO_BITS         = BliObjBits(224)
const BLIS_UPPER_BIT         = BliObjBits(32)
const BLIS_DIAG_BIT          = BliObjBits(64)
const BLIS_LOWER_BIT         = BliObjBits(128)
const BLIS_UNIT_DIAG_BIT     = BliObjBits(256)
const BLIS_INVERT_DIAG_BIT   = BliObjBits(512)
const BLIS_TARGET_DT_BITS    = BliObjBits(7168)
const BLIS_TARGET_DOMAIN_BIT = BliObjBits(1024)
const BLIS_TARGET_PREC_BIT   = BliObjBits(2048)
const BLIS_EXEC_DT_BITS      = BliObjBits(57344)
const BLIS_EXEC_DOMAIN_BIT   = BliObjBits(8192)
const BLIS_EXEC_PREC_BIT     = BliObjBits(16384)
const BLIS_STRUC_BITS        = BliObjBits(402653184)

#-------------------------
# Export enumerate types.
export BLIS_NO_TRANSPOSE, BLIS_TRANSPOSE, BLIS_CONJ_NO_TRANSPOSE, BLIS_CONJ_TRANSPOSE
export BLIS_NO_CONJUGATE, BLIS_CONJUGATE
export BLIS_ZEROS, BLIS_LOWER, BLIS_UPPER, BLIS_DENSE
export BLIS_LEFT, BLIS_RIGHT
export BLIS_NONUNIT_DIAG, BLIS_UNIT_DIAG
export BLIS_NO_INVERT_DIAG, BLIS_INVERT_DIAG
export BLIS_GENERAL, BLIS_HERMITIAN, BLIS_SYMMETRIC, BLIS_TRIANGULAR
export BLIS_FLOAT, BLIS_DOUBLE, BLIS_SCOMPLEX, BLIS_DCOMPLEX, BLIS_IN, BLIS_CONSTANT
export BLIS_DT_LO, BLIS_DT_HI

#------------------------------------------
# Lookup table for traditional BLAS chars.
char_to_trans = Dict( 'N' => BLIS_NO_TRANSPOSE,
                      'n' => BLIS_NO_TRANSPOSE,
                      'T' => BLIS_TRANSPOSE,
                      't' => BLIS_TRANSPOSE,
                      'C' => BLIS_CONJ_TRANSPOSE,
                      'c' => BLIS_CONJ_TRANSPOSE )

char_to_conj = Dict( 'N' => BLIS_NO_CONJUGATE,
                     'n' => BLIS_NO_CONJUGATE,
                     'C' => BLIS_CONJUGATE,
                     'c' => BLIS_CONJUGATE )

char_to_side = Dict( 'L' => BLIS_LEFT,
                     'l' => BLIS_LEFT,
                     'R' => BLIS_RIGHT,
                     'r' => BLIS_RIGHT )

char_to_uplo = Dict( 'U' => BLIS_UPPER,
                     'u' => BLIS_UPPER,
                     'L' => BLIS_LOWER,
                     'l' => BLIS_LOWER )

char_to_diag = Dict( 'N' => BLIS_NONUNIT_DIAG,
                     'n' => BLIS_NONUNIT_DIAG,
                     'U' => BLIS_UNIT_DIAG,
                     'u' => BLIS_UNIT_DIAG )

#-----------------------------------------------------------
# Lookup table from Julia datatype to BLIS type indicators.
ctype_to_bli_num = Dict( Float32    => BLIS_FLOAT,
                         Float64    => BLIS_DOUBLE,
                         ComplexF32 => BLIS_SCOMPLEX,
                         ComplexF64 => BLIS_DCOMPLEX,
                         BliInt     => BLIS_INT )

