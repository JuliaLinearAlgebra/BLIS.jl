/**
 * This file prints out enum members of BLIS.
 * It will be used to (manually) write enumerate sections of types.jl
 */
#include "blis.h"
#include <stdio.h>


void enum_as_jl_construct_member(unsigned long val, char *varname, char *typejl)
{
    printf("const %s = %s(%lu)\n", varname, typejl, val);
}

void enum_as_jl_construct_define(char *varname)
{
    printf("struct %s enum::BliInt end\n", varname);
}

int main(void)
{
    printf("#--------------------------------------------------\n");
    printf("# Enumerate types generated from list_blis_enums.c \n");

    // trans_t
    enum_as_jl_construct_define("BliTrans");
    enum_as_jl_construct_member(BLIS_NO_TRANSPOSE     , "BLIS_NO_TRANSPOSE     ", "BliTrans");
    enum_as_jl_construct_member(BLIS_TRANSPOSE        , "BLIS_TRANSPOSE        ", "BliTrans");
    enum_as_jl_construct_member(BLIS_CONJ_NO_TRANSPOSE, "BLIS_CONJ_NO_TRANSPOSE", "BliTrans");
    enum_as_jl_construct_member(BLIS_CONJ_TRANSPOSE   , "BLIS_CONJ_TRANSPOSE   ", "BliTrans");
    putchar('\n');

    // conj_t
    enum_as_jl_construct_define("BliConj");
    enum_as_jl_construct_member(BLIS_NO_CONJUGATE     , "BLIS_NO_CONJUGATE     ", "BliConj");
    enum_as_jl_construct_member(BLIS_CONJUGATE        , "BLIS_CONJUGATE        ", "BliConj");
    putchar('\n');

    // uplo_t
    enum_as_jl_construct_define("BliUpLo");
    enum_as_jl_construct_member(BLIS_ZEROS            , "BLIS_ZEROS            ", "BliUpLo");
    enum_as_jl_construct_member(BLIS_LOWER            , "BLIS_LOWER            ", "BliUpLo");
    enum_as_jl_construct_member(BLIS_UPPER            , "BLIS_UPPER            ", "BliUpLo");
    enum_as_jl_construct_member(BLIS_DENSE            , "BLIS_DENSE            ", "BliUpLo");
    putchar('\n');

    // side_t
    enum_as_jl_construct_define("BliSide");
    enum_as_jl_construct_member(BLIS_LEFT             , "BLIS_LEFT             ", "BliSide");
    enum_as_jl_construct_member(BLIS_RIGHT            , "BLIS_RIGHT            ", "BliSide");
    putchar('\n');

    // diag_t
    enum_as_jl_construct_define("BliDiag");
    enum_as_jl_construct_member(BLIS_NONUNIT_DIAG     , "BLIS_NONUNIT_DIAG     ", "BliDiag");
    enum_as_jl_construct_member(BLIS_UNIT_DIAG        , "BLIS_UNIT_DIAG        ", "BliDiag");
    putchar('\n');

    // invdiag_t
    enum_as_jl_construct_define("BliInvDiag");
    enum_as_jl_construct_member(BLIS_NO_INVERT_DIAG   , "BLIS_NO_INVERT_DIAG   ", "BliInvDiag");
    enum_as_jl_construct_member(BLIS_INVERT_DIAG      , "BLIS_INVERT_DIAG      ", "BliInvDiag");
    putchar('\n');
    
    // struc_t
    enum_as_jl_construct_define("BliStruc");
    enum_as_jl_construct_member(BLIS_GENERAL          , "BLIS_GENERAL          ", "BliStruc");
    enum_as_jl_construct_member(BLIS_HERMITIAN        , "BLIS_HERMITIAN        ", "BliStruc");
    enum_as_jl_construct_member(BLIS_SYMMETRIC        , "BLIS_SYMMETRIC        ", "BliStruc");
    enum_as_jl_construct_member(BLIS_TRIANGULAR       , "BLIS_TRIANGULAR       ", "BliStruc");
    putchar('\n');

    return 0;
}

