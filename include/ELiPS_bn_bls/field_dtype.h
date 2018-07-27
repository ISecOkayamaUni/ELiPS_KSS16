//
//  field_dtype.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 1/31/18.

/*
 * ELiPS is an Efficient Library for Pairing-based Systems
 * Copyright (C) 2008-2018 ELiPS Authors. Please refer to the COPYRIGHT file
 * for contact information.
 *
 * This file is part of ELiPS Project. ELiPS is legal property of its
 * developers, whose names are not listed here. 
 *
 * ELiPS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * ELiPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ELiPS. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file
 *
 * Interaface of finite field data types.
 * 
 * @ingroup ff
 */

#ifndef field_dtype_h
#define field_dtype_h

#include <ELiPS_bn_bls/Commont_headers.h>

/*============================================================================*/
/* Struct definitions                                                          */
/*============================================================================*/

/**
 * Fp is the basic type that represent prime field element. Consist of one member element of type mpz_t
 */
typedef struct Fp Fp;
struct Fp{
    mpz_t x0;
};

/**
 * Fp2 is degree 2 extension over Fp. Consist of two Fp element.
 */
typedef struct Fp2 Fp2;
struct Fp2{
    struct Fp x0, x1;
};

/**
 * Fp4 is degree 2 extension over Fp2. Consist of two Fp2 element.
 */
typedef struct Fp4 Fp4;
struct Fp4{
    struct Fp2 x0,x1;
};

/**
 * Fp6 is degree 3 extension over Fp2. Consist of three Fp2 element.
 */
typedef struct Fp6 Fp6;
struct Fp6{
    Fp2 x0,x1,x2;
};

/**
 * Fp8 is degree 2 extension over Fp4. Consist of three Fp4 element.
 */
typedef struct Fp8 Fp8;
struct Fp8{
    struct Fp4 x0,x1;
};

/**
 * Fp12 is degree 2 extension over Fp6. Consist of three Fp2 element.
 */
typedef struct Fp12 Fp12;
struct Fp12{
    Fp6 x0,x1;
};

/**
 * Fp16 is degree 2 extension over Fp8. Consist of three Fp8 element.
 */
typedef struct Fp16 Fp16;
struct Fp16{
    struct Fp8 x0,x1;
};


#endif /* field_dtype_h */
