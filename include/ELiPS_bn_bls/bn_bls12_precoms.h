//
//  bn_bls12_precoms.h
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
 * Interaface for BLS12's pre-computationa and constants
 *
 * @ingroup bls12
 */

#ifndef bn_precoms_h
#define bn_precoms_h

#include <ELiPS_bn_bls/bn_fp6.h>
#include <ELiPS_bn_bls/curve_settings.h>

/*============================================================================*/
/* Constant definitions                                                       */
/*============================================================================*/
#define d12 12
#define d24 24

/*============================================================================*/
/* Enum definitions                                                          */
/*============================================================================*/

/**
 * Enumate the p-th power
 */
enum state{
    f_p1,f_p2,f_p3,f_p4,f_p5,f_p6,f_p7,f_p8,f_p9,f_p10,f_p11,f_p12
};

/*============================================================================*/
/* Global variable definitions                                                */
/*============================================================================*/

/**
 * Basis element in prime field.
 */
extern struct Fp Fp_basis;

/**
 * Basis element in Fp2 extension field.
 */
extern struct Fp2 Fp2_basis;

/**
 * Inverted Fp2 basis element.
 */
extern struct Fp2 Fp2_basis_inv;

/**
 * Basis element in Fp6 extension field.
 */
extern struct Fp6 Fp6_basis;

/**
 * Primitive Cube root of 1.
 */
extern mpz_t epsilon1,epsilon2;

/**
 * Constant values for Frobenius map (FM) to multiplied during FM calculations.
 */
extern  Fp2 d12_frobenius_constant[d12][6];

/**
 * Constant values for skew Frobenius map (FM) to multiplied during skew FM calculations.
 */
extern  Fp2 d12_skew_frobenius_constant[d12][2];

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/
//TODO:
/**
 * Initialize BLS12 curve pre-computations
 * 
 *  @param[in] 	curvetype		    -  send 2 for BLS12 and 1 for BN 
 */
extern void init_precoms(int curvetype);

/**
 * Calculate primitive cubic root of 1 in Fp.
 */
extern void get_epsilon(void);

/**
 * Set basis elements
 */
extern void set_basis(void);

/**
 * Set Frobenius map constants
 */
extern void set_frobenius_constant(void);
#endif /* bn_precoms_h */
