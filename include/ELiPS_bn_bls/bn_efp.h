//
//  bn_efp.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/1/18.

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
 * Interaface for ellipti cuve operation in BN and  BLS12 curve
 *
 * @ingroup bn12
 */

#ifndef bn_efp_h
#define bn_efp_h

#include <ELiPS_bn_bls/curve_dtypes.h>

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

/**
 * Initialize an EFp type structure 
 *
 * @param[in] P			        -  input P in EFp.
 */
extern void EFp_init(EFp *P);

/**
 * clear of memory an EFp type structure 
 *
 * @param[in] P			        -  input P in EFp.
 */
extern void EFp_clear(EFp *P);

/**
 * Prints an EFp type structure 
 *
 * @param[in] P			        -  input P in EFp.
 * @param[in] str			    -  any string to print.
 */
extern void EFp_printf(EFp *P,char *str);

/**
 * Sets an EFp type structure 
 *
 * @param[out] ANS			    -  output ANS is set as P in EFp.
 * @param[in] P			        -  input P in EFp.
 */
extern void EFp_set(EFp *ANS,EFp *P);

/**
 * Sets an EFp type structure 
 *
 * @param[out] ANS			    -  output ANS is set an unsigned long int.
 * @param[in] UI			    -  input in unsigned long int 
 */
extern void EFp_set_ui(EFp *ANS,unsigned long int UI);

/**
 * Sets an EFp type structure 
 *
 * @param[out] ANS			    -  output ANS is set as EFp.
 * @param[in] UI			    -  input in GMP mpz_t 
 */
extern void EFp_set_mpz(EFp *ANS,mpz_t A);

/**
 * Negate EFp 
 *
 * @param[out] ANS			    -  output pointer.
 * @param[in] UI			    -  input pointer in EFp 
 */
extern void EFp_set_neg(EFp *ANS,EFp *P);

/**
 * Generate rational point in EFp for BN curve
 *
 * @param[in] P			    -  passed input pointer in EFp and obtain genterated point
 */
extern void EFp_rational_point_bn(EFp *P);

/**
 * Generate rational point in EFp for BLS12 curve
 *
 * @param[in] P			    -  passed input pointer in EFp and obtain genterated point
 */
extern void EFp_rational_point_bls12(EFp *P);

/**
 * Elliptiptic curve doubling  ANS = 2[P] in EFp for BN and BLS curve
 *
 * @param[out] ANS			-  output
 * @param[in] P			    -  passed input pointer in EFp 
 */
extern void EFp_ECD(EFp *ANS,EFp *P);

/**
 * Elliptiptic curve addition in EFp ANS = P1+P2 for BN and BLS curve
 *
 * @param[out] ANS			    -  output
 * @param[in] P1			    -  passed input pointer in EFp 
 * @param[in] P1			    -  passed input pointer in EFp 
 */
extern void EFp_ECA(EFp *ANS,EFp *P1,EFp *P2);

/**
 * Elliptiptic Scalar Multiplication (SCM) in EFp ANS = [scalar]P for BN and BLS curve
 *
 * @param[out] ANS			    -  output
 * @param[in] P			        -  passed input pointer in EFp 
 * @param[in] P1			    -  input integer  
 */
extern void EFp_SCM(EFp *ANS,EFp *P,mpz_t scalar);

#endif /* bn_efp_h */
