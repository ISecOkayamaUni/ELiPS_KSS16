//
//  bn_efp12.h
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
 * Interaface for elliptic cuve operation in BN and  BLS12 curve in Fp2 extension field
 *
 * @ingroup bn12
 */

#ifndef bn_efp12_h
#define bn_efp12_h

#include <ELiPS_bn_bls/bn_efp6.h>
/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

/**
 * Initialize an EFp12 type structure 
 *
 * @param[in] P			        -  input P in EFp12.
 */
extern void EFp12_init(EFp12 *P);

/**
 * Clear memory of an EFp12 type structure 
 *
 * @param[in] P			        -  input P in EFp12.
 */
extern void EFp12_clear(EFp12 *P);

/**
 * Prints an EFp12 rational point
 *
 * @param[in] P			        -  input P as pointer in EFp12.
 * @param[in] str		        -  input str as string.
 */
extern void EFp12_printf(EFp12 *P,char *str);

/**
 * Sets EFp12 point an EFp12 type structure 
 *
 * @param[in] ANS			    -   output.
 * @param[in] P			        -   input P in EFp12.
 */
extern void EFp12_set(EFp12 *ANS,EFp12 *P);

/**
 * Sets unsigned int an EFp12 type structure 
 *
 * @param[in] ANS			    -   output.
 * @param[in] u			        -   input u in us integer.
 */
extern void EFp12_set_ui(EFp12 *ANS,unsigned long int UI);

/**
 * Sets mpz_t in an EFp12 type structure 
 *
 * @param[in] ANS			    -   output.
 * @param[in] str			    -   input A in mpz_t int.
 */
extern void EFp12_set_mpz(EFp12 *ANS,mpz_t A);

/**
 * Negate an EFp12 type structure 
 *
 * @param[in] ANS			    -   output.
 * @param[in] P			        -   input P in EFp12.
 */
extern void EFp12_set_neg(EFp12 *ANS,EFp12 *P);

/**
 * Generate an EFp12 rational point for BN curve
 *
 * @param[in] P			        -   input rational point P and.
 */
extern void EFp12_rational_point_bn(EFp12 *P);

/**
 * Generate an EFp12 rational point BLS12
 *
 * @param[in] P			        -   input rational point P and.
 */
extern void EFp12_rational_point_bls12(EFp12 *P);

/**
 * ECD in EFp12 rational point 
 *
 * @param[in] ANS			    -   output.
 * @param[in] P			        -   input P in EFp12.
 */
extern void EFp12_ECD(EFp12 *ANS,EFp12 *P);

/**
 * ECA in EFp12 rational point 
 *
 * @param[in] ANS			        -   output.
 * @param[in] P1			        -   input P1 in EFp12.
 * @param[in] P2			        -   input P2 in EFp12.
 */
extern void EFp12_ECA(EFp12 *ANS,EFp12 *P1,EFp12 *P2);

/**
 * SCM in EFp12 rational point 
 *
 * @param[in] ANS			        -   output.
 * @param[in] P1			        -   input P1 in EFp12.
 * @param[in] scalar		        -   input scalar in mpz_t integer.
 */
extern void EFp12_SCM(EFp12 *ANS,EFp12 *P,mpz_t scalar);

#endif /* bn_efp12_h */
