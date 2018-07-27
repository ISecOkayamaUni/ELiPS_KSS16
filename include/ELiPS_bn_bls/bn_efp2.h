//
//  bn_efp2.h
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
#ifndef bn_efp2_h
#define bn_efp2_h

#include <ELiPS_bn_bls/bn_efp.h>

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

/**
 * Initialize an EFp type structure 
 *
 * @param[in] P			        -  input P in EFp2.
 */

extern void EFp2_init(EFp2 *P);

/**
 * Clears an EFp2 type structure 
 *
 * @param[in] P			        -  input P in EFp2.
 */
extern void EFp2_clear(EFp2 *P);

/**
 * Prints an EFp2 type structure 
 *
 * @param[in] P			        -  input P in EFp2.
 * @param[in] str			    -  any string to print.
 */
extern void EFp2_printf(EFp2 *P,char *str);
/**
 * Sets EFp2 point an EFp2 type structure 
 *
 * @param[in] ANS			    -   output.
 * @param[in] P			        -   input P in EFp2.
 */
extern void EFp2_set(EFp2 *ANS,EFp2 *P);

/**
 * Sets unsigned int an EFp2 type structure 
 *
 * @param[in] ANS			    -   output.
 * @param[in] u			        -   input u in us integer.
 */
extern void EFp2_set_ui(EFp2 *ANS,unsigned long int u);

/**
 * Sets mpz_t in an EFp2 type structure 
 *
 * @param[in] ANS			    -   output.
 * @param[in] str			    -   input A in mpz_t int.
 */
extern void EFp2_set_mpz(EFp2 *ANS,mpz_t A);

/**
 * Negate an EFp2 type structure 
 *
 * @param[in] ANS			    -   output.
 * @param[in] P			        -   input P in EFp2.
 */
extern void EFp2_set_neg(EFp2 *ANS,EFp2 *P);

/**
 * Generate an EFp2 rational point
 *
 * @param[in] P			        -   input rational point P and.
 */
extern void EFp2_rational_point(EFp2 *P);

/**
 * ECD in EFp2 rational point 
 *
 * @param[in] ANS			    -   output.
 * @param[in] P			        -   input P in EFp2.
 */
extern void EFp2_ECD(EFp2 *ANS,EFp2 *P);

/**
 * ECA in EFp2 rational point 
 *
 * @param[in] ANS			        -   output.
 * @param[in] P1			        -   input P1 in EFp2.
 * @param[in] P2			        -   input P2 in EFp2.
 */
extern void EFp2_ECA(EFp2 *ANS,EFp2 *P1,EFp2 *P2);

/**
 * ECA in EFp2 rational point 
 *
 * @param[in] ANS			        -   output.
 * @param[in] P1			        -   input P1 in EFp2.
 * @param[in] scalar		        -   input scalar in mpz_t integer.
 */
extern void EFp2_SCM(EFp2 *ANS,EFp2 *P,mpz_t scalar);

#endif /* bn_efp2_h */
