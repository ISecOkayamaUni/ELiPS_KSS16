//
//  bn_efp6.h
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

#ifndef bn_efp6_h
#define bn_efp6_h

#include <ELiPS_bn_bls/bn_efp2.h>
/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/
/**
 * Initialize an EFp6 type structure 
 *
 * @param[in] P			        -  input P in EFp6.
 */
extern void EFp6_init(EFp6 *P);

/**
 * Clears an EFp6 type structure 
 *
 * @param[in] P			        -  input P in EFp6.
 */
extern void EFp6_clear(EFp6 *P);

/**
 * Prints an EFp6 type structure 
 *
 * @param[in] P			        -  input P in EFp6.
 * @param[in] str			    -  any string to print.
 */
extern void EFp6_printf(EFp6 *P,char *str);

/**
 * Sets EFp6 point an EFp6 type structure 
 *
 * @param[in] ANS			    -   output.
 * @param[in] P			        -   input P in EFp6.
 */
extern void EFp6_set(EFp6 *ANS,EFp6 *P);

/**
 * Sets unsigned int an EFp6 type structure 
 *
 * @param[in] ANS			    -   output.
 * @param[in] u			        -   input u in us integer.
 */
extern void EFp6_set_ui(EFp6 *ANS,unsigned long int UI);

/**
 * Sets mpz_t in an EFp6 type structure 
 *
 * @param[in] ANS			    -   output.
 * @param[in] str			    -   input A in mpz_t int.
 */
extern void EFp6_set_mpz(EFp6 *ANS,mpz_t A);

/**
 * Negate an EFp6 type structure 
 *
 * @param[in] ANS			    -   output.
 * @param[in] P			        -   input P in EFp6.
 */
extern void EFp6_set_neg(EFp6 *ANS,EFp6 *P);

/**
 * Generate an EFp6 rational point
 *
 * @param[in] P			        -   input rational point P and.
 */
extern void EFp6_rational_point(EFp6 *P);

/**
 * ECD in EFp6 rational point 
 *
 * @param[in] ANS			    -   output.
 * @param[in] P			        -   input P in EFp6.
 */
extern void EFp6_ECD(EFp6 *ANS,EFp6 *P);

/**
 * ECA in EFp6 rational point 
 *
 * @param[in] ANS			        -   output.
 * @param[in] P1			        -   input P1 in EFp6.
 * @param[in] P2			        -   input P2 in EFp6.
 */
extern void EFp6_ECA(EFp6 *ANS,EFp6 *P1,EFp6 *P2);

/**
 * SCM in EFp6 rational point 
 *
 * @param[in] ANS			        -   output.
 * @param[in] P1			        -   input P1 in EFp6.
 * @param[in] scalar		        -   input scalar in mpz_t integer.
 */
extern void EFp6_SCM(EFp6 *ANS,EFp6 *P,mpz_t scalar);



#endif /* bn_efp6_h */
