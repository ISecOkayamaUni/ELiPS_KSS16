//
//  bls12_frobenius.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/5/18.

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
 * Header of the Frobenius mapping implementation for BLS12 over Fp12 extension field.
 *
 * @ingroup bls12
 */

#ifndef bls12_frobenius_h
#define bls12_frobenius_h

#include <ELiPS_bn_bls/bn_efp12.h>
#include <ELiPS_bn_bls/bn_bls12_precoms.h>

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

/**
 * Calculate Frobenius map of A in Fp12 extension field as A^p where p is the BLS12 prime.
 *
 * @param[out] ANS				- the result.
 * @param[in] A				    - the input vecotr.
 */
extern void bls12_Fp12_frobenius_map_p1(Fp12 *ANS,Fp12 *A);

/**
 * Calculate Frobenius map of A in Fp12 extension field as A^p^2 where p is the BLS12 prime.
 *
 * @param[out] ANS				- the result.
 * @param[in] A				    - the input vecotr.
 */
extern void bls12_Fp12_frobenius_map_p2(Fp12 *ANS,Fp12 *A);

/**
 * Calculate Frobenius map of A in Fp12 extension field as A^p^3 where p is the BLS12 prime.
 *
 * @param[out] ANS				- the result.
 * @param[in] A				    - the input vecotr.
 */
extern void bls12_Fp12_frobenius_map_p3(Fp12 *ANS,Fp12 *A);

/**
 * Calculate Frobenius map of A in Fp12 extension field as A^p^4 where p is the BLS12 prime.
 *
 * @param[out] ANS				- the result.
 * @param[in] A				    - the input vecotr.
 */
extern void bls12_Fp12_frobenius_map_p4(Fp12 *ANS,Fp12 *A);

/**
 * Calculate Frobenius map of A in Fp12 extension field as A^p^6 where p is the BLS12 prime.
 *
 * @param[out] ANS				- the result.
 * @param[in] A				    - the input vecotr.
 */
extern void bls12_Fp12_frobenius_map_p6(Fp12 *ANS,Fp12 *A);

/**
 * Calculate Frobenius map of A in Fp12 extension field as A^p^8 where p is the BLS12 prime.
 *
 * @param[out] ANS				- the result.
 * @param[in] A				    - the input vecotr.
 */
extern void bls12_Fp12_frobenius_map_p8(Fp12 *ANS,Fp12 *A);

/**
 * Calculate Frobenius map of A in Fp12 extension field as A^p^10 where p is the BLS12 prime.
 *
 * @param[out] ANS				- the result.
 * @param[in] A				    - the input vecotr.
 */
extern void bls12_Fp12_frobenius_map_p10(Fp12 *ANS,Fp12 *A);

#endif /* bls12_frobenius_h */
