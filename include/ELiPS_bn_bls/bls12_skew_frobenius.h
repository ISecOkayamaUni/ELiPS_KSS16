//
//  bls12_skew_frobenius.h
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
 * Interaface for BLS12's Skew Frobenius map implementations
 *
 * @ingroup bls12
 */

#ifndef bls12_skew_frobenius_h
#define bls12_skew_frobenius_h

#include <ELiPS_bn_bls/bn_bls12_precoms.h>
#include <ELiPS_bn_bls/bn_efp12.h>

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

/**
 * Calculate Skew Frobenius map for p^2 in G1 of BLS12 curve
 *
 * @param[out] ANS			    -  output in G1
 * @param[in] A			        -  A in EFp.
 */
extern void bls12_EFp_skew_frobenius_map_p2(EFp *ANS,EFp *A);

/**
 * Calculate Skew Frobenius map for pwer of p in G2 of BLS12 curve
 *
 * @param[out] ANS			    -  output in G2'
 * @param[in] A			        -  A in EFp2.
 */
extern void bls12_EFp2_skew_frobenius_map_p1(EFp2 *ANS,EFp2 *A);

/**
 * Calculate Skew Frobenius map for pwer of p^2 in G2 of BLS12 curve
 *
 * @param[out] ANS			    -  output in G2'
 * @param[in] A			        -  A in EFp2.
 */
extern void bls12_EFp2_skew_frobenius_map_p2(EFp2 *ANS,EFp2 *A);

/**
 * Calculate Skew Frobenius map for pwer of p^3 in G2 of BLS12 curve
 *
 * @param[out] ANS			    -  output in G2'
 * @param[in] A			        -  A in EFp2.
 */
extern void bls12_EFp2_skew_frobenius_map_p3(EFp2 *ANS,EFp2 *A);

/**
 * Calculate Skew Frobenius map for pwer of p^10 in G2 of BLS12 curve
 *
 * @param[out] ANS			    -  output in G2'
 * @param[in] A			        -  A in EFp2.
 */
extern void bls12_EFp2_skew_frobenius_map_p10(EFp2 *ANS,EFp2 *A);


#endif /* bls12_skew_frobenius_h */
