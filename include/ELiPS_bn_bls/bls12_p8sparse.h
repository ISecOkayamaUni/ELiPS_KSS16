//
//  bls12_p8sparse.h
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
 * Interaface for BLS12 Pseudo 8-sparse multplication in Miller's algo of tate pairing https://eprint.iacr.org/2017/1174.pdf
 *
 * @ingroup bls12
 */

#ifndef bls12_p8sparse_h
#define bls12_p8sparse_h

#include <ELiPS_bn_bls/bn_efp2.h>
/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

/**
 * Calculate Precomputations for Pseudo 8-sparse multplication and Millers algo.
 *
 * @param[in] P			        - P in G1 
 * @param[in] Q			        - Q in G2.
 * @param[in] L		            - L in Fp from Miller's algo.
 */
extern void bls12_Pseudo_8_sparse_mapping(EFp *P,EFp2 *Q,Fp *L);

/**
 * Calculate Precomputations for Pseudo 8-sparse multplication and Millers algo.
 *
 * @param[out] ANS			    - output of Pseudo 8-sparse multplication
 * @param[in] A			        - Fp12 vector.
 * @param[in] B		            - Fp12 vector.
 */
extern void bls12_Pseudo_8_sparse_mul(Fp12 *ANS,Fp12 *A,Fp12 *B);

#endif /* bls12_p8sparse_h */
