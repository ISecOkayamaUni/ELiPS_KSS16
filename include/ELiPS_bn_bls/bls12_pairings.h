//
//  bls12_pairings.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/6/18.

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
 * Interaface for BLS12 Pairing implementations
 *
 * @ingroup bls12
 */

#ifndef bls12_pairings_h
#define bls12_pairings_h

#include <ELiPS_bn_bls/bls12_miller_ate.h>
#include <ELiPS_bn_bls/bls12_miller_tate.h>
#include <ELiPS_bn_bls/bls12_miller_optate.h>
#include <ELiPS_bn_bls/bls12_finalexp.h>
#include <ELiPS_bn_bls/bls12_timeprint.h>
#include <ELiPS_bn_bls/bn_utils.h>

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

/**
 * Calculate tate pairing in BLS12 curve.
 *
 * @param[out] ANS			    -  output 
 * @param[in] P			        - P in  G1 mapped to Fp12.
 * @param[in] Q		            - Q in G2.
 */
extern void bls12_tate(Fp12 *ANS,EFp12 *P, EFp12 *Q);

/**
 * Calculate ate pairing in BLS12 curve.
 *
 * @param[out] ANS			    -  output 
 * @param[in] P			        - P in  G1 mapped to Fp12.
 * @param[in] Q		            - Q in G2.
 */
extern void bls12_plain_ate(Fp12 *ANS,EFp12 *P,EFp12 *Q);

/**
 * Calculate optimal ate pairing in BLS12 curve.
 *
 * @param[out] ANS			    -  output 
 * @param[in] P			        - P in  G1 mapped to Fp12.
 * @param[in] Q		            - Q in G2.
 */
extern void bls12_opt_ate(Fp12 *ANS,EFp12 *P, EFp12 *Q);

#endif /* bls12_pairings_h */
