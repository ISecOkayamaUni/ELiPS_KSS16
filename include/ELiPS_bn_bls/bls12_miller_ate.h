//
//  bls12_miller_ate.h
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
 * Interaface for BLS12 Miller's algo for unoptimized version of ate pairing.
 *
 * @ingroup bls12
 */
#ifndef bls12_miller_ate_h
#define bls12_miller_ate_h

#include <ELiPS_bn_bls/bls12_twist.h>
#include <ELiPS_bn_bls/bls12_p8sparse.h>
#include <ELiPS_bn_bls/bls12_line_ate.h>

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/
/**
 * Calculate Miller's loop for ate pairing in BLS12 curve.
 *
 * @param[out] ANS			    - output in Fp12
 * @param[in] Q			        - Input rational point Q in G2.
 * @param[in] P		            - Input rational point P in G1.
 */
extern void bls12_Miller_algo_for_plain_ate(Fp12 *ANS,EFp12 *Q,EFp12 *P);


#endif /* bls12_miller_ate_h */
