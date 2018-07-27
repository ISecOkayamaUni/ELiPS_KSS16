//
//  bn_final_exp.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/2/18.

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
 * Interaface for final exponentiation in BN curve 
 *
 * @ingroup bn12
 */
#ifndef bn_final_exp_h
#define bn_final_exp_h

#include <ELiPS_bn_bls/bn_frobenius.h>
/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

/**
 * Un-optimized final exp for BN curve
 *
 * @param[in] P			        -  input P in EFp12.
 */
extern void bn_final_exp_plain(Fp12 *ANS,Fp12 *A);

/**
 *   Optimized final exp for BN curve
 *
 * @param[out] ANS			    -  output in Fp12.
 * @param[in] A			        -  input P in Fp12 as obtained from Miller's algo output.
 */
extern void bn_final_exp_optimal(Fp12 *ANS,Fp12 *A);

/**
 * Efficient exponentiation by mother parameter 
 *
 * @param[out] ANS			    -  output in Fp12.
 * @param[in] A			        -  input P in EFp12.
 */
extern void bn_fp12_power_motherparam(Fp12 *ANS,Fp12 *A);

#endif /* bn_final_exp_h */
