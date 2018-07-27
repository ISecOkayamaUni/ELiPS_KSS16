//
//  bls12_line_ate.h
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
 * Interaface for BLS12 line evaluation for ate pairing.
 *
 * @ingroup bls12
 */

#ifndef bls12_line_ate_h
#define bls12_line_ate_h

#include <ELiPS_bn_bls/bls12_p8sparse.h>

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/
/**
 * Calculate line equation ltt that goes through T and P in Miller's algo.
 *
 * @param[out] f			    - output in Fp12 
 * @param[in] T			        - Input rational point T in Fp2.
 * @param[in] P		            - Input rational point P in mapped Fp.
 * @param[in] L			        - Input element L in mapped Fp after precomputation in Miller's algo.
 */
extern void bls12_ff_ltt(Fp12 *f,EFp2 *T,EFp *P,Fp *L);

/**
 * Calculate line equation ltq that goes through T,Q and P in Miller's algo.
 *
 * @param[out] f			    - output in Fp12 
 * @param[in] T			        - Input rational point T in Fp2.
 * @param[in] Q		            - Input rational point Q mapped in Fp2.
 * @param[in] P		            - Input rational point P mapped in Fp.
 * @param[in] L			        - Input element L in mapped Fp after precomputation in Miller's algo.
 */
extern void bls12_f_ltq(Fp12 *f,EFp2 *T,EFp2 *Q,EFp *P,Fp *L);

#endif /* bls12_line_ate_h */
