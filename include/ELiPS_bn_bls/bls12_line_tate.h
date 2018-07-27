//
//  bls12_line_tate.h
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
 * Interaface for BLS12 line evaluation for tate pairing.
 *
 * @ingroup bls12
 */

#ifndef bls12_line_tate_h
#define bls12_line_tate_h

#include <ELiPS_bn_bls/bn_efp12.h>

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/
/**
 * Calculate line equation for tate pairing in BLS12 curve.
 *
 * @param[out] f			    - output in Fp12 
 * @param[in] T			        - Input rational point T in Fp.
 * @param[in] Q		            - Input rational point Q mapped in Fp22.
 */
extern void bls12_ff_ltt_vtt_for_tate(Fp12 *f,EFp *T,EFp12 *Q);

/**
 * Calculate line equation for tate pairing in BLS12 curve.
 *
 * @param[out] f			    - output in Fp12 
 * @param[in] T			        - Input rational point T in Fp.
 * @param[in] P			        - Input rational point P in Fp.
 * @param[in] Q		            - Input rational point Q mapped in Fp22.
 */
extern void bls12_f_ltp_vtp_for_tate(Fp12 *f,EFp *T,EFp *P,EFp12 *Q);

#endif /* bls12_line_tate_h */
