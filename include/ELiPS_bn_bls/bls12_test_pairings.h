//
//  bls12_test_pairings.h
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
 * Interaface for BLS12's Usefull pairing and scm functionality
 *
 * @ingroup bls12
 */

#ifndef bls12_test_pairings_h
#define bls12_test_pairings_h



#include <ELiPS_bn_bls/bls12_pairings.h>
#include <ELiPS_bn_bls/bls12_generate_points.h>
#include <ELiPS_bn_bls/bls12_timeprint.h>
#include <ELiPS_bn_bls/bls12_scm.h>
#include <ELiPS_bn_bls/bls12_G3_exp.h>

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

/**
 * Test tate pairing in BLS12 curve
 */
extern void bls12_test_tate_pairing(void);

/**
 * Test ate pairing in BLS12 curve
 */
extern void bls12_test_plain_ate_pairing(void);

/**
 * Test opt-ate pairing in BLS12 curve
 */
extern void bls12_test_opt_ate_pairing(void);

/**
 * Test G1 scm  in BLS12 curve
 */
extern void bls12_test_G1_scm(void);

/**
 * Test G2 scm  in BLS12 curve
 */
extern void bls12_test_G2_scm(void);

/**
 * Test G3 exponentiation in BLS12 curve
 */
extern void bls12_test_G3_exp(void);

#endif /* bls12_test_pairings_h */
