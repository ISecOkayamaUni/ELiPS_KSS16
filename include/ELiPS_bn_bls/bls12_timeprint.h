//
//  bls12_timeprint.h
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
 * Interaface for BLS12's time profiling 
 *
 * @ingroup bls12
 */

#ifndef bls12_timeprint_h
#define bls12_timeprint_h

#include <stdio.h>
/*============================================================================*/
/* Macro definitions                                                          */
/*============================================================================*/

/**
 * Time types of different operations.
 */
extern double bls12_MILLER_TATE,bls12_MILLER_PLAINATE,bls12_MILLER_OPTATE;
extern double bls12_FINALEXP_PLAIN,bls12_FINALEXP_OPT;
extern double bls12_G1SCM_PLAIN,bls12_G1SCM_2SPLIT;
extern double bls12_G2SCM_PLAIN,bls12_G2SCM_2SPLIT,bls12_G2SCM_4SPLIT;
extern double bls12_G3EXP_PLAIN,bls12_G3EXP_2SPLIT,bls12_G3EXP_4SPLIT;

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

/**
 * Prints curve parameter: BLS12 curve
 */
extern void bls12_print_parameters(void);

/**
 * Prints G1 rational point: BLS12 curve
 */
extern void bls12_print_G1_point(void);

/**
 * Prints G2 rational point: BLS12 curve
 */
extern void bls12_print_G2_point(void);

/**
 * Prints Tate paring time: BLS12 curve
 */
extern void bls12_print_tate_time(void);

/**
 * Prints ate paring time: BLS12 curve
 */
extern void bls12_print_plain_ate_time(void);


/**
 * Prints opt-ate paring time: BLS12 curve
 */
extern void bls12_print_opt_ate_time(void);

/**
 * Prints plain final exp time: BLS12 curve
 */
extern void bls12_print_final_exp_plain_time(void);

/**
 * Prints efficient final exp time: BLS12 curve
 */
extern void bls12_print_final_exp_optimal_time(void);

/**
 * Prints plain G1 scm time: BLS12 curve
 */
extern void bls12_print_plain_G1_scm_time(void);

/**
 * Prints 2 split of scalr G1 scm time: BLS12 curve
 */
extern void bls12_print_2split_G1_scm_time(void);

/**
 * Prints plain G2 scm time: BLS12 curve
 */
extern void bls12_print_plain_G2_scm_time(void);

/**
 * Prints 2 split of scalr G2 scm time: BLS12 curve
 */
extern void bls12_print_2split_G2_scm_time(void);
/**
 * Prints 4 split of scalr G2 scm time: BLS12 curve
 */
extern void bls12_print_4split_G2_scm_time(void);

/**
 * Prints no split of scalr G3 exponentiation time: BLS12 curve
 */
extern void bls12_print_plain_G3_exp_time(void);

/**
 * Prints 2  split of scalr G3 exponentiation time: BLS12 curve
 */
extern void bls12_print_2split_G3_exp_time(void);

/**
 * Prints 2  split of scalr G3 exponentiation time: BLS12 curve
 */
extern void bls12_print_4split_G3_exp_time(void);

#endif /* bls12_timeprint_h */
