//
//  curve_settings.h
//
//  Created by Khandaker Md. Al-Amin on 1/25/18.
//  Copyright © 2018 Khandaker Md. Al-Amin. All rights reserved.
//
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
 * @defgroup elips Curve parameter generation and settings
 */

/**
 * @file
 *
 * Interface for curve parameters.
 *
 * @ingroup elips
 */

#ifndef bn_settings_h
#define bn_settings_h

#include <ELiPS_bn_bls/Commont_headers.h>
//#include "bn_fp.h"
//#include "bn_fp2.h"

extern int isCleared; // 1 is cleared. 0 is not cleared

/*============================================================================*/
/* Struct definitions                                                          */
/*============================================================================*/


/**
 * @brief Curve Parameters
 *
 * This curve_params structure contains member variable related to pairing friendly curves.
 */
struct curve_params {
    mpz_t prime;        /**< Prime number genereted using polynomial formula of X curve_params#prime. */
    mpz_t X;            /**< Mother parameter curve_params#X. */
    int sign;
    mpz_t trace_t;       /**< Curves  Frobenius trace  curve_params#trace_t. */
    mpz_t order;         /**< Curves  order r given by polynomial formula of X  curve_params#order. */
    mpz_t EFp_total;     /**<Number of rational points in curve defined in Fp  curve_params#EFp_total. */
    mpz_t EFp2_total;    /**<Number of rational points in curve defined in Fp2  curve_params#EFp2_total. */
    mpz_t EFp6_total;    /**<Number of rational points in curve defined in Fp6  curve_params#EFp6_total. */
    mpz_t EFp12_total;   /**<Number of rational points in curve defined in Fp12  curve_params#EFp12_total. */
    mpz_t EFpd_total;    /**<Number of rational points in twisted curve defined in Fp^k/d  curve_params#EFp_total. */
    mpz_t curve_a;       /**<Curves coeffient y^2=x^3+ax+b curve_params#a */
    mpz_t curve_b;       /**<Curves coeffient y^2=x^3+ax+b curve_params#b */
};
/**
 * @brief bn curve's systematically obtained parameters.
 *
 * It's a global variable that give access to bn public parameters.
 */
extern struct curve_params curve_parameters;

/**
 *  BN curve's mother parameter length in bit.
 *
 */
#define bn_X_length 114

/**
 * Global random state to generate random element
 *
 */
extern gmp_randstate_t state;

/**
 * Arrary to hold mother parameter of BN curve
 *
 */
extern char X_binary[bn_X_length+1];

/**
 * Arrary to hold loop length Miller's algo for opt-ate pairing of BN curve
 *
 */
extern char X_binary_opt[bn_X_length+3];

//bls12
/**
 *  BLS12 curve's mother parameter length in bit.
 *
 */
extern int bls12_X_length;

extern mpz_t bls12_X;
/**
 *  Array to hold BLS curve's mother parameter.
 *
 */
extern int bls12_X_binary[78];



/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/


//init
/**
 *  Initialize BN curves settings.
 *
 */
extern void init_bn_settings(void);

/**
 *  Initialize BN curves parameter generation.
 *
 */
extern void init_bn_parameters( void);

/**
 *  Generate BN curves mother parameter X.
 *
 */
extern void generate_bn_mother_parameter(void);

/**
 *  Generate BN curve prime.
 *
 */
extern int  generate_bn_prime(void);

/**
 *  Generate BN curve order.
 *
 */
extern int  generate_bn_order(void);

/**
 *  Generate BN curve Frobenius trace.
 *
 */
extern void generate_bn_trace(void);

/**
 *  Set BN curve coefficients.
 *
 */
extern void set_bn_curve_parameter(void);

/**
 *  Generate BN curve's  number of elements in EFp.
 *
 */
extern void weil(void);

/**
 *  Prints curves public parameters.
 *
 */
extern void print_curve_parameters(void);


/**
 *  Initialize BLS12 curve settings.
 *
 */
extern void init_bls12_settings(void);

/**
 *  Initialize BLS12 curve parameters.
 *
 */
extern void init_bls12_parameters( void);

/**
 *  Initialize BLS12 curve mother parameter.
 *
 */
extern void bls12_generate_X(void);

/**
 *  Initialize BLS12 curve prime.
 *
 */
extern int  bls12_generate_prime(void);

/**
 *  Initialize BLS12 curve order.
 *
 */
extern int  bls12_generate_order(void);

/**
 *  Initialize BLS12 curve Frobenius trace.
 *
 */
extern void bls12_generate_trace(void);

/**
 *  Sets BLS12 curve coefficients b.
 *
 */
extern void bls12_set_curve_parameter(void);

/**
 *  Calculate BLS12 curve field size EFp.
 *
 */
extern void bls12_weil(void);

/**
 *  Prints BLS12 curve parameters.
 *
 */
extern void bls12_print_parameters(void);


//KSS-16
#define TRUE 1
#define FALSE 0
#define KSS16_X_length 35

extern char X_bit_binary_kss16[KSS16_X_length+1];

extern mpz_t C1_INV; //c=2 its iverse value
//mpz_t PRIME_P,order_r,trace_t, order_EFp, a_x;
//mpz_t tmp_a;

extern void init_kss16_settings(void);
extern void generate_kss16_motherparam(void);
extern void generate_kss16_parameters(void);
extern void init_kss16_parameters(void);
#endif /* bn_settings_h */
