/*
 * ELiPS is an Efficient Library for Pairing-based Systems
 * Copyright (C) 2008-2017 ELiPS Authors. Please refer to the COPYRIGHT file
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


/** @file Settings.h
 *  @brief Brief description about the source file/header file.
 *
 *  Then here write about the small details what this source file is intended for.
 *  If it calculated prime field arithmetic then write what are the major funcion
 *  ality it provides though the public functions.  For example what type of
 *  public funtions it will give or global variables you will need.
 *
 *  @author Khandaker Md. Al-Amin(khandaker@s.okayama-u.ac.jp)
 *  @date 4/26/17
 *  @bug No know bug yet
 */

#ifndef ELiPS_KSS16_Settings_h
#define ELiPS_KSS16_Settings_h

/* -- Includes -- */
/* inlcudes 1st level header of project configuration*/
#include <ELiPS_KSS16/Common_Headers.h>

/* inlcudes system time header*/
#include <sys/time.h>

/*============================================================================*/
/* Macro definitions                                                          */
/*============================================================================*/

/**
 * @brief Defines TRUE or FALSE using integer
 *
 * The TRUE_1 defines true as int 1 when boolean decesion in true.
 * FALSE_0 defines the oppsite
 */
#define TRUE_1 1
#define FALSE_0 0

/**
 * @brief Defines the constant values of KSS16 parameters that are pre-calculated.
 *
 * Fixed the prime number, group order r, Frobenius trace t, number of rational point EFp
 */
#define PRIME "615623382030675150502066218751443438064107566348210118507940234835256709422634902533028653925239565581"
#define ORDER_R "10897499371578763791778093615151768824360936005521891580808300080405508061745073"
#define TRACE_T "1403565040305261127593483468292232497718027806453082"
#define ORDER_EFp "615623382030675150502066218751443438064107566348208714942899929574129115939166610300530935897433112500"

/**
 * @brief maximum degree of mother parameter X's polynomial
 * 
 * 2^35-2^32-2^18+2^8+1 is the polynomial used to find X.
 */
#define x_bit 35

/*============================================================================*/
/* Data structure declarations                                                */
/*============================================================================*/

/**
 * @brief Structure to hold the KSS16 curve settings.
 *
 * #KSS16_params hold the KSS16 curves parametets.
 */
struct kss16_params {
    mpz_t prime;        /**< Prime number or Characteristics of KSS16 curve */
    mpz_t order_r;      /**< Group order r of KSS16 curve */
    mpz_t trace_t;      /**< Frobenius trace of KSS16 curve */
    mpz_t order_EFp;    /**< Number of rational point of KSS16 curve over the Prime field Fp*/
    mpz_t X;            /**< Mother parameter of KSS16 curve*/
};

/**
 * @brief KSS16 curve's constant values.
 *
 * KSS16 curve y^2 = x^3+ a * x has one constant value a. But at first we need to get the value of a by putting 
 * arbitary values of x.
 */
struct KSS16_constants{
    mpz_t a;            /**< Constant a in KSS16 pairing-friedly elliptic curve y^2 = x^3+ a * x KSS16 curve*/
    mpz_t tmp_a;        /**< Temoporary variable for constant a*/
};

/**
 * @brief KSS16 curve's systematically obtained parameters.
 *
 * It's a global variable that give access to KSS16 public parameters.
 */
extern struct kss16_params params;

/**
 * @brief KSS16 curve's constants.
 *
 * It's a global variable that give access to KSS16 constant term a.
 */
extern struct KSS16_constants kss_curve_const;

/**
 * @brief The qudratic non residue over prime filed Fp.
 *
 * X^2-c1 is irreducible over the Fp
 */
extern int c1;

/**
 * @brief Checks if all the parameters of KSS16 curve is initialize.
 *
 */
extern int is_params_set;

extern int TRUE;
extern int FALSE;

/**
 * @brief Chracter arry of size [x_bit+1] that holds the singned binary representation of the mother parameter X.
 * 
 * This x_signed_binary is important since it is used in the Miller's algorithm for efficeint calculation using the singned 
 * binary representation of the mother parameter X.
 */
extern char x_signed_binary [x_bit+1];

/*============================================================================*/
/* Function declarations                                                      */
/*============================================================================*/

/**
 * @brief Initialize and sets the parameters of KSS16 curve y^2=x^3+1x.
 *
 */
extern void set_kss16_params (void);

/**
 * @brief Initialize and sets the constant a = 1 of KSS16 curve y^2=x^3+1x.
 *
 */
extern void set_kss16_curve_const (void);

/**
 * @brief Generate mother parameter X=2^35-2^32-2^18+2^8+1 given in http://eprint.iacr.org/2017/334
 *
 */
extern void generate_mother_parameter (void);

/**
 * @brief Utility function to calculate the time difference.
 *
 */
extern float timedifference_msec (struct timeval t0, struct timeval t1);

#endif /* ELiPS_KSS16_Settings_h */
