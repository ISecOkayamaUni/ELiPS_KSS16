//
//  bn_fp2.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 1/29/18.

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
 * Interaface of degree 2 extension field (Fp2) operations
 *
 * @ingroup ff
 */

#ifndef bn_fp2_h
#define bn_fp2_h


#include <ELiPS_bn_bls/bn_fp.h>
#include <ELiPS_bn_bls/field_dtype.h>
#include <ELiPS_bn_bls/bn_bls12_precoms.h>

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

/**
 * Initializes an Fp2 type struct
 *
 * @param[in] A			        -  input A is a pointer of Fp2 type struct.
 */
extern void Fp2_init(Fp2 *A);

/**
 * Clears an Fp2 type struct from memory
 *
 * @param[in] A			        -  input A is a pointer of Fp2 type struct.
 */
extern void Fp2_clear(Fp2 *A);

/**
 * Prints an Fp2 type struct
 *
 * @param[in] A			        -  input A is a pointer of Fp2 type struct.
 * @param[in] str			    -  input str is a string pointer.
 */
extern void Fp2_printf(Fp2 *A,char *str);

/**
 * Sets an Fp2 type struct to another Fp2 struct
 *
 * @param[out] ANS			-  input A is a pointer of Fp2 type struct.
 * @param[in] A			    -  input A is a pointer of Fp2 type struct.
 */
extern void Fp2_set(Fp2 *ANS,Fp2 *A);

/**
 * Sets an unsigned long int type element to another Fp2 struct
 *
 * @param[out] ANS			-  input A is a pointer of Fp2 type struct.
 * @param[in] A			    -  input A is a unsigned long integer.
 */
extern void Fp2_set_ui(Fp2 *ANS,unsigned long int A);

/**
 * Sets an mpz_t type integer to another Fp2 struct. 
 * It calls the similar functions in Fp for each elements in Fp.
 *
 * @param[out] ANS			-  input A is a pointer of Fp2 type struct.
 * @param[in] A			    -  input str is a string pointer.
 */
extern void Fp2_set_mpz(Fp2 *ANS,mpz_t A);

/**
 * Negate A and sets it in ANS as a Fp2 vector. It calls Fp_set_neg fuctions for each elements.
 * 
 *
 * @param[out] ANS			    -  set 
 * @param[in] A                 -  send A as Fp.
 */
extern void Fp2_set_neg(Fp2 *ANS,Fp2 *A);

/**
 * Sets a random Fp2 element in ANS where each Fp coefficient < prime
 *
 * @param[out] ANS			    -  output an random Fp2 vector.
 * @param[in] state                 -  send a gmp_randstate_t.
 */
extern void Fp2_set_random(Fp2 *ANS,gmp_randstate_t state);

/**
 * Fp2 extension field multiplication with modular reduction for.
 * Modular polynomial reduction x^2-c whre c has no sqare root in prime field Fp, 
 * needs to be defined before to use Fp_mul_basis function. 
 *
 * @param[out] ANS			    -  output ANS as pinter in Fp2.
 * @param[in] ANS			    -  send ANS as pinter in Fp2. 
 *                                 It will be returned as result of multiplication.
 * @param[in] A                 -  send pointer A in Fp2.
 * @param[in] B                 -  send pointer B in Fp2.
 */
extern void Fp2_mul(Fp2 *ANS,Fp2 *A,Fp2 *B);
extern void Fp2_mul_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);
extern void Fp2_mul_mpz(Fp2 *ANS,Fp2 *A,mpz_t B);
extern void Fp2_mul_Fp(Fp2 *ANS, Fp2 *A, Fp *B);
extern void Fp2_mul_basis(Fp2 *ANS,Fp2 *A);
extern void Fp2_inv_basis(Fp2 *ANS,Fp2 *A);
extern void Fp2_sqr(Fp2 *ANS,Fp2 *A);
extern void Fp2_add(Fp2 *ANS,Fp2 *A,Fp2 *B);
extern void Fp2_add_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);
extern void Fp2_add_mpz(Fp2 *ANS,Fp2 *A,mpz_t B);
extern void Fp2_sub(Fp2 *ANS,Fp2 *A,Fp2 *B);
extern void Fp2_sub_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI);
extern void Fp2_sub_mpz(Fp2 *ANS,Fp2 *A,mpz_t B);
extern void Fp2_inv(Fp2 *ANS,Fp2 *A);
extern int  Fp2_legendre(Fp2 *A);
extern int  Fp2_isCNR(Fp2 *A);
extern void Fp2_sqrt(Fp2 *ANS,Fp2 *A);
extern void Fp2_pow(Fp2 *ANS,Fp2 *A,mpz_t scalar);
extern int  Fp2_cmp(Fp2 *A,Fp2 *B);
extern int  Fp2_cmp_ui(Fp2 *A,unsigned long int UI);
extern int  Fp2_cmp_mpz(Fp2 *A,mpz_t B);
extern int  Fp2_cmp_zero(Fp2 *A);
extern int  Fp2_cmp_one(Fp2 *A);
extern void Fp2_neg(Fp2 *ans,Fp2 *a);

extern void Fp2_mul_KSS16(Fp2 *ANS,Fp2 *A,Fp2 *B);
extern void Fp2_mul_basis_KSS16(Fp2 *ANS, Fp2 *A);
extern void Fp2_invert_kss16(Fp2 *ANS,Fp2 *A);


// #pragma mark Fp2 methods
//void Fp2_sqr(struct Fp2 *ANS,struct Fp2 *A);
//void Fp2_sqr_kss16(struct Fp2 *ANS,struct Fp2 *A);
//void Fp2_sqr_complex(struct Fp2 *ANS,struct Fp2 *A);
//void Fp2_mul_i(struct Fp2 *ANS,struct Fp2 *A);
//void Fp2_mul_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int B);
//void Fp2_mul_mpz(struct Fp2 *ANS,struct Fp2 *A,mpz_t B);
//void Fp2_mul_Fp(struct Fp2 *ANS,struct Fp2 *A,struct Fp *B);
//void Fp2_neg(struct Fp2 *ANS,struct Fp2 *A);
extern void Fp2_inv(Fp2 *ANS, Fp2 *A);
//void Fp2_div(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B);
//void Fp2_pow(struct Fp2 *ANS,struct Fp2 *A,mpz_t B);
//void Fp2_sqrt(struct Fp2 *ANS,struct Fp2 *A);//x^2=a mod p
//int  Fp2_legendre(struct Fp2 *A);
//int  Fp2_cmp(struct Fp2 *A,struct Fp2 *B);
//int  Fp2_cmp_mpz(struct Fp2 *A,mpz_t B);
//void Fp2_neg(struct Fp2 *ANS,struct Fp2 *A);

#endif /* bn_fp2_h */
