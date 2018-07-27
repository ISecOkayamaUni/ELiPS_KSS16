//
//  bn_fp.h
//  bn_Single_File_Elips
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
 * @defgroup ff Finite field construction
 */

/**
 * @file
 *
 * Interaface of prime field operations. Primarily targeted for BN and BLS12 curve
 *
 * @ingroup ff
 */

#ifndef bn_fp_h
#define bn_fp_h

#include <ELiPS_bn_bls/curve_settings.h>
#include <ELiPS_bn_bls/field_dtype.h>

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/

/**
 * Initializes an Fp type struct
 *
 * @param[in] A			        -  input A is a pointer of Fp type struct.
 */
extern void Fp_init(Fp *A);

/**
 * Clears memory an Fp type struct
 *
 * @param[in] A			        -  input A is a pointer of Fp type struct.
 */
extern void Fp_clear(Fp *A);

/**
 * Prints an Fp type struct
 *
 * @param[out] A			        -  input A is a pointer of Fp type struct.
 * @param[in] str               -  any string pointer.
 */
extern void Fp_printf(Fp *A,char *str);

/**
 * Set an Fp type struct to ANS
 *
 * @param[out] ANS			    -  set A in ANS.
 * @param[in] A                 -  send A as pointer.
 */
extern void Fp_set(Fp *ANS,Fp *A);

/**
 * Set an unsigned int in ANS Fp type struct.
 *
 * @param[out] ANS			    -  set A in ANS.
 * @param[in] A                 -  send A as UI.
 */
extern void Fp_set_ui(Fp *ANS,unsigned long int A);

/**
 * Set an gmp int mpz in ANS Fp type struct.
 *
 * @param[out] ANS			    -  set A in ANS.
 * @param[in] A                 -  send A as mpz_t.
 */
extern void Fp_set_mpz(Fp *ANS,mpz_t B);

/**
 * Negate A and set it in ANS
 *
 * @param[out] ANS			    -  set prime-A in ANS.
 * @param[in] A                 -  send A as Fp.
 */
extern void Fp_set_neg(Fp *ANS,Fp *A);

/**
 * Set random Fp of size <prime in ANS
 *
 * @param[out] ANS			    -  output an Fp < prime.
 * @param[in] A                 -  send a gmp_randstate_t.
 */
extern void Fp_set_random(Fp *ANS,gmp_randstate_t state);

/**
 * Prime field multiplication with reduction as ANS= A*B mod prime
 *
 * @param[out] ANS			    -  output ANS < prime in Fp.
 * @param[in] A                 -  send pointer A in Fp.
 * @param[in] B                 -  send pointer B in Fp.
 */
extern void Fp_mul(Fp *ANS,Fp *A,Fp *B);

/**
 * Prime field multiplication with reduction as ANS= A*B mod prime
 *
 * @param[out] ANS			    -  output ANS < prime in Fp.
 * @param[in] A                 -  send pointer A in Fp.
 * @param[in] B                 -  send unsigned int.
 */
extern void Fp_mul_ui(Fp *ANS,Fp *A,unsigned long int B);

/**
 * Prime field multiplication with reduction as ANS= A*B mod prime
 *
 * @param[out] ANS			    -  output ANS < prime in Fp.
 * @param[in] A                 -  send pointer A in Fp.
 * @param[in] B                 -  send mpz_t .
 */
extern void Fp_mul_mpz(Fp *ANS,Fp *A,mpz_t B);

/**
 * Prime field multiplication with reduction as ANS= A*B mod prime
 *
 * @param[out] ANS			    -  output ANS < prime in Fp.
 * @param[in] A                 -  send pointer A in Fp.
 * @param[in] B                 -  send basis element \alpha^2 = some int
 */
extern void Fp_mul_basis(Fp *ANS,Fp *A);

/**
 * Prime field multiplication with reduction as ANS= A*B mod prime
 *
 * @param[out] ANS			    -  output ANS < prime in Fp.
 * @param[in] A                 -  send pointer A in Fp.
 * @param[in] B                 -  send basis element \alpha^2 = c is a int
 */
extern void Fp_mul_basis_KSS16(Fp *ANS,Fp *A);

/**
 * Prime field addition with reduction as ANS= A+B mod prime
 *
 * @param[out] ANS			    -  output ANS < prime in Fp.
 * @param[in] A                 -  send pointer A in Fp.
 * @param[in] B                 -  send pointer B in Fp.
 */
extern void Fp_add(Fp *ANS,Fp *A,Fp *B);

/**
 * Prime field addition with reduction as ANS= A+B mod prime
 *
 * @param[out] ANS			    -  output ANS < prime in Fp.
 * @param[in] A                 -  send pointer A in Fp.
 * @param[in] B                 -  send unsigned int B.
 */
extern void Fp_add_ui(Fp *ANS,Fp *A,unsigned long int B);

/**
 * Prime field addition with reduction as ANS= A+B mod prime
 *
 * @param[out] ANS			    -  output ANS < prime in Fp.
 * @param[in] A                 -  send pointer A in Fp.
 * @param[in] B                 -  send mpz_t int B.
 */
extern void Fp_add_mpz(Fp *ANS,Fp *A,mpz_t B);

/**
 * Prime field subtraction with reduction as ANS= A-B mod prime
 *
 * @param[out] ANS			    -  output ANS < prime in Fp.
 * @param[in] A                 -  send pointer A in Fp.
 * @param[in] B                 -  send pointer  B.
 */
extern void Fp_sub(Fp *ANS,Fp *A,Fp *B);

/**
 * Prime field subtraction with reduction as ANS= A-B mod prime
 *
 * @param[out] ANS			    -  output ANS < prime in Fp.
 * @param[in] A                 -  send pointer A in Fp.
 * @param[in] B                 -  send  int B.
 */
extern void Fp_sub_ui(Fp *ANS,Fp *A,unsigned long int B);

/**
 * Prime field subtraction with reduction as ANS = A-B mod prime
 *
 * @param[out] ANS			    -  output ANS < prime in Fp.
 * @param[in] A                 -  send pointer A in Fp.
 * @param[in] B                 -  send mpz_t int B.
 */
extern void Fp_sub_mpz(Fp *ANS,Fp *A,mpz_t B);

/**
 * Prime field inversion A^-1 mod prime
 *
 * @param[out] ANS			    -  output A^-1 mod prime in Fp.
 * @param[in] A                 -  send pointer A in Fp.
 */
extern void Fp_inv(Fp *ANS,Fp *A);

/**
 * Calculate legendre symbol  (A/prime) to determine A has any sqare root in prime field.
 *
 * @param[out] int			    -  return 1 if A has QR and -1 for QNR. 
 * @param[in] A                 -  send pointer A in Fp.
 */
extern int  Fp_legendre(Fp *A);

/**
 * Check if A has qubic root or not in prime field.
 *
 * @param[out] int			    -  return 1 if A has CR and -1  CNR. 
 * @param[in] A                 -  send pointer A in Fp.
 */
extern int  Fp_isCNR(Fp *A);

/**
 * Calculate  sqrt(A) mod prime using Toneli-Shanks alogrithm in prime field.
 *
 * @param[out] ANS			    -  return 1 if A has CR and -1  CNR. 
 * @param[in] A                 -  send pointer A in Fp.
 */
extern void Fp_sqrt(Fp *ANS,Fp *A);

/**
 * Calculate  A^mpz_t mod prime in prime field.
 *
 * @param[out] ANS			    -  return 1 if A has CR and -1  CNR. 
 * @param[in] A                 -  send pointer A in Fp.
 */
extern void Fp_pow(Fp *ANS,Fp *A,mpz_t scalar);

/**
 * Compares A and B in Fp
 *
 * @param[out] int			    -  return 0 if A==B and 1 otherwise. 
 * @param[in] A                 -  send pointer A in Fp.
 * @param[in] B                 -  send pointer B in Fp.
 */
extern int  Fp_cmp(Fp *A,Fp *B);

/**
 * Compares A and B in Fp
 *
 * @param[out] int			    -  return 0 if A==B and 1 otherwise. 
 * @param[in] A                 -  send pointer A in Fp.
 * @param[in] B                 -  send pointer B in int.
 */
extern int  Fp_cmp_ui(Fp *A,unsigned long int UI);

/**
 * Compares A and B in Fp
 *
 * @param[out] int			    -  return 0 if A==B and 1 otherwise. 
 * @param[in] A                 -  send pointer A in Fp.
 * @param[in] B                 -  send pointer B in mpz_t.
 */
extern int  Fp_cmp_mpz(Fp *A,mpz_t B);

/**
 * Compares A == 0 in Fp
 *
 * @param[out] int			    -  return 0 if A==0 and 1 otherwise. 
 * @param[in] A                 -  send pointer A in Fp.
 * @param[in] B                 -  send pointer B in int.
 */
extern int  Fp_cmp_zero(Fp *A);

/**
 * Compares A == 1 in Fp
 *
 * @param[out] int			    -  return 0 if A==1 and 1 otherwise. 
 * @param[in] A                 -  send pointer A in Fp.
 */
extern int  Fp_cmp_one(Fp *A);

/**
 *  Calculate -A in Fp
 *
 * @param[out] ANS			    -  return -A. 
 * @param[in] A                 -  send pointer A in Fp.
 */
extern void Fp_neg(struct Fp *ANS,struct Fp *A);

#endif /* bn_fp_h */
