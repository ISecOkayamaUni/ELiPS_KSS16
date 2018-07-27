//
//  bn_fp6.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 1/31/18.

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

#ifndef bn_fp6_h
#define bn_fp6_h

#include <ELiPS_bn_bls/bn_fp2.h>


extern void Fp6_init(Fp6 *A);
extern void Fp6_clear(Fp6 *A);
extern void Fp6_printf(Fp6 *A,char *str);
extern void Fp6_set(Fp6 *ANS,Fp6 *A);
extern void Fp6_set_ui(Fp6 *ANS,unsigned long int UI);
extern void Fp6_set_mpz(Fp6 *ANS,mpz_t A);
extern void Fp6_set_neg(Fp6 *ANS,Fp6 *A);
extern void Fp6_set_random(Fp6 *ANS,gmp_randstate_t state);
extern void Fp6_mul(Fp6 *ANS,Fp6 *A,Fp6 *B);
extern void Fp6_mul_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI);
extern void Fp6_mul_mpz(Fp6 *ANS,Fp6 *A,mpz_t B);
extern void Fp6_mul_basis(Fp6 *ANS,Fp6 *A);
extern void Fp6_sqr(Fp6 *ANS,Fp6 *A);
extern void Fp6_add(Fp6 *ANS,Fp6 *A,Fp6 *B);
extern void Fp6_add_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI);
extern void Fp6_add_mpz(Fp6 *ANS,Fp6 *A,mpz_t B);
extern void Fp6_sub(Fp6 *ANS,Fp6 *A,Fp6 *B);
extern void Fp6_sub_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI);
extern void Fp6_sub_mpz(Fp6 *ANS,Fp6 *A,mpz_t B);
extern void Fp6_inv(Fp6 *ANS,Fp6 *A);
extern int  Fp6_legendre(Fp6 *A);
extern int  Fp6_isCNR(Fp6 *A);
extern void Fp6_sqrt(Fp6 *ANS,Fp6 *A);
extern void Fp6_pow(Fp6 *ANS,Fp6 *A,mpz_t scalar);
extern int  Fp6_cmp(Fp6 *A,Fp6 *B);
extern int  Fp6_cmp_ui(Fp6 *A,unsigned long int UI);
extern int  Fp6_cmp_mpz(Fp6 *A,mpz_t B);
extern int  Fp6_cmp_zero(Fp6 *A);
extern int  Fp6_cmp_one(Fp6 *A);
#endif /* bn_fp6_h */
