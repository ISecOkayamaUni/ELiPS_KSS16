//
//  fp8.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/26/18.

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

#ifndef fp8_h
#define fp8_h

#include <ELiPS_bn_bls/fp4.h>

void Fp8_init( Fp8 *A);
void Fp8_set( Fp8 *ANS, Fp8 *A);
void Fp8_set_ui( Fp8 *A,signed long int B);
void Fp8_random( Fp8 *A);
void Fp8_clear( Fp8 *A);
void Fp8_printf( Fp8 *A);
void Fp8_add( Fp8 *ANS, Fp8 *A, Fp8 *B);
void Fp8_add_ui( Fp8 *ANS, Fp8 *A,unsigned long int B);
void Fp8_sub( Fp8 *ANS, Fp8 *A, Fp8 *B);
void Fp8_mul( Fp8 *ANS, Fp8 *A, Fp8 *B);
void Fp8_sqr( Fp8 *ANS, Fp8 *A);
void Fp8_sqr_complex( Fp8 *ANS, Fp8 *A);
void Fp8_mul_v( Fp8 *ANS, Fp8 *A);
void Fp8_mul_ui( Fp8 *ANS, Fp8 *A,unsigned long int B);
void Fp8_mul_mpz( Fp8 *ANS, Fp8 *A,mpz_t B);
void Fp8_mul_Fp( Fp8 *ANS, Fp8 *A, Fp *B);
void Fp8_neg( Fp8 *ANS, Fp8 *A);
void Fp8_invert( Fp8 *ANS, Fp8 *A);
void Fp8_div( Fp8 *ANS, Fp8 *A, Fp8 *B);
void Fp8_pow( Fp8 *ANS, Fp8 *A,mpz_t B);
void Fp8_sqrt( Fp8 *ANS, Fp8 *A);//x^2=a mod p
int  Fp8_legendre( Fp8 *A);
int  Fp8_cmp( Fp8 *A, Fp8 *B);
int  Fp8_cmp_mpz( Fp8 *A,mpz_t B);

#endif /* fp8_h */
