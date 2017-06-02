//
//  Fp8.h
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/26/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#ifndef ELiPS_Fp8_h
#define ELiPS_Fp8_h

#include <ELiPS_KSS16/ELiPS_Fp4.h>

struct Fp8{
    struct Fp4 x0,x1;
};

#pragma mark Fp8 methods declaration
extern void Fp8_init(struct Fp8 *A);
extern void Fp8_set(struct Fp8 *ANS,struct Fp8 *A);
extern void Fp8_set_ui(struct Fp8 *A,signed long int B);
extern void Fp8_random(struct Fp8 *A);
extern void Fp8_clear(struct Fp8 *A);
extern void Fp8_printf(struct Fp8 *A);
extern void Fp8_add(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B);
extern void Fp8_add_ui(struct Fp8 *ANS,struct Fp8 *A,unsigned long int B);
extern void Fp8_sub(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B);
extern void Fp8_mul(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B);
extern void Fp8_mul_v(struct Fp8 *ANS,struct Fp8 *A);
extern void Fp8_mul_ui(struct Fp8 *ANS,struct Fp8 *A,unsigned long int B);
extern void Fp8_mul_mpz(struct Fp8 *ANS,struct Fp8 *A,mpz_t B);
extern void Fp8_mul_Fp(struct Fp8 *ANS,struct Fp8 *A,struct Fp *B);
extern void Fp8_neg(struct Fp8 *ANS,struct Fp8 *A);
extern void Fp8_invert(struct Fp8 *ANS,struct Fp8 *A);
extern void Fp8_div(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B);
extern void Fp8_pow(struct Fp8 *ANS,struct Fp8 *A,mpz_t B);
extern void Fp8_sqrt(struct Fp8 *ANS,struct Fp8 *A);//x^2=a mod p
extern int  Fp8_legendre(struct Fp8 *A);
extern int  Fp8_cmp(struct Fp8 *A,struct Fp8 *B);
extern int  Fp8_cmp_mpz(struct Fp8 *A,mpz_t B);
extern void Fp8_frobenius_map(struct Fp8 *ANS, struct Fp8 *A);
#endif /* ELiPS_Fp8_h */
