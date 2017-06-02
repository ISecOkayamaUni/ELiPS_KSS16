//
//  Fp2.h
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/26/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#ifndef ELiPS_Fp2
#define ELiPS_Fp2

#include <ELiPS_KSS16/ELiPS_Fp.h>

struct Fp2{
    struct Fp x0,x1;
};

#pragma mark Fp2 method delaration

extern void Fp2_init(struct Fp2 *A);

extern void Fp2_set(struct Fp2 *ANS,
                    struct Fp2 *A);

extern void Fp2_set_ui(struct Fp2 *A,
                       signed long int B);

extern void Fp2_random(struct Fp2 *A);

extern void Fp2_clear(struct Fp2 *A);

extern void Fp2_printf(struct Fp2 *A);

extern void Fp2_add(struct Fp2 *ANS,
                    struct Fp2 *A,
                    struct Fp2 *B);

extern void Fp2_add_ui(struct Fp2 *ANS,
                       struct Fp2 *A,
                       unsigned long int B);

extern void Fp2_sub(struct Fp2 *ANS,
                    struct Fp2 *A,
                    struct Fp2 *B);

extern void Fp2_mul(struct Fp2 *ANS,
                    struct Fp2 *A,
                    struct Fp2 *B);

extern void Fp4_mul_basis(struct Fp2 *ANS,
                      struct Fp2 *A);

extern void Fp2_mul_ui(struct Fp2 *ANS,
                       struct Fp2 *A,
                       unsigned long int B);

extern void Fp2_mul_mpz(struct Fp2 *ANS,
                        struct Fp2 *A,
                        mpz_t B);

extern void Fp2_mul_Fp(struct Fp2 *ANS,
                       struct Fp2 *A,
                       struct Fp *B);

extern void Fp2_neg(struct Fp2 *ANS,
                    struct Fp2 *A);

extern void Fp2_invert(struct Fp2 *ANS,
                       struct Fp2 *A);

extern void Fp2_div(struct Fp2 *ANS,
                    struct Fp2 *A,
                    struct Fp2 *B);

extern void Fp2_pow(struct Fp2 *ANS,
                    struct Fp2 *A,
                    mpz_t B);

extern void Fp2_sqrt(struct Fp2 *ANS,
                     struct Fp2 *A);//x^2=a mod p

extern int  Fp2_legendre(struct Fp2 *A);

extern int  Fp2_cmp(struct Fp2 *A,
                    struct Fp2 *B);

extern int  Fp2_cmp_mpz(struct Fp2 *A,
                        mpz_t B);

extern void Fp2_neg(struct Fp2 *ANS,
                    struct Fp2 *A);

extern void Fp2_frobenius_map(struct Fp2 *ANS,
                              struct Fp2 *A);

#endif /* ELiPS_Fp2 */
