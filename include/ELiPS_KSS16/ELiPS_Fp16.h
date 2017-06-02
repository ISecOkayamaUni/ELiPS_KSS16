//
//  Fp16.h
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/26/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#ifndef ELiPS_Fp16_h
#define ELiPS_Fp16_h

#include <ELiPS_KSS16/ELiPS_Fp8.h>

struct Fp16{
    struct Fp8 x0,x1;
};

#pragma mark Fp16 methods declarations

extern void Fp16_init(struct Fp16 *A);

extern void Fp16_set(struct Fp16 *ANS,
                     struct Fp16 *A);

extern void Fp16_set_ui(struct Fp16 *A,
                        signed long int B);

extern void Fp16_random(struct Fp16 *A);

extern void Fp16_clear(struct Fp16 *A);

extern void Fp16_printf(struct Fp16 *A);

extern void Fp16_add(struct Fp16 *ANS,
                     struct Fp16 *A,
                     struct Fp16 *B);

extern void Fp16_add_ui(struct Fp16 *ANS,
                        struct Fp16 *A,
                        unsigned long int B);

extern void Fp16_sub(struct Fp16 *ANS,
                     struct Fp16 *A,
                     struct Fp16 *B);

extern void Fp16_mul(struct Fp16 *ANS,
                     struct Fp16 *A,
                     struct Fp16 *B);

extern void Fp16_mul_ui(struct Fp16 *ANS,
                        struct Fp16 *A,
                        unsigned long int B);

extern void Fp16_mul_mpz(struct Fp16 *ANS,
                         struct Fp16 *A,
                         mpz_t B);

extern void Fp16_mul_Fp(struct Fp16 *ANS,
                        struct Fp16 *A,
                        struct Fp *B);

extern void Fp16_neg(struct Fp16 *ANS,
                     struct Fp16 *A);

extern void Fp16_invert(struct Fp16 *ANS,
                        struct Fp16 *A);

extern void Fp16_div(struct Fp16 *ANS,
                     struct Fp16 *A,
                     struct Fp16 *B);

extern void Fp16_pow(struct Fp16 *ANS,
                     struct Fp16 *A,
                     mpz_t B);

extern void Fp16_sqrt(struct Fp16 *ANS,
                      struct Fp16 *A);//x^2=a mod p

extern int  Fp16_legendre(struct Fp16 *A);

extern int  Fp16_cmp(struct Fp16 *A,
                     struct Fp16 *B);

extern int  Fp16_cmp_mpz(struct Fp16 *A,
                         mpz_t B);

extern void Fp16_neg(struct Fp16 *ANS,
                     struct Fp16 *A);

extern void Fp16_frobenius_map(struct Fp16 *ANS,
                               struct Fp16 *A);

#endif /* ELiPS_Fp16_h */
