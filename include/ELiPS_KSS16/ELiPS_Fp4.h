//
//  Fp4.h
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/26/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#ifndef ELiPS_Fp4_h
#define ELiPS_Fp4_h

#include <ELiPS_KSS16/ELiPS_Fp2.h>

/**
 
 */
struct Fp4{
    struct Fp2 x0,x1;
};

extern struct Fp4 z_inv2;

#pragma mark Fp4 methods

/**
 * @brief Initializer for Fp4 vector
 * 
 * A Fp4 vector pointer is passed to the function which is initialized to all zero coefficients.
 *
 * @param[in] A Fp4 vector pointer as input to initialize.
 */
extern void Fp4_init(struct Fp4 *A);


extern void Fp4_set(struct Fp4 *ANS,
                    struct Fp4 *A);

extern void Fp4_set_ui(struct Fp4 *A,
                       signed long int B);

extern void Fp4_random(struct Fp4 *A);

extern void Fp4_clear(struct Fp4 *A);

extern void Fp4_printf(struct Fp4 *A);

extern void Fp4_add(struct Fp4 *ANS,
                    struct Fp4 *A,
                    struct Fp4 *B);

extern void Fp4_add_ui(struct Fp4 *ANS,
                       struct Fp4 *A,
                       unsigned long int B);

extern void Fp4_sub(struct Fp4 *ANS,
                    struct Fp4 *A,
                    struct Fp4 *B);

extern void Fp4_mul(struct Fp4 *ANS,
                    struct Fp4 *A,
                    struct Fp4 *B);

extern void Fp4_mul_v(struct Fp4 *ANS,
                      struct Fp4 *A);

extern void Fp4_mul_ui(struct Fp4 *ANS,
                       struct Fp4 *A,
                       unsigned long int B);

extern void Fp4_mul_mpz(struct Fp4 *ANS,
                        struct Fp4 *A,
                        mpz_t B);

extern void Fp4_mul_Fp(struct Fp4 *ANS,
                       struct Fp4 *A,
                       struct Fp *B);

extern void Fp4_neg(struct Fp4 *ANS,
                    struct Fp4 *A);

extern void Fp4_invert(struct Fp4 *ANS,
                       struct Fp4 *A);

extern void Fp4_div(struct Fp4 *ANS,
                    struct Fp4 *A,
                    struct Fp4 *B);

extern void Fp4_pow(struct Fp4 *ANS,
                    struct Fp4 *A,
                    mpz_t B);

extern void Fp4_sqrt(struct Fp4 *ANS,
                     struct Fp4 *A);//x^2=a mod p

extern int  Fp4_legendre(struct Fp4 *A);

extern int  Fp4_cmp(struct Fp4 *A,
                    struct Fp4 *B);

extern int  Fp4_cmp_mpz(struct Fp4 *A,
                        mpz_t B);

extern void Fp4_neg(struct Fp4 *ANS,
                    struct Fp4 *A);

extern void Fp4_frobenius_map(struct Fp4 *ANS,
                              struct Fp4 *A);

extern void Fp4_mul_beta_inv(struct Fp4 *ANS);
#endif /* ELiPS_Fp4_h */
