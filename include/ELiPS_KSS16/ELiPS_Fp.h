//
//  ELiPS_Fp.h
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/26/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#ifndef ELiPS_Fp
#define ELiPS_Fp

#include <ELiPS_KSS16/ELiPS_KSS16_Settings.h>

struct Fp{
    mpz_t x0;
};

extern struct Fp qnr_c1,c_inv,pm1d4,pm5d8,pm13d16;
extern struct Fp m_pm1d4,m_cpm13d16;
extern struct Fp m_cpm5d8,m_cpm1d4pm5d8,pm1d4pm5d8,cpm1d4pm5d8,pm5d8pm13d16,m_cpm5d8pm13d16,cpm5d8pm13d16,ccpm5d8pm13d16,m_ccpm1d4pm5d8p13d16,cpm1d4pm5d6pm13d16,ccpm1d4pm5d8pm13d16;
extern mpz_t p8p1dr;


extern void Fp_init(struct Fp *A);

extern void Fp_set(struct Fp *ANS,
                   struct Fp *A);

extern void Fp_set_ui(struct Fp *A,
                      signed long int B);

extern void Fp_random(struct Fp *A);

extern void Fp_clear(struct Fp *A);

extern void Fp_printf(struct Fp *A);

extern void Fp_add(struct Fp *ans,
                   struct Fp *a,
                   struct Fp *b);//ans=a+b mod p

extern void Fp_add_ui(struct Fp *ans,
                      struct Fp *a,
                      unsigned long int b);//ans=a+b mod p

extern void Fp_add_mpz(struct Fp *ans,
                       struct Fp *a,
                       mpz_t b);//ans=a+b mod p

extern void Fp_sub(struct Fp *ans,
                   struct Fp *a,
                   struct Fp *b);//ans=a-b mod p

extern void Fp_sub_ui(struct Fp *ans,
                      struct Fp *a,
                      unsigned long int b);//ans=a+b mod p

extern void Fp_mul(struct Fp *ans,
                   struct Fp *a,
                   struct Fp *b);//ans=a*b mod p

extern void Fp_mul_ui(struct Fp *ans,
                      struct Fp *a,
                      unsigned long int b);//ans=a*b mod p

extern void Fp_div(struct Fp *ans,
                   struct Fp *a,
                   struct Fp *b);//ans=a/b mod p

extern void Fp_pow(struct Fp *ans,
                   struct Fp *a,
                   mpz_t b);

extern void Fp_sqrt(struct Fp *ans,
                    struct Fp *a);//x^2=a mod p

extern int  Fp_cmp_mpz(struct Fp *A,
                       mpz_t B);

extern void Fp_mul_mpz(struct Fp *ANS,
                       struct Fp *A,
                       mpz_t B);

extern void Fp_neg(struct Fp *ANS,
                   struct Fp *A);

extern int  Fp_cmp(struct Fp *A,
                   struct Fp *B);

#pragma mark util methods declaration
extern void dealloc_constants (void);
extern void pre_calculate (void);

#endif /* ELiPS_Fp */
