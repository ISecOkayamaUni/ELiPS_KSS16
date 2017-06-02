//
//  EFp16.h
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#ifndef ELiPS_EFp16_h
#define ELiPS_EFp16_h

#include <ELiPS_KSS16/ELiPS_Fp16.h>
#include <ELiPS_KSS16/ELiPS_EFp.h>
#include <ELiPS_KSS16/ELiPS_EFp4.h>

struct EFp16{
    struct Fp16 x,y;
    int PoI;
};

extern void EFp16_init(struct EFp16 *A);

extern void EFp16_set(struct EFp16 *A,
                      struct EFp16 *B);

extern void EFp16_set_poi(struct EFp16 *A);

extern void EFp16_set_EFp(struct EFp16 *A,
                          struct EFp *B);

extern void EFp16_clear(struct EFp16 *A);

extern void EFp16_printf(struct EFp16 *A);

extern void EFp16_ecd(struct EFp16 *ANS,
                      struct EFp16 *P);//ANS=2*P

extern void EFp16_eca(struct EFp16 *ANS,
                      struct EFp16 *P1,
                      struct EFp16 *P2);//ANS=P1+P2

extern int  EFp16_cmp(struct EFp16 *A,
                      struct EFp16 *B);

extern void EFp16_scm_bin(struct EFp16 *ANS,
                          struct EFp16 *P,
                          mpz_t j);

extern void EFp16_scm_win(struct EFp16 *ANS,
                          struct EFp16 *P,
                          mpz_t scalar);

extern void EFp16_frobenius_map(struct EFp16 *ANS,
                                struct EFp16 *A);

extern void EFp16_random_set(struct EFp16 *ANS);

extern void EFp16_random_set_G2(struct EFp16 *ANS);

extern void EFp16_to_EFp4_map(struct EFp4 *ANS,
                              struct EFp16 *A);

extern void EFp4_to_EFp16_map(struct EFp16 *ANS,
                              struct EFp4 *A);

#endif /* ELiPS_EFp16_h */
