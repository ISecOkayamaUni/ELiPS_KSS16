//
//  EFp2.h
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#ifndef ELiPS_EFp2_h
#define ELiPS_EFp2_h

#include <ELiPS_KSS16/ELiPS_EFp.h>
#include <ELiPS_KSS16/ELiPS_Fp2.h>

struct EFp2{
    struct Fp2 x,y;
    int PoI;
};


// #pragma mark EFp2 methods
extern void EFp2_init(struct EFp2 *A);
extern void EFp2_set(struct EFp2 *A,struct EFp2 *B);
extern void EFp2_set_PoI(struct EFp2 *A);
extern void EFp2_clear(struct EFp2 *A);
extern void EFp2_printf(struct EFp2 *A);
extern void EFp2_ecd(struct EFp2 *ANS, struct EFp2 *P);//ANS=2*P
extern void EFp2_eca(struct EFp2 *ANS, struct EFp2 *P1, struct EFp2 *P2);//ANS=P1+P2
extern int  EFp2_cmp(struct EFp2 *A,struct EFp2 *B);
extern void EFp2_set_EFp(struct EFp2 *A,struct EFp *B);
extern void EFp2_random_set(struct EFp2 *ANS);
extern void EFp2_scm_bin(struct EFp2 *ANS, struct EFp2 *P, mpz_t j);
extern void EFp2_scm_naf(struct EFp2 *ANS, struct EFp2 *P, mpz_t scalar);
extern void EFp2_scm_win(struct EFp2 *ANS, struct EFp2 *P, mpz_t scalar);
extern void EFp2_neg(struct EFp2 *ANS, struct EFp2 *A);
extern void EFp2_rational_point_check(struct EFp2 *A);
#endif /* ELiPS_EFp2_h */
