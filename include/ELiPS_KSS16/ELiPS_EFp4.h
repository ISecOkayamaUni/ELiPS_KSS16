//
//  EFp4.h
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#ifndef EFp4_h
#define EFp4_h

#include <ELiPS_KSS16/ELiPS_Fp4.h>
#include <ELiPS_KSS16/ELiPS_EFp.h>

struct EFp4{
    struct Fp4 x,y;
    int PoI;
};

#pragma mark EFp4 methods declarations
extern void EFp4_init(struct EFp4 *A);
extern void EFp4_set(struct EFp4 *A,struct EFp4 *B);
extern void EFp4_set_PoI(struct EFp4 *A);
extern void EFp4_set_EFp(struct EFp4 *ANS,struct EFp *A);
extern void EFp4_clear(struct EFp4 *A);
extern void EFp4_printf(struct EFp4 *A);
extern void EFp4_ecd(struct EFp4 *ANS, struct EFp4 *P);//ANS=2*P
extern void EFp4_eca(struct EFp4 *ANS, struct EFp4 *P1, struct EFp4 *P2);//ANS=P1+P2
extern int  EFp4_cmp(struct EFp4 *A,struct EFp4 *B);
extern void EFp4_scm_bin(struct EFp4 *ANS, struct EFp4 *P, mpz_t j);
extern void EFp4_scm_win(struct EFp4 *ANS, struct EFp4 *P, mpz_t scalar);
extern void EFp4_scm_ml(struct EFp4 *RES, struct EFp4 *P,mpz_t scalar);
extern void EFp4_neg(struct EFp4 *ANS, struct EFp4 *A);
extern void EFp4_scm_bin_sparse(struct EFp4 *ANS,struct EFp4 *P,mpz_t j);
extern void EFp4_ecd_sparse(struct EFp4 *ANS, struct EFp4 *P);
extern void kss16_skew_frobenius_map(struct EFp4 *ANS, struct EFp4 *Qt);

//extern void EFp4_random_set(struct EFp4 *ANS);
//extern void EFp4_SCM_BIN_Pseudo_Sparse(struct EFp4 *ANS,struct EFp4 *P,mpz_t j);
//extern void EFp4_ECD_Pseudo_Sparse(struct EFp4 *ANS, struct EFp4 *P);
#endif /* EFp4_h */
