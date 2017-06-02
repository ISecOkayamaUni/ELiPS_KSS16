//
//  EFp.h
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/26/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#ifndef ELiPS_EFp_h
#define ELiPS_EFp_h

#include <ELiPS_KSS16/ELiPS_Fp.h>

struct EFp{
    struct Fp x,y;
    int PoI;
};

#pragma mark EFp methods
extern void EFp_init(struct EFp *A);
extern void EFp_set(struct EFp *A,struct EFp *B);
extern void EFp_set_PoI(struct EFp *A);
extern void EFp_clear(struct EFp *A);
extern void EFp_printf(struct EFp *A);
extern void EFp_scm_bin(struct EFp *ANS, struct EFp *P,mpz_t j);
extern void EFp_ecd(struct EFp *ANS, struct EFp *P);//ANS=2*P
extern void EFp_eca(struct EFp *ANS, struct EFp *P1, struct EFp *P2);//ANS=P1+P2
extern int  EFp_cmp(struct EFp *A,struct EFp *B);
extern void EFp_random_set(struct EFp *ANS);//random set EFp on curve
extern void EFp_neg(struct EFp *ANS, struct EFp *A);
extern void rational_point_check(struct EFp *A);
#endif /* ELiPS_EFp_h */
