//
//  KSS16_Sparse7_Miller.h
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#ifndef ELiPS_KSS16_Sparse_Miller_h
#define ELiPS_KSS16_Sparse_Miller_h

#include <ELiPS_KSS16/ELiPS_KSS16_Sparse_Line_Calc.h>

extern void sparse_miller_kss16(struct Fp16 *ANS,
                                struct EFp4 *P,
                                struct EFp4 *Q,
                                mpz_t loop);

#endif /* ELiPS_KSS16_Sparse_Miller_h */
