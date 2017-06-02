//
//  tate_pairing_kss16.h
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#ifndef ELiPS_KSS16_Tate_Pairing_h
#define ELiPS_KSS16_Tate_Pairing_h

#include <ELiPS_KSS16/ELiPS_KSS16_Line_Evaluation.h>

extern void tate_pairing_kss16(struct Fp16 *ANS,
                               struct EFp16 *P,
                               struct EFp16 *Q);
#endif /* ELiPS_KSS16_Tate_Pairing_h */
