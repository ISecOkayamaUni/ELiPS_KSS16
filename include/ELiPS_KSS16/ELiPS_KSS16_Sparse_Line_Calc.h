//
//  ELiPS_KSS16_Sparse_Line_Calc.h
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#ifndef ELiPS_KSS16_Sparse_Line_Calc_h
#define ELiPS_KSS16_Sparse_Line_Calc_h

#include <ELiPS_KSS16/ELiPS_EFp16.h>

extern void sparse_add_line_kss16(struct Fp16 *l_ANS,
                                  struct EFp4 *T_ANS,
                                  struct EFp4 *T,
                                  struct EFp4 *P,
                                  struct EFp4 *Q,
                                  struct Fp4 *Px_neg);

extern void sparse_dbl_line_kss16(struct Fp16 *l_ANS,
                                  struct EFp4 *T_ANS,
                                  struct EFp4 *T,
                                  struct EFp4 *Q,
                                  struct Fp4 *Px_neg);
#endif /* ELiPS_KSS16_Sparse_Line_Calc_h */
