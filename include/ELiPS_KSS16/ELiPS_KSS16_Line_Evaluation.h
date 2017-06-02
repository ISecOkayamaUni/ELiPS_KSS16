//
//  Line_evaluation.h
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#ifndef ELiPS_KSS16_Line_Evaluation_h
#define ELiPS_KSS16_Line_Evaluation_h

#include <ELiPS_KSS16/ELiPS_EFp16.h>

extern void ltt_q_kss16(struct Fp16  *ANS,
                        struct EFp16 *T,
                        struct EFp16 *Q);

extern void v2t_q_kss16(struct Fp16  *ANS,
                        struct EFp16 *T,
                        struct EFp16 *Q);

extern void ltp_q_kss16(struct Fp16  *ANS,
                        struct EFp16 *T,
                        struct EFp16 *P,
                        struct EFp16 *Q);

extern void vtp_q_kss16(struct Fp16  *ANS,
                        struct EFp16 *T,
                        struct EFp16 *Q);

extern void add_line_kss16(struct Fp16  *l_ANS,
                           struct EFp16 *T_ANS,
                           struct EFp16 *T,
                           struct EFp16 *P,
                           struct EFp16 *Q,
                           struct Fp16  *Qx_neg);

extern void dbl_line_kss16(struct Fp16  *l_ANS,
                           struct EFp16 *T_ANS,
                           struct EFp16 *T,
                           struct EFp16 *Q,
                           struct Fp16  *Qx_neg);

#endif /* ELiPS_KSS16_Line_Evaluation_h */
