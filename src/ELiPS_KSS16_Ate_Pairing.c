//
//  Ate_Pairing.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_KSS16_Ate_Pairing.h>
#include <ELiPS_KSS16/ELiPS_KSS16_Millers_Algo.h>
#include <ELiPS_KSS16/ELiPS_KSS16_Final_Exp.h>

void ate_pairing_kss16(struct Fp16 *ANS,struct EFp16 *G1,struct EFp16 *G2){
    struct Fp16 t_ans;
    Fp16_init(&t_ans);
    
    mpz_t tm1;
    mpz_init(tm1);
    mpz_sub_ui(tm1,params.trace_t,1);
    
    millers_algo_kss16(&t_ans,G2,G1,tm1);
    final_exp_kss16(&t_ans,&t_ans);
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&t_ans);
}
