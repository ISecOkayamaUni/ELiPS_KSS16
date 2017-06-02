//
//  KSS16_Opt_Ate.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_KSS16_Opt_Ate.h>
#include <ELiPS_KSS16/ELiPS_KSS16_Optimal_Millers_Algo.h>
#include <ELiPS_KSS16/ELiPS_KSS16_Final_Exp.h>

void opt_ate_kss16(struct Fp16 *ANS,struct EFp16 *G1,struct EFp16 *G2){
    struct Fp16 Miller_X, t_ans;
    Fp16_init(&Miller_X);
    Fp16_init(&t_ans);
    
    optimal_miller_algo_kss16(&Miller_X,G1,G2,params.X);
    
    final_exp_kss16(&t_ans,&Miller_X);
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&Miller_X);
    Fp16_clear(&t_ans);
}
