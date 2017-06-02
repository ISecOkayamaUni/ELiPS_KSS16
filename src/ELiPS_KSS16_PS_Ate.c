//
//  KSS16_PS_Ate.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_KSS16_PS_Ate.h>
#include <ELiPS_KSS16/ELiPS_KSS16_PS_Miller.h>
#include <ELiPS_KSS16/ELiPS_KSS16_Final_Exp.h>

void ps_ate_kss16(struct Fp16 *ANS,struct EFp *G1,struct EFp16 *G2){
    struct Fp16 t_ans;
    Fp16_init(&t_ans);
    struct EFp4 EFp4_G1,EFp4_G2;
    EFp4_init(&EFp4_G1);
    EFp4_init(&EFp4_G2);
    mpz_t tm1;
    mpz_init(tm1);
    mpz_sub_ui(tm1,params.trace_t,1);
    
    EFp4_set_EFp(&EFp4_G1,G1);//get P of FP in Fp4
    EFp16_to_EFp4_map(&EFp4_G2,G2); //non isomorphic map of Q to Q'
    ps_miller_kss16(&t_ans,&EFp4_G1,&EFp4_G2,tm1);
    
    final_exp_kss16(&t_ans,&t_ans);
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&t_ans);
    mpz_clear(tm1);
}
