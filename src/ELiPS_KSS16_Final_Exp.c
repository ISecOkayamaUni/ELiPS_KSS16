//
//  KSS16_Final_Exp.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_KSS16_Final_Exp.h>

void final_exp_kss16(struct Fp16 *ANS,struct Fp16 *A){
    struct Fp16 temp_Fp16,temp_A;
    Fp16_init(&temp_Fp16);
    Fp16_init(&temp_A);
    Fp16_set(&temp_A, A);
    
    Fp16_frobenius_map(&temp_Fp16,&temp_A);
    Fp16_frobenius_map(&temp_A,&temp_Fp16);
    Fp16_frobenius_map(&temp_Fp16,&temp_A);
    Fp16_frobenius_map(&temp_A,&temp_Fp16);
    Fp16_frobenius_map(&temp_Fp16,&temp_A);
    Fp16_frobenius_map(&temp_A,&temp_Fp16);
    Fp16_frobenius_map(&temp_Fp16,&temp_A);
    Fp16_frobenius_map(&temp_A,&temp_Fp16);
    
    Fp16_div(&temp_Fp16, &temp_A, A);
    Fp16_pow(ANS,&temp_Fp16,p8p1dr);
    
    Fp16_clear(&temp_Fp16);
    Fp16_clear(&temp_A);
}
