//
//  millers_algo_kss16.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_KSS16_Millers_Algo.h>

void millers_algo_kss16(struct Fp16 *ANS,struct EFp16 *P,struct EFp16 *Q, mpz_t loop){
    struct Fp16 l_sum,v_sum;
    Fp16_init(&l_sum);
    Fp16_init(&v_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    Fp_set_ui(&v_sum.x0.x0.x0.x0,1);
    
    struct EFp16 T;
    EFp16_init(&T);
    EFp16_set(&T,P);
    
    
    struct Fp16 ltt,ltp,v2t,vtp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    Fp16_init(&v2t);
    Fp16_init(&vtp);
    
    int i;
    struct Fp16 tmp1;
    Fp16_init(&tmp1);
    //    Fp16_init(&lambda);
    int r_bit;//bit数
    
    r_bit= (int)mpz_sizeinbase(loop,2);
    
    for(i=r_bit-2;i>=0;i--){
        Fp16_mul(&l_sum,&l_sum,&l_sum);
        Fp16_mul(&v_sum,&v_sum,&v_sum);
        
        ltt_q_kss16(&ltt,&T,Q);
        Fp16_mul(&l_sum,&l_sum,&ltt);
        
        EFp16_ecd(&T,&T);
        v2t_q_kss16(&v2t,&T,Q);
        Fp16_mul(&v_sum,&v_sum,&v2t);
        
        if(mpz_tstbit(loop,i)==1){
            ltp_q_kss16(&ltp,&T,P,Q);
            Fp16_mul(&l_sum,&l_sum,&ltp);
            
            EFp16_eca(&T,&T,P);
            vtp_q_kss16(&vtp,&T,Q);
            Fp16_mul(&v_sum,&v_sum,&vtp);
        }
    }
    
    
    // EFp16_printf(&T);
    Fp16_div(ANS,&l_sum,&v_sum);
    
    Fp16_clear(&l_sum);
    Fp16_clear(&v_sum);
    EFp16_clear(&T);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
    Fp16_clear(&v2t);
    Fp16_clear(&vtp);
    Fp16_clear(&tmp1);
}
