//
//  Optimal_Millers_Algo.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_KSS16_Optimal_Millers_Algo.h>

void optimal_miller_algo_kss16(struct Fp16 *ANS,struct EFp16 *P,struct EFp16 *Q, mpz_t loop){
    
    struct EFp16 T,EFp16_tmp;
    EFp16_init(&T);
    EFp16_init(&EFp16_tmp);
    
    struct Fp16 l_sum;
    Fp16_init(&l_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    
    struct Fp16 Px_neg;
    Fp16_init(&Px_neg);
    Fp16_neg(&Px_neg,&P->x);//TODO Why neg Px?
    
    
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    struct EFp16 Q_neg;
    EFp16_init(&Q_neg);
    Fp16_neg(&Q_neg.y,&Q->y);
    Fp16_set(&Q_neg.x,&Q->x);
    
    if(x_signed_binary[x_bit]==-1){
        EFp16_set(&T,&Q_neg);
    }else{
        EFp16_set(&T,Q);
    }
    int i;
    for(i=x_bit-1;i>=0;i--){
        switch (x_signed_binary[i]){
            case 0:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                dbl_line_kss16(&ltt,&T,&T,P,&Px_neg);
                Fp16_mul(&l_sum,&l_sum,&ltt);
                break;
                
            case 1:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                
                dbl_line_kss16(&ltt,&T,&T,P,&Px_neg);
                add_line_kss16(&ltp,&T,&T,Q,P,&Px_neg);
                
                Fp16_mul(&l_sum,&l_sum,&ltt);
                Fp16_mul(&l_sum,&l_sum,&ltp);
                break;
            case -1:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                
                dbl_line_kss16(&ltt,&T,&T,P,&Px_neg);
                add_line_kss16(&ltp,&T,&T,&Q_neg,P,&Px_neg);
                
                Fp16_mul(&l_sum,&l_sum,&ltt);
                Fp16_mul(&l_sum,&l_sum,&ltp);
                break;
        }
    }
    
    //  EFp16_scm_bin(&EFp_tmp,Q,prime);
    struct EFp4 Q_bar,EFp4_tmp;
    EFp4_init(&Q_bar);
    EFp4_init(&EFp4_tmp);
    EFp16_to_EFp4_map(&Q_bar, Q);
    kss16_skew_frobenius_map(&EFp4_tmp,&Q_bar);
    EFp4_to_EFp16_map(&EFp16_tmp, &EFp4_tmp);
    //  EFp16_frobenius_map(&EFp_tmp, Q);
    
    
    ltp_q_kss16(&ltp,&T,&EFp16_tmp,P);
    Fp16_mul(&l_sum,&l_sum,&ltp);
    
    
    struct Fp16 tmp_f;
    Fp16_init(&tmp_f);
    
    
    //    Fp16_pow(&tmp_f, &l_sum,prime);
    Fp16_frobenius_map(&tmp_f, &l_sum);
    //    Fp16_pow(&l_sum, &tmp_f,prime);
    Fp16_frobenius_map(&l_sum, &tmp_f);
    //    Fp16_pow(&tmp_f, &l_sum,prime);
    Fp16_frobenius_map(&tmp_f, &l_sum);
    
    ltt_q_kss16(&ltt,Q,P);
    
    
    Fp16_mul(&l_sum,&tmp_f,&ltt);
    Fp16_set(ANS,&l_sum);
    
    
    Fp16_clear(&l_sum);
    EFp16_clear(&T);
    EFp16_clear(&EFp16_tmp);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
    EFp4_clear(&Q_bar);
    EFp4_clear(&EFp4_tmp);
}
