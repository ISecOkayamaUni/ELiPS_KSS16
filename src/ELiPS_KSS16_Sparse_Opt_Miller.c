//
//  KSS16_Sparse7_Opt_Miller.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_KSS16_Sparse_Opt_Miller.h>

void sparse_opt_miller_kss16(struct Fp16 *ANS, struct EFp4 *P, struct EFp4 *Q, mpz_t loop){
    struct EFp4 T, EFp4_tmp;
    EFp4_init (&T);
    EFp4_init (&EFp4_tmp);
    
    struct EFp4 Q_neg;
    EFp4_init(&Q_neg);
    Fp4_neg(&Q_neg.y,&Q->y);
    Fp4_set(&Q_neg.x,&Q->x);
    
    struct Fp16 l_sum, ltt, ltp;
    Fp16_init (&l_sum);
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    Fp_set_ui (&l_sum.x0.x0.x0.x0, 1);
    
    struct Fp4 Px_neg;
    Fp4_init(&Px_neg);
    Fp4_neg(&Px_neg,&P->x);
    
    if(x_signed_binary[x_bit]==-1){
        EFp4_set(&T,&Q_neg);
    }
    else{
        EFp4_set(&T,Q);
    }
    
    int i;
    for(i=x_bit-1;i>=0;i--){
        switch (x_signed_binary[i]){
            case 0:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                sparse_dbl_line_kss16(&ltt,&T,&T,P,&Px_neg);
                Fp16_mul(&l_sum,&l_sum,&ltt);
                break;
                
            case 1:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                
                sparse_dbl_line_kss16(&ltt,&T,&T,P,&Px_neg);
                sparse_add_line_kss16(&ltp,&T,&T,Q,P,&Px_neg);
                
                Fp16_mul(&l_sum,&l_sum,&ltt);
                Fp16_mul(&l_sum,&l_sum,&ltp);
                break;
            case -1:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                
                sparse_dbl_line_kss16(&ltt,&T,&T,P,&Px_neg);
                sparse_add_line_kss16(&ltp,&T,&T,&Q_neg,P,&Px_neg);
                
                Fp16_mul(&l_sum,&l_sum,&ltt);
                Fp16_mul(&l_sum,&l_sum,&ltp);
                break;
        }
    }
    
    kss16_skew_frobenius_map(&EFp4_tmp,Q);
    sparse_add_line_kss16(&ltp,&T,&T,&EFp4_tmp,P,&Px_neg);
    Fp16_mul(&l_sum,&l_sum,&ltp);
    
    
    struct Fp16 tmp_f;
    Fp16_init(&tmp_f);
    
    
    //    Fp16_pow(&tmp_f, &l_sum,prime);
    Fp16_frobenius_map(&tmp_f, &l_sum);
    //    Fp16_pow(&l_sum, &tmp_f,prime);
    Fp16_frobenius_map(&l_sum, &tmp_f);
    //    Fp16_pow(&tmp_f, &l_sum,prime);
    Fp16_frobenius_map(&tmp_f, &l_sum);
    
    
    //ltt_q_kss16(&ltt,Q,P);
    sparse_dbl_line_kss16(&ltt,&T,Q,P,&Px_neg);
    
    Fp16_mul(&l_sum,&tmp_f,&ltt);
    Fp16_set(ANS,&l_sum);
    
    
    Fp16_clear(&l_sum);
    EFp4_clear(&T);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
    EFp4_clear(&EFp4_tmp);
    
}
