//
//  KSS16_Sparse7_Miller.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_KSS16_Sparse_Miller.h>

void sparse_miller_kss16(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop){
    struct Fp16 l_sum;
    Fp16_init(&l_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    
    struct EFp4 T,EFp_tmp;
    EFp4_init(&T);
    EFp4_init(&EFp_tmp);
    
    struct Fp4 Px_neg;
    Fp4_init(&Px_neg);
    Fp4_neg(&Px_neg,&P->x);
    
    EFp4_set(&T,Q);
    
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    int i;
    int r_bit;
    r_bit= (int)mpz_sizeinbase(loop,2);
    
    for(i=r_bit-2;i>=0;i--){
        if(mpz_tstbit(loop,i)==1){
            Fp16_mul(&l_sum,&l_sum,&l_sum);
            sparse_dbl_line_kss16(&ltt,&T,&T,P,&Px_neg);
            sparse_add_line_kss16(&ltp,&T,&T,Q,P,&Px_neg);
            Fp16_mul(&l_sum,&l_sum,&ltt);
            Fp16_mul(&l_sum,&l_sum,&ltp);
        }else{
            Fp16_mul(&l_sum,&l_sum,&l_sum);
            sparse_dbl_line_kss16(&ltt,&T,&T,P,&Px_neg);
            Fp16_mul(&l_sum,&l_sum,&ltt);
        }
    }
    Fp16_set(ANS,&l_sum);
    
    Fp16_clear(&l_sum);
    EFp4_clear(&T);
    EFp4_clear(&EFp_tmp);
    Fp4_clear(&Px_neg);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
}
