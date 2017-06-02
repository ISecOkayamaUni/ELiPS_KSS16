//
//  KSS16_Sparse7_Line_Evaluation.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_KSS16_Sparse_Line_Calc.h>

void sparse_add_line_kss16(struct Fp16 *l_ANS,struct EFp4 *T_ANS,struct EFp4 *T,struct EFp4 *P,struct EFp4 *Q,struct Fp4 *Qx_neg){
    struct Fp4 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    Fp4_init(&tmp4);
    Fp4_init(&lambda);
    Fp4_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp4 x,y,tmp;
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&tmp);
    
    struct EFp4 x3_tmp;
    EFp4_init(&x3_tmp);
    struct Fp4 A,B,C,D,E,F;
    Fp4_init(&A);
    Fp4_init(&B);
    Fp4_init(&C);
    Fp4_init(&D);
    Fp4_init(&E);
    Fp4_init(&F);
    
    Fp4_sub(&A,&P->x,&T->x);//xt-xp
    Fp4_sub(&B,&P->y,&T->y);//yt-yp
    Fp4_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)
    
    
    Fp4_add(&D,&T->x,&P->x);
    Fp4_mul(&tmp1,&C,&C);
    Fp4_sub(&x3_tmp.x,&tmp1,&D);
    
    Fp4_mul(&tmp2,&C,&T->x);
    Fp4_sub(&E,&tmp2,&T->y);
    
    Fp4_mul(&tmp3,&C,&x3_tmp.x);
    Fp4_sub(&x3_tmp.y,&E,&tmp3);
    
    Fp4_set(&l_tmp.x0.x0,&Q->y);
    
    Fp4_set(&l_tmp.x1.x1,&E);
    
    Fp4_mul(&F,&C,Qx_neg);
    Fp4_set(&l_tmp.x1.x0,&F);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp4_set(T_ANS,&x3_tmp);
    
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp4_clear(&tmp4);
    Fp4_clear(&lambda);
    Fp4_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&x3_tmp);
    Fp4_clear(&A);
    Fp4_clear(&B);
    Fp4_clear(&C);
    Fp4_clear(&D);
    Fp4_clear(&E);
    Fp4_clear(&F);
}

void sparse_dbl_line_kss16(struct Fp16 *l_ANS,struct EFp4 *T_ANS,struct EFp4 *T,struct EFp4 *Q,struct Fp4 *Qx_neg){
    struct Fp4 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    Fp4_init(&tmp4);
    Fp4_init(&lambda);
    Fp4_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp4 x,y,tmp;
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&tmp);
    
    struct EFp4 x3_tmp;
    EFp4_init(&x3_tmp);
    struct Fp4 A,B,C,D,E,F;
    Fp4_init(&A);
    Fp4_init(&B);
    Fp4_init(&C);
    Fp4_init(&D);
    Fp4_init(&E);
    Fp4_init(&F);
    
    
    Fp4_add(&A,&T->y,&T->y);//xt-xp
    Fp4_mul(&B,&T->x,&T->x);
    Fp4_mul_ui(&B,&B,3);
    struct Fp4 ac_inv;
    Fp4_init(&ac_inv);
    Fp4_mul_beta_inv(&ac_inv);
    Fp4_add(&B,&B,&ac_inv);
    Fp4_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)
    
    Fp4_add(&D,&T->x,&T->x);
    Fp4_mul(&tmp1,&C,&C);
    Fp4_sub(&x3_tmp.x,&tmp1,&D);
    
    Fp4_mul(&tmp2,&C,&T->x);
    Fp4_sub(&E,&tmp2,&T->y);
    
    Fp4_mul(&tmp3,&C,&x3_tmp.x);
    Fp4_sub(&x3_tmp.y,&E,&tmp3);
    
    Fp4_set(&l_tmp.x0.x0,&Q->y);
    
    Fp4_set(&l_tmp.x1.x1,&E);
    
    Fp4_mul(&F,&C,Qx_neg);
    Fp4_set(&l_tmp.x1.x0,&F);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp4_set(T_ANS,&x3_tmp);
    
    if(T->PoI==TRUE){
        EFp4_set(T_ANS,T);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp4_cmp_mpz(&T->y,cmp)==0){//P.y==0
        EFp4_set_PoI(T_ANS);
        return;
    }
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp4_clear(&tmp4);
    Fp4_clear(&lambda);
    Fp4_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&x3_tmp);
    Fp4_clear(&A);
    Fp4_clear(&B);
    Fp4_clear(&C);
    Fp4_clear(&D);
    Fp4_clear(&E);
    Fp4_clear(&F);
    mpz_clear(cmp);
}
