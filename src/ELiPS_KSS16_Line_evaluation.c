//
//  Line_evaluation.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_KSS16_Line_Evaluation.h>

void add_line_kss16(struct Fp16 *l_ANS,struct EFp16 *T_ANS,struct EFp16 *T,struct EFp16 *P,struct EFp16 *Q,struct Fp16 *Qx_neg){
    struct Fp16 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp16_init(&tmp1);
    Fp16_init(&tmp2);
    Fp16_init(&tmp3);
    Fp16_init(&tmp4);
    Fp16_init(&lambda);
    Fp16_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp16 x,y,tmp;
    Fp16_init(&x);
    Fp16_init(&y);
    Fp16_init(&tmp);
    
    struct EFp16 x3_tmp;
    EFp16_init(&x3_tmp);
    struct Fp16 A,B,C,D,E,F;
    Fp16_init(&A);
    Fp16_init(&B);
    Fp16_init(&C);
    Fp16_init(&D);
    Fp16_init(&E);
    Fp16_init(&F);
    
    
    Fp16_sub(&A,&P->x,&T->x);//xt-xp
    Fp16_sub(&B,&P->y,&T->y);//yt-yp
    Fp16_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)
    
    Fp16_add(&D,&T->x,&P->x);
    Fp16_mul(&tmp1,&C,&C);
    Fp16_sub(&x3_tmp.x,&tmp1,&D);
    
    Fp16_mul(&tmp2,&C,&T->x);
    Fp16_sub(&E,&tmp2,&T->y);
    
    Fp16_mul(&tmp3,&C,&x3_tmp.x);
    Fp16_sub(&x3_tmp.y,&E,&tmp3);
    
    Fp16_set(&l_tmp,&Q->y);
    
    Fp16_add(&l_tmp,&l_tmp,&E);
    
    Fp16_mul(&F,&C,Qx_neg);
    Fp16_add(&l_tmp,&l_tmp,&F);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp16_set(T_ANS,&x3_tmp);
    
    Fp16_clear(&tmp1);
    Fp16_clear(&tmp2);
    Fp16_clear(&tmp3);
    Fp16_clear(&tmp4);
    Fp16_clear(&lambda);
    Fp16_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp16_clear(&x);
    Fp16_clear(&y);
    Fp16_clear(&tmp);
    EFp16_clear(&x3_tmp);
    Fp16_clear(&A);
    Fp16_clear(&B);
    Fp16_clear(&C);
    Fp16_clear(&D);
    Fp16_clear(&E);
    Fp16_clear(&F);
}
//TODO *Q = *P, *Qx_neg = * Px_neg
void dbl_line_kss16(struct Fp16 *l_ANS,struct EFp16 *T_ANS,struct EFp16 *T,struct EFp16 *P,struct Fp16 *Px_neg){
    struct Fp16 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp16_init(&tmp1);
    Fp16_init(&tmp2);
    Fp16_init(&tmp3);
    Fp16_init(&tmp4);
    Fp16_init(&lambda);
    Fp16_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp16 x,y,tmp;
    Fp16_init(&x);
    Fp16_init(&y);
    Fp16_init(&tmp);
    
    struct EFp16 x3_tmp;
    EFp16_init(&x3_tmp);
    struct Fp16 A,B,C,D,E,F;
    Fp16_init(&A);
    Fp16_init(&B);
    Fp16_init(&C);
    Fp16_init(&D);
    Fp16_init(&E);
    Fp16_init(&F);
    
    
    Fp16_add(&A,&T->y,&T->y);//2y
    Fp16_mul(&B,&T->x,&T->x);//x^2
    Fp16_mul_ui(&B,&B,3);//3x^2
    Fp_add_mpz(&B.x0.x0.x0.x0,&B.x0.x0.x0.x0,kss_curve_const.a);//lambda=3x^2+a
    Fp16_div(&C,&B,&A);//lambda=3x^2+a/2y
    
    Fp16_add(&D,&T->x,&T->x); //D=2x
    Fp16_mul(&tmp1,&C,&C);// lamda^2
    Fp16_sub(&x3_tmp.x,&tmp1,&D); //x3.x=lamda^2-x2t
    
    Fp16_mul(&tmp2,&C,&T->x);//xt.lamda
    Fp16_sub(&E,&tmp2,&T->y);//xt.lamda-yt
    
    Fp16_mul(&tmp3,&C,&x3_tmp.x); //x3.lamda
    Fp16_sub(&x3_tmp.y,&E,&tmp3); // x3.y = xt.lamda-yt - x3.lamda
    
    Fp16_set(&l_tmp,&P->y);
    // Fp_set_ui(&l_tmp.x0.x0.x0,1);
    
    Fp16_add(&l_tmp,&l_tmp,&E);
    
    Fp16_mul(&F,&C,Px_neg);
    Fp16_add(&l_tmp,&l_tmp,&F);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp16_set(T_ANS,&x3_tmp);
    
    if(T->PoI==TRUE){
        EFp16_set(T_ANS,T);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp16_cmp_mpz(&T->y,cmp)==0){//P.y==0
        EFp16_set_poi(T_ANS);
        return;
    }
    Fp16_clear(&tmp1);
    Fp16_clear(&tmp2);
    Fp16_clear(&tmp3);
    Fp16_clear(&tmp4);
    Fp16_clear(&lambda);
    Fp16_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp16_clear(&x);
    Fp16_clear(&y);
    Fp16_clear(&tmp);
    EFp16_clear(&x3_tmp);
    Fp16_clear(&A);
    Fp16_clear(&B);
    Fp16_clear(&C);
    Fp16_clear(&D);
    Fp16_clear(&E);
    Fp16_clear(&F);
    mpz_clear(cmp);
}
//-------------------
void ltt_q_kss16(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q){
    struct Fp16 tmp1,tmp2,tmp3,lambda,ltt;
    Fp16_init(&tmp1);
    Fp16_init(&tmp2);
    Fp16_init(&tmp3);
    Fp16_init(&lambda);
    Fp16_init(&ltt);
    
    Fp16_mul(&tmp1,&T->x,&T->x);//xt^2
    Fp16_add(&tmp2,&tmp1,&tmp1);
    Fp16_add(&tmp1,&tmp2,&tmp1);//3xt^2
    Fp_add_mpz(&tmp1.x0.x0.x0.x0,&tmp1.x0.x0.x0.x0,kss_curve_const.a);//TODO
    Fp16_add(&tmp2,&T->y,&T->y);//2yt
    
    Fp16_div(&lambda,&tmp1,&tmp2);//lambda=3xt^2+a/2yt
    Fp16_sub(&tmp3,&Q->x,&T->x);//tmp3=xq-xt
    Fp16_mul(&tmp3,&tmp3,&lambda);//tmp3=lambda(xq-xt)
    
    Fp16_sub(&ltt,&Q->y,&T->y);//yq-yt
    Fp16_sub(&ltt,&ltt,&tmp3);//ltt=yq-yt-lambda(xq-xt)
    
    Fp16_set(ANS,&ltt);
    
    Fp16_clear(&tmp1);
    Fp16_clear(&tmp2);
    Fp16_clear(&tmp3);
    Fp16_clear(&lambda);
    Fp16_clear(&ltt);
}
void v2t_q_kss16(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q){
    struct Fp16 v2t;
    Fp16_init(&v2t);
    
    Fp16_sub(&v2t,&Q->x,&T->x);//v2t=xq-xt
    Fp16_set(ANS,&v2t);
    
    Fp16_clear(&v2t);
}
void ltp_q_kss16(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *P,struct EFp16 *Q){
    struct Fp16 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp16_init(&tmp1);
    Fp16_init(&tmp2);
    Fp16_init(&tmp3);
    Fp16_init(&tmp4);
    Fp16_init(&lambda);
    Fp16_init(&ltp);
    
    if((Fp16_cmp(&T->x,&P->x))==0&&(Fp16_cmp(&T->y,&P->y))!=0){//xt==xp&&yt!=yp
        Fp16_sub(&ltp,&Q->x,&T->x);
        Fp16_set(ANS,&ltp);
        
        return;
    }
    
    Fp16_sub(&tmp1,&T->x,&P->x);//xt-xp
    Fp16_sub(&tmp2,&T->y,&P->y);//yt-yp
    Fp16_div(&lambda,&tmp2,&tmp1);//lambda=(yt-tp)/(xt-xp)
    
    Fp16_sub(&tmp3,&Q->x,&T->x);//tmp3=(xq-xt)
    Fp16_mul(&tmp4,&tmp3,&lambda);//tmp4=lambda(xq-xt)
    
    Fp16_sub(&ltp,&Q->y,&T->y);//ltp=yq-yt
    Fp16_sub(&ltp,&ltp,&tmp4);//ltp=yq-yt-lambda(xq-xt)
    
    Fp16_set(ANS,&ltp);
    
    Fp16_clear(&tmp1);
    Fp16_clear(&tmp2);
    Fp16_clear(&tmp3);
    Fp16_clear(&tmp4);
    Fp16_clear(&lambda);
    Fp16_clear(&ltp);
}
void vtp_q_kss16(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q){
    struct Fp16 vtp;
    Fp16_init(&vtp);
    if(T->PoI==1){//if T is PoI
        Fp16_set_ui(ANS,0);
        Fp_set_ui(&ANS->x0.x0.x0.x0,1);
        return;
    }
    
    Fp16_sub(&vtp,&Q->x,&T->x);
    Fp16_set(ANS,&vtp);
    
    Fp16_clear(&vtp);
}
