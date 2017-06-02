//
//  KSS16_PS_Line_Evaluation.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_KSS16_PS_Line_Calc.h>

void ps_add_line_kss16(struct Fp16 *l_ANS,struct EFp4 *T_ANS,struct EFp4 *T,struct EFp4 *P,struct EFp4 *Q,struct Fp4 *L){
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
    struct Fp4 A,B,C,D,E;
    Fp4_init(&A);
    Fp4_init(&B);
    Fp4_init(&C);
    Fp4_init(&D);
    Fp4_init(&E);
    
    
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
    
    Fp_set_ui(&l_tmp.x0.x0.x0.x0,1);
    
    Fp4_mul(&l_tmp.x1.x1,&E,L);
    
    Fp4_neg(&l_tmp.x1.x0,&C);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp4_set(T_ANS,&x3_tmp);
    // if((Fp4_cmp(&T->x,&Q->x))==0&&(Fp4_cmp(&T->y,&Q->y))!=0){//xt==xp&&yt!=yp
    // 	Fp4_sub(&ltp,&P->x,&T->x);
    // 	Fp4_set(&l_ANS->x0.x0,&ltp);
    // }
    // if(T->PoI==TRUE){//if P2==inf
    // 	EFp4_set(T_ANS,P);
    // 	return;
    // }
    // else if(P->PoI==TRUE){//if P1==inf
    // 	EFp4_set(T_ANS,T);
    // 	return;
    // }
    // else if(Fp4_cmp(&T->x,&P->x)==0&&Fp4_cmp(&T->y,&P->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
    // 	EFp4_set_PoI(T_ANS);
    // 	return;
    // }
    // else if(EFp4_cmp(T,P)==0){ // P=P
    // 	EFp4_ecd(T_ANS,T);
    // 	return;
    // }
    
    
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
}
void ps_dbl_line_kss16(struct Fp16 *l_ANS,struct EFp4 *T_ANS,struct EFp4 *T,struct EFp4 *Q,struct Fp4 *L){
    struct Fp4 tmp1,tmp2,tmp3,tmp4,lambda,ltt;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    Fp4_init(&tmp4);
    Fp4_init(&lambda);
    Fp4_init(&ltt);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp4 x,y,tmp;
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&tmp);
    
    struct EFp4 x3_tmp;
    EFp4_init(&x3_tmp);
    struct Fp4 A,B,C,D,E;
    Fp4_init(&A);
    Fp4_init(&B);
    Fp4_init(&C);
    Fp4_init(&D);
    Fp4_init(&E);
    
    mpz_t ac_inv;
    mpz_init(ac_inv);
    
    Fp4_add(&A,&T->y,&T->y);//2yt
    Fp4_mul(&B,&T->x,&T->x);//xt^2
    Fp4_mul_ui(&B,&B,3);//3xt^2
    Fp4_mul_beta_inv(&tmp);
    Fp4_mul(&tmp, &tmp, &z_inv2);
    Fp4_add(&B, &B, &tmp);
    
    Fp4_div(&C,&B,&A);//
    
    Fp4_add(&D,&T->x,&T->x);//D=2xt
    Fp4_mul(&tmp1,&C,&C);//C^2
    Fp4_sub(&x3_tmp.x,&tmp1,&D);//x3.x = C^2-D
    
    Fp4_mul(&tmp2,&C,&T->x); // C*xt
    Fp4_sub(&E,&tmp2,&T->y); //E=C*xt-yt
    
    Fp4_mul(&tmp3,&C,&x3_tmp.x);//C*x3.x
    Fp4_sub(&x3_tmp.y,&E,&tmp3);//x3.y = E-C*x3.x
    
    Fp_set_ui(&l_tmp.x0.x0.x0.x0,1);
    
    Fp4_mul(&l_tmp.x1.x1,&E,L);//l_tmp=
    
    Fp4_neg(&l_tmp.x1.x0,&C);
    Fp16_set(l_ANS,&l_tmp);
    
    EFp4_set(T_ANS,&x3_tmp);
    
    // if(T->PoI==TRUE){
    // 	EFp4_set(T_ANS,T);
    // 	return;
    // }
    // mpz_t cmp;
    // mpz_init(cmp);
    // mpz_set_ui(cmp,0);
    // if(Fp4_cmp_mpz(&T->y,cmp)==0){//P.y==0
    // 	EFp4_set_PoI(T_ANS);
    // 	return;
    // }
    
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp4_clear(&tmp4);
    Fp4_clear(&lambda);
    Fp4_clear(&ltt);
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
}
void ps_mul_line_kss16(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B){
    
    struct Fp16 a,b;
    Fp16_init(&a);
    Fp16_init(&b);
    
    struct Fp4 a0,a1,a2,a3,b2,b3,C_0,C_1,C_2,C_3,T0,T1,T2,T3,T4,a2a3,tmp2,tmp3;
    Fp4_init(&a0);
    Fp4_init(&a1);
    Fp4_init(&a2);
    Fp4_init(&a3);
    Fp4_init(&b2);
    Fp4_init(&b3);
    
    Fp4_init(&C_0);
    Fp4_init(&C_1);
    Fp4_init(&C_2);
    Fp4_init(&C_3);
    
    Fp4_init(&T0);
    Fp4_init(&T1);
    Fp4_init(&T2);
    Fp4_init(&T3);
    Fp4_init(&T4);
    
    Fp4_init(&a2a3);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    
    
    Fp16_set(&a, A);
    Fp16_set(&b, B);
    
    Fp4_set(&a0, &a.x0.x0);
    Fp4_set(&a1, &a.x0.x1);
    Fp4_set(&a2, &a.x1.x0);
    Fp4_set(&a3, &a.x1.x1);
    
    Fp4_set(&b2, &b.x1.x0);
    Fp4_set(&b3, &b.x1.x1);
    
    
    Fp4_add(&a2a3, &a2, &a3);   /**< (a2+a3) */
    Fp4_add(&T4, &b2, &b3);     /**< t4=(b2+b3) */
    
    Fp4_mul(&T1, &a2, &b2);     /**< t1=(a2*b2) */
    Fp4_mul(&T0, &a3, &b3);     /**< t0=(a3*b3) */
    //    Fp4_mul_v(&T0, &tmp2);/**< t0=(a3*b3)*beta */
    
    Fp4_mul(&tmp2, &a2a3, &T4);
    Fp4_sub(&tmp3, &tmp2, &T1);
    Fp4_sub(&C_0, &tmp3, &T0);
    Fp4_mul_v(&tmp3, &C_0);
    Fp4_add(&C_0, &tmp3, &a0);
    
    Fp4_mul_v(&tmp2, &T0);
    Fp4_add(&C_1, &T1, &tmp2);
    Fp4_add(&C_1, &C_1, &a1);
    
    Fp4_mul(&T3, &a0, &b2);
    Fp4_mul(&T2, &a1, &b3);
    Fp4_mul_v(&tmp2, &T2);
    Fp4_add(&C_2, &T3, &tmp2);
    Fp4_add(&C_2, &C_2, &a2);
    
    
    Fp4_add(&tmp2, &a0, &a1);
    Fp4_mul(&tmp3, &tmp2, &T4);
    Fp4_sub(&tmp2, &tmp3, &T3);
    Fp4_sub(&C_3, &tmp2, &T2);
    Fp4_add(&C_3, &C_3, &a3);
    
    Fp4_set(&a.x0.x0,&C_0);
    Fp4_set(&a.x0.x1,&C_1);
    Fp4_set(&a.x1.x0,&C_2);
    Fp4_set(&a.x1.x1,&C_3);
    
    Fp16_set(ANS, &a);
    
    
    Fp16_clear(&a);
    Fp16_clear(&b);
    Fp4_clear(&a0);
    Fp4_clear(&a1);
    Fp4_clear(&a2);
    Fp4_clear(&a3);
    Fp4_clear(&b2);
    Fp4_clear(&b3);
    Fp4_clear(&C_0);
    Fp4_clear(&C_1);
    Fp4_clear(&C_2);
    Fp4_clear(&C_3);
    Fp4_clear(&T0);
    Fp4_clear(&T1);
    Fp4_clear(&T2);
    Fp4_clear(&T3);
    Fp4_clear(&T4);
}
