//
//  EFp2.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_EFp2.h>
// #pragma mark EFp2 methods

void EFp2_init(struct EFp2 *A){
    Fp2_init(&A->x);
    Fp2_init(&A->y);
    A->PoI=FALSE;
}
void EFp2_set(struct EFp2 *A,struct EFp2 *B){
    Fp2_set(&A->x,&B->x);
    Fp2_set(&A->y,&B->y);
    A->PoI=B->PoI;
}
void EFp2_set_PoI(struct EFp2 *A){
    Fp2_set_ui(&A->x,0);
    Fp2_set_ui(&A->y,0);
    A->PoI=TRUE;
}
void EFp2_clear(struct EFp2 *A){
    Fp2_clear(&A->x);
    Fp2_clear(&A->y);
}
void EFp2_printf(struct EFp2 *A){
    gmp_printf("(%Zd,%Zd)\n(%Zd,%Zd)\n",A->x.x0.x0,A->x.x1.x0,A->y.x0.x0,A->y.x1.x0);
}
void EFp2_scm_bin(struct EFp2 *ANS,struct EFp2 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp2 Q,R;
    EFp2_init(&Q);
    EFp2_set(&Q,P);
    EFp2_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp2_ecd(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp2_eca(&Q,&Q,P);
        }
    }
    EFp2_set(ANS,&Q);
    
    EFp2_clear(&Q);
    EFp2_clear(&R);
    return;
}
void EFp2_ecd(struct EFp2 *ANS, struct EFp2 *P){
    if(P->PoI==TRUE){
        EFp2_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp2_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp2_set_PoI(ANS);
        return;
    }
    
    struct Fp2 x,y,lambda,tmp;
    struct EFp2 t_ans;
    Fp2_init(&x);
    Fp2_init(&lambda);
    Fp2_init(&tmp);
    Fp2_init(&y);
    EFp2_init(&t_ans);
    
    Fp2_mul(&x,&P->x,&P->x);
    Fp2_add(&tmp,&x,&x);
    Fp2_add(&x,&tmp,&x);
    Fp_add_mpz(&x.x0,&x.x0,kss_curve_const.a);
    Fp2_add(&y,&P->y,&P->y);
    Fp2_div(&lambda,&x,&y);
    Fp2_mul(&tmp,&lambda,&lambda);
    Fp2_add(&x,&P->x,&P->x);
    Fp2_sub(&x,&tmp,&x);
    Fp2_sub(&tmp,&P->x,&x);
    Fp2_set(&t_ans.x,&x);
    Fp2_mul(&tmp,&tmp,&lambda);
    Fp2_sub(&t_ans.y,&tmp,&P->y);
    
    EFp2_set(ANS,&t_ans);
    
    Fp2_clear(&x);
    Fp2_clear(&lambda);
    Fp2_clear(&y);
    Fp2_clear(&tmp);
    EFp2_clear(&t_ans);
}
void EFp2_eca(struct EFp2 *ANS, struct EFp2 *P1, struct EFp2 *P2){
    if(P2->PoI==TRUE){//if P2==inf
        EFp2_set(ANS,P1);
        return;
    }
    else if(P1->PoI==TRUE){//if P1==inf
        EFp2_set(ANS,P2);
        return;
    }
    else if(Fp2_cmp(&P1->x,&P2->x)==0&&Fp2_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp2_set_PoI(ANS);
        return;
    }
    else if(EFp2_cmp(P1,P2)==0){ // P=Q
        EFp2_ecd(ANS,P1);
        return;
    }
    
    struct Fp2 x,y,lambda,tmp;
    struct EFp2 t_ans;
    
    Fp2_init(&x);
    Fp2_init(&y);
    Fp2_init(&lambda);
    Fp2_init(&tmp);
    EFp2_init(&t_ans);
    
    Fp2_sub(&x,&P2->x,&P1->x);
    Fp2_sub(&y,&P2->y,&P1->y);
    Fp2_div(&lambda,&y,&x);
    Fp2_mul(&tmp,&lambda,&lambda);
    Fp2_add(&x,&P1->x,&P2->x);
    Fp2_sub(&x,&tmp,&x);
    Fp2_sub(&tmp,&P1->x,&x);
    Fp2_set(&t_ans.x,&x);
    Fp2_mul(&tmp,&tmp,&lambda);
    Fp2_sub(&t_ans.y,&tmp,&P1->y);
    
    EFp2_set(ANS,&t_ans);
    
    Fp2_clear(&x);
    Fp2_clear(&y);
    Fp2_clear(&lambda);
    Fp2_clear(&tmp);
    EFp2_clear(&t_ans);
}
int EFp2_cmp(struct EFp2 *A,struct EFp2 *B){
    if(Fp2_cmp(&A->x,&B->x)==0 && Fp2_cmp(&A->y,&B->y)==0){
        return 0;
    }
    return 1;
}
void EFp2_random_set(struct EFp2 *ANS){
    struct EFp2 P;
    EFp2_init(&P);
    
    struct Fp2 x,a,tmp_fp;
    Fp2_init(&a);
    Fp2_init(&x);
    Fp2_init(&tmp_fp);
    
    mpz_t t2,p2,p22,tmp,r2;
    mpz_t set_3;
    mpz_init(set_3);
    mpz_set_ui(set_3,3);
    
    mpz_init(t2);
    mpz_init(p2);
    mpz_init(p22);
    mpz_init(tmp);
    mpz_init(r2);
    
    mpz_pow_ui(t2,params.trace_t,2);
    mpz_mul_ui(p2,params.prime,2);
    mpz_sub(tmp,t2,p2);
    
    mpz_pow_ui(p22,params.prime,2);
    mpz_add_ui(p22,p22,1);
    mpz_sub(p22,p22,tmp);
    
    do{
        Fp2_random(&x);
        Fp2_pow(&a,&x,set_3);
        Fp2_mul_mpz(&tmp_fp, &x, kss_curve_const.a);
        Fp2_add(&a, &a, &tmp_fp);
    }while(Fp2_legendre(&a)!=1);
    // Fp2_printf(&a);
    Fp2_sqrt(&P.y,&a);
    Fp2_set(&P.x,&x);
    
    mpz_t r12_div_r2;
    mpz_init(r12_div_r2);
    mpz_div(r12_div_r2,p22,params.order_r);
    mpz_div(r12_div_r2,r12_div_r2,params.order_r);
    
    EFp2_scm_bin(ANS,&P,r12_div_r2);

    EFp2_clear(&P);
    Fp2_clear(&a);
    Fp2_clear(&x);
    Fp2_clear(&tmp_fp);
    mpz_clear(tmp);
    mpz_clear(t2);
    mpz_clear(p22);
    mpz_clear(p2);
}

void EFp2_rational_point_check(struct EFp2 *A){
    struct EFp2 SS,QQ,TT;
    EFp2_init(&SS);
    EFp2_init(&QQ);
    EFp2_init(&TT);
    EFp2_set(&QQ,A);
    Fp2_mul(&SS.y,&QQ.y,&QQ.y);//y^2
    Fp2_mul(&SS.x,&QQ.x,&QQ.x);
    Fp2_mul(&SS.x,&SS.x,&QQ.x);//x^3
    
    Fp2_mul_mpz(&TT.x,&QQ.x,kss_curve_const.a); //a*x
    Fp2_add(&SS.x,&SS.x,&TT.x);
    if(Fp2_cmp(&SS.x,&SS.y)!=0){
        printf("\nNot Rational point\n");
        EFp2_printf(&SS);
        EFp2_printf(A);
    }
}
