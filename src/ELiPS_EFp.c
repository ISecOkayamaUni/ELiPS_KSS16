//
//  EFp.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/26/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_EFp.h>

//--------------------------
// #pragma mark EFp methods
void EFp_init(struct EFp *A) {
    Fp_init (&A->x);
    Fp_init (&A->y);
    A->PoI = FALSE;
}
void EFp_set(struct EFp *A,struct EFp *B) {
    Fp_set (&A->x, &B->x);
    Fp_set (&A->y, &B->y);
    A->PoI = B->PoI;
}
void EFp_set_PoI(struct EFp *A) {
    Fp_set_ui(&A->x,0);
    Fp_set_ui(&A->y,0);
    A->PoI=TRUE;
}
void EFp_clear (struct EFp *A) {
    Fp_clear (&A->x);
    Fp_clear (&A->y);
}
void EFp_printf(struct EFp *A) {
    gmp_printf ("(%Zd, %Zd)\n", A->x.x0, A->y.x0);
}
void EFp_scm_bin(struct EFp *ANS, struct EFp *P,mpz_t j){
    int i;
    int r;//bit数
    r = (int)mpz_sizeinbase (j, 2);
    
    struct EFp Q;
    EFp_init(&Q);
    EFp_set(&Q,P);
    
    for (i = r-2; i >= 0; i--) {
        if (mpz_tstbit(j, i) == 1) {
            EFp_ecd (&Q,&Q);
            EFp_eca (&Q,&Q,P);
        }
        else {
            EFp_ecd (&Q,&Q);
        }
    }
    
    EFp_set (ANS,&Q);
    EFp_clear (&Q);
    return;
}


void EFp_ecd(struct EFp *ANS, struct EFp *P){
    if(P->PoI==TRUE){
        EFp_set(ANS,P);
        return;
    }
    if(mpz_sgn(P->y.x0)==0){//P.y==0
        EFp_set_PoI(ANS);
        return;
    }
    
    struct Fp x,y,lambda,tmp;
    struct EFp t_ans;
    Fp_init(&x);
    Fp_init(&lambda);
    Fp_init(&tmp);
    Fp_init(&y);
    EFp_init(&t_ans);
    
    Fp_mul(&x,&P->x,&P->x);
    Fp_add(&tmp,&x,&x);
    Fp_add(&x,&tmp,&x);//3x^2+a
    Fp_add_mpz(&x,&x,kss_curve_const.tmp_a);
    Fp_add(&y,&P->y,&P->y);//2y
    
    Fp_div(&lambda,&x,&y);
    Fp_mul(&tmp,&lambda,&lambda);
    Fp_add(&x,&P->x,&P->x);
    Fp_sub(&x,&tmp,&x);
    Fp_sub(&tmp,&P->x,&x);
    
    
    Fp_set(&t_ans.x,&x);
    
    Fp_mul(&tmp,&tmp,&lambda);
    Fp_sub(&t_ans.y,&tmp,&P->y);
    
    EFp_set(ANS,&t_ans);
    
    Fp_clear(&x);
    Fp_clear(&lambda);
    Fp_clear(&y);
    Fp_clear(&tmp);
    EFp_clear(&t_ans);
}
void EFp_eca(struct EFp *ANS, struct EFp *P1, struct EFp *P2){
    if(P2->PoI==TRUE){//if P2==inf
        EFp_set(ANS,P1);
        return;
    }
    else if(P1->PoI==TRUE){//if P1==inf
        EFp_set(ANS,P2);
        return;
    }
    else if(Fp_cmp(&P1->x,&P2->x)==0&&Fp_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp_set_PoI(ANS);
        return;
    }
    else if(EFp_cmp(P1,P2)==0){ // P=Q
        EFp_ecd(ANS,P1);
        return;
    }
    
    struct Fp x,y,lambda,tmp;
    struct EFp t_ans;
    
    Fp_init(&x);
    Fp_init(&y);
    Fp_init(&lambda);
    Fp_init(&tmp);
    EFp_init(&t_ans);
    
    Fp_sub(&x,&P2->x,&P1->x);
    Fp_sub(&y,&P2->y,&P1->y);
    Fp_div(&lambda,&y,&x);
    Fp_mul(&tmp,&lambda,&lambda);
    Fp_add(&x,&P1->x,&P2->x);
    Fp_sub(&x,&tmp,&x);
    Fp_sub(&tmp,&P1->x,&x);
    Fp_set(&t_ans.x,&x);
    Fp_mul(&tmp,&tmp,&lambda);
    Fp_sub(&t_ans.y,&tmp,&P1->y);
    
    EFp_set(ANS,&t_ans);
    
    Fp_clear(&x);
    Fp_clear(&y);
    Fp_clear(&lambda);
    Fp_clear(&tmp);
    EFp_clear(&t_ans);
}
int EFp_cmp(struct EFp *A,struct EFp *B){
    if(Fp_cmp(&A->x,&B->x)==0 && Fp_cmp(&A->y,&B->y)==0){
        return 0;
    }
    return 1;
}
void EFp_neg(struct EFp *ANS, struct EFp *A){
    struct EFp tmp;
    EFp_init(&tmp);
    Fp_neg(&tmp.y,&A->y);
    Fp_set(&tmp.x,&A->x);
    EFp_set(ANS,&tmp);
    EFp_clear(&tmp);
}
void EFp_random_set(struct EFp *ANS){
    struct Fp a,x,tmp;
    Fp_init(&a);
    Fp_init(&x);
    Fp_init(&tmp);
    
    struct EFp P,Q;
    EFp_init(&P);
    EFp_init(&Q);
    
    
    do{
        Fp_random(&x);
        Fp_mul(&a,&x,&x);
        Fp_mul(&a,&a,&x);
        Fp_mul_mpz(&tmp, &x, kss_curve_const.a);
        Fp_add(&a, &a, &tmp);
    }while(mpz_legendre(a.x0,params.prime)!=1);
    Q.PoI=0;
    Fp_sqrt(&P.y,&a);
    Fp_set(&P.x,&x);
    EFp_set(ANS,&P);
    rational_point_check(ANS);
    
    Fp_clear(&a);
    Fp_clear(&x);
    EFp_clear(&P);
}

void rational_point_check(struct EFp *A){
    struct EFp SS,QQ,TT;
    EFp_init(&SS);
    EFp_init(&QQ);
    EFp_init(&TT);
    EFp_set(&QQ,A);
    Fp_mul(&SS.y,&QQ.y,&QQ.y);//y^2
    Fp_mul(&SS.x,&QQ.x,&QQ.x);
    Fp_mul(&SS.x,&SS.x,&QQ.x);//x^3
    
    Fp_mul_mpz(&TT.x,&QQ.x,kss_curve_const.a); //a*x
    Fp_add(&SS.x,&SS.x,&TT.x);
    if(Fp_cmp(&SS.x,&SS.y)!=0){
        printf("\nNot Rational point\n");
        EFp_printf(&SS);
        EFp_printf(A);
    }
}


