//
//  EFp8.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_EFp8.h>

void EFp8_init(struct EFp8 *A){
    Fp8_init(&A->x);
    Fp8_init(&A->y);
    A->PoI=FALSE;
}

void EFp8_set(struct EFp8 *A,struct EFp8 *B){
    Fp8_set(&A->x,&B->x);
    Fp8_set(&A->y,&B->y);
    A->PoI=B->PoI;
}

void EFp8_set_poi(struct EFp8 *A){
    Fp8_set_ui(&A->x,0);
    Fp8_set_ui(&A->y,0);
    A->PoI=TRUE;
}
void EFp8_clear(struct EFp8 *A){
    Fp8_clear(&A->x);
    Fp8_clear(&A->y);
}
void EFp8_printf(struct EFp8 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd",A->x.x0.x0.x0.x0,A->x.x0.x0.x1.x0,A->x.x0.x1.x0.x0,A->x.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd)\n",A->x.x1.x0.x0.x0,A->x.x1.x0.x1.x0,A->x.x1.x1.x0.x0,A->x.x1.x1.x1.x0);
    
    gmp_printf("(%Zd,%Zd,%Zd,%Zd",A->y.x0.x0.x0.x0,A->y.x0.x0.x1.x0,A->y.x0.x1.x0.x0,A->y.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd)\n",A->y.x1.x0.x0.x0,A->y.x1.x0.x1.x0,A->y.x1.x1.x0.x0,A->y.x1.x1.x1.x0);
    
}
void EFp8_scm_bin(struct EFp8 *ANS,struct EFp8 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp8 Q,R;
    EFp8_init(&Q);
    EFp8_set(&Q,P);
    EFp8_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp8_ecd(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp8_eca(&Q,&Q,P);
        }
    }
    EFp8_set(ANS,&Q);
    
    EFp8_clear(&Q);
    EFp8_clear(&R);
    return;
}
void EFp8_ecd(struct EFp8 *ANS, struct EFp8 *P){
    if(P->PoI==TRUE){
        EFp8_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp8_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp8_set_poi(ANS);
        return;
    }
    
    struct Fp8 x,y,lambda,tmp;
    struct EFp8 t_ans;
    Fp8_init(&x);
    Fp8_init(&lambda);
    Fp8_init(&tmp);
    Fp8_init(&y);
    EFp8_init(&t_ans);
    
    Fp8_mul(&x,&P->x,&P->x);
    Fp8_add(&tmp,&x,&x);
    Fp8_add(&x,&tmp,&x);
    Fp8_add(&y,&P->y,&P->y);
    Fp8_div(&lambda,&x,&y);
    Fp8_mul(&tmp,&lambda,&lambda);
    Fp8_add(&x,&P->x,&P->x);
    Fp8_sub(&x,&tmp,&x);
    Fp8_sub(&tmp,&P->x,&x);
    Fp8_set(&t_ans.x,&x);
    Fp8_mul(&tmp,&tmp,&lambda);
    Fp8_sub(&t_ans.y,&tmp,&P->y);
    
    EFp8_set(ANS,&t_ans);
    
    Fp8_clear(&x);
    Fp8_clear(&lambda);
    Fp8_clear(&y);
    Fp8_clear(&tmp);
    EFp8_clear(&t_ans);
}
void EFp8_eca(struct EFp8 *ANS, struct EFp8 *P1, struct EFp8 *P2){
    if(P2->PoI==TRUE){//if P2==inf
        EFp8_set(ANS,P1);
        return;
    }
    else if(P1->PoI==TRUE){//if P1==inf
        EFp8_set(ANS,P2);
        return;
    }
    else if(Fp8_cmp(&P1->x,&P2->x)==0&&Fp8_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp8_set_poi(ANS);
        return;
    }
    else if(EFp8_cmp(P1,P2)==0){ // P=Q
        EFp8_ecd(ANS,P1);
        return;
    }
    
    struct Fp8 x,y,lambda,tmp;
    struct EFp8 t_ans;
    
    Fp8_init(&x);
    Fp8_init(&y);
    Fp8_init(&lambda);
    Fp8_init(&tmp);
    EFp8_init(&t_ans);
    
    Fp8_sub(&x,&P2->x,&P1->x);
    Fp8_sub(&y,&P2->y,&P1->y);
    Fp8_div(&lambda,&y,&x);
    Fp8_mul(&tmp,&lambda,&lambda);
    Fp8_add(&x,&P1->x,&P2->x);
    Fp8_sub(&x,&tmp,&x);
    Fp8_sub(&tmp,&P1->x,&x);
    Fp8_set(&t_ans.x,&x);
    Fp8_mul(&tmp,&tmp,&lambda);
    Fp8_sub(&t_ans.y,&tmp,&P1->y);
    
    EFp8_set(ANS,&t_ans);
    
    Fp8_clear(&x);
    Fp8_clear(&y);
    Fp8_clear(&lambda);
    Fp8_clear(&tmp);
    EFp8_clear(&t_ans);
}
int EFp8_cmp(struct EFp8 *A,struct EFp8 *B){
    if(Fp8_cmp(&A->x,&B->x)==0 && Fp8_cmp(&A->y,&B->y)==0){
        return 0;
    }
    return 1;
}
