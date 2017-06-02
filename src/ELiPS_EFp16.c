//
//  EFp16.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_EFp16.h>

#pragma mark EFp16 methods

void EFp16_init(struct EFp16 *A){
    Fp16_init(&A->x);
    Fp16_init(&A->y);
    A->PoI=FALSE;
}
void EFp16_set(struct EFp16 *A,struct EFp16 *B){
    Fp16_set(&A->x,&B->x);
    Fp16_set(&A->y,&B->y);
    A->PoI=B->PoI;
}
void EFp16_set_poi(struct EFp16 *A){
    Fp16_set_ui(&A->x,0);
    Fp16_set_ui(&A->y,0);
    A->PoI=TRUE;
}
void EFp16_set_EFp(struct EFp16 *A,struct EFp *B){
    Fp16_set_ui(&A->x,0);
    Fp16_set_ui(&A->y,0);
    
    Fp_set(&A->x.x0.x0.x0.x0,&B->x);
    Fp_set(&A->y.x0.x0.x0.x0,&B->y);
    A->PoI=B->PoI;
}
void EFp16_clear(struct EFp16 *A){
    Fp16_clear(&A->x);
    Fp16_clear(&A->y);
}
void EFp16_printf(struct EFp16 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd,\n",A->x.x0.x0.x0.x0.x0,A->x.x0.x0.x0.x1.x0,A->x.x0.x0.x1.x0.x0,A->x.x0.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->x.x0.x1.x0.x0.x0,A->x.x0.x1.x0.x1.x0,A->x.x0.x1.x1.x0.x0,A->x.x0.x1.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->x.x1.x0.x0.x0.x0,A->x.x1.x0.x0.x1.x0,A->x.x1.x0.x1.x0.x0,A->x.x1.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd)\n",A->x.x1.x1.x0.x0.x0,A->x.x1.x1.x0.x1.x0,A->x.x1.x1.x1.x0.x0,A->x.x1.x1.x1.x1.x0);
    
    gmp_printf("(%Zd,%Zd,%Zd,%Zd,\n",A->y.x0.x0.x0.x0.x0,A->y.x0.x0.x0.x1.x0,A->y.x0.x0.x1.x0.x0,A->y.x0.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->y.x0.x1.x0.x0.x0,A->y.x0.x1.x0.x1.x0,A->y.x0.x1.x1.x0.x0,A->y.x0.x1.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->y.x1.x0.x0.x0.x0,A->y.x1.x0.x0.x1.x0,A->y.x1.x0.x1.x0.x0,A->y.x1.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd)\n",A->y.x1.x1.x0.x0.x0,A->y.x1.x1.x0.x1.x0,A->y.x1.x1.x1.x0.x0,A->y.x1.x1.x1.x1.x0);
}
void EFp16_scm_bin(struct EFp16 *ANS,struct EFp16 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp16 Q,R;
    EFp16_init(&Q);
    EFp16_set(&Q,P);
    EFp16_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp16_ecd(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp16_eca(&Q,&Q,P);
        }
    }
    EFp16_set(ANS,&Q);
    
    EFp16_clear(&Q);
    EFp16_clear(&R);
    return;
}
void EFp16_ecd(struct EFp16 *ANS, struct EFp16 *P){
    if(P->PoI==TRUE){
        EFp16_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp16_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp16_set_poi(ANS);
        return;
    }
    
    struct Fp16 x,y,lambda,tmp;
    struct EFp16 t_ans;
    Fp16_init(&x);
    Fp16_init(&lambda);
    Fp16_init(&tmp);
    Fp16_init(&y);
    EFp16_init(&t_ans);
    
    Fp16_mul(&x,&P->x,&P->x);
    Fp16_add(&tmp,&x,&x);
    Fp16_add(&x,&tmp,&x);
    Fp_add_mpz(&x.x0.x0.x0.x0,&x.x0.x0.x0.x0,kss_curve_const.a);
    Fp16_add(&y,&P->y,&P->y);
    Fp16_div(&lambda,&x,&y);
    Fp16_mul(&tmp,&lambda,&lambda);
    Fp16_add(&x,&P->x,&P->x);
    Fp16_sub(&x,&tmp,&x);
    Fp16_sub(&tmp,&P->x,&x);
    Fp16_set(&t_ans.x,&x);
    Fp16_mul(&tmp,&tmp,&lambda);
    Fp16_sub(&t_ans.y,&tmp,&P->y);
    
    EFp16_set(ANS,&t_ans);
    
    Fp16_clear(&x);
    Fp16_clear(&lambda);
    Fp16_clear(&y);
    Fp16_clear(&tmp);
    EFp16_clear(&t_ans);
}
void EFp16_eca(struct EFp16 *ANS, struct EFp16 *P1, struct EFp16 *P2){
    if(P2->PoI==TRUE){//if P2==inf
        EFp16_set(ANS,P1);
        return;
    }
    else if(P1->PoI==TRUE){//if P1==inf
        EFp16_set(ANS,P2);
        return;
    }
    else if(Fp16_cmp(&P1->x,&P2->x)==0&&Fp16_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp16_set_poi(ANS);
        return;
    }
    else if(EFp16_cmp(P1,P2)==0){ // P=Q
        EFp16_ecd(ANS,P1);
        return;
    }
    
    struct Fp16 x,y,lambda,tmp;
    struct EFp16 t_ans;
    
    Fp16_init(&x);
    Fp16_init(&y);
    Fp16_init(&lambda);
    Fp16_init(&tmp);
    EFp16_init(&t_ans);
    
    Fp16_sub(&x,&P2->x,&P1->x);
    Fp16_sub(&y,&P2->y,&P1->y);
    Fp16_div(&lambda,&y,&x);
    Fp16_mul(&tmp,&lambda,&lambda);
    Fp16_add(&x,&P1->x,&P2->x);
    Fp16_sub(&x,&tmp,&x);
    Fp16_sub(&tmp,&P1->x,&x);
    Fp16_set(&t_ans.x,&x);
    Fp16_mul(&tmp,&tmp,&lambda);
    Fp16_sub(&t_ans.y,&tmp,&P1->y);
    
    EFp16_set(ANS,&t_ans);
    
    Fp16_clear(&x);
    Fp16_clear(&y);
    Fp16_clear(&lambda);
    Fp16_clear(&tmp);
    EFp16_clear(&t_ans);
}
int EFp16_cmp(struct EFp16 *A,struct EFp16 *B){
    if(Fp16_cmp(&A->x,&B->x)==0 && Fp16_cmp(&A->y,&B->y)==0){
        return 0;
    }
    return 1;
}
void EFp16_random_set(struct EFp16 *ANS){
    struct EFp16 P,ans_temp;
    EFp16_init(&P);
    EFp16_init(&ans_temp);
    EFp16_init(&P);
    
    struct Fp16 x,a,tmp16;
    Fp16_init(&a);
    Fp16_init(&x);
    Fp16_init(&tmp16);
    
    //t16=a^16+b^16=((t^2-2p)^2-2p^2)^2-2p^4)^2-2p^8
    mpz_t t2,p2,p22,p4,p8,tmp1,tmp2,t16;
    mpz_init(t2);
    mpz_init(p2);
    mpz_init(p22);
    mpz_init(p4);
    mpz_init(p8);
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init(t16);
    
    mpz_pow_ui(tmp1,params.trace_t,2);//t^2
    mpz_mul_ui(p2,params.prime,2);//2p
    mpz_sub(t2,tmp1,p2); //t2=t^2-2p
    mpz_pow_ui(t2,t2,2);//t2=(t^2-2p)^2
    
    mpz_pow_ui(tmp1,params.prime,2); //p^2
    mpz_mul_ui(p22,tmp1,2);//2p^2
    mpz_sub(tmp1,t2,p22); // (t^2-2p)^2-2p^2
    mpz_pow_ui(tmp2,tmp1,2);//tmp2=((t^2-2p)^2-2p^2)^2
    
    mpz_pow_ui(tmp1,params.prime,4); //p^4
    mpz_mul_ui(p4,tmp1,2);//2p^4
    mpz_sub(tmp1,tmp2,p4); //(((t^2-2p)^2-2p^2)^2-2p^4)
    mpz_pow_ui(tmp2,tmp1,2);//tmp2=(((t^2-2p)^2-2p^2)^2-2p^4)^2
    
    mpz_pow_ui(tmp1,params.prime,8); //p^8
    mpz_mul_ui(p8,tmp1,2);//2p^8
    mpz_sub(t16,tmp2,p8);
    
    
    mpz_t r16_div_r2,sEFp_16;
    mpz_init(r16_div_r2);
    mpz_init(sEFp_16);
    mpz_pow_ui(tmp1,params.prime,16);
    mpz_add_ui(tmp1,tmp1,1);
    mpz_sub(sEFp_16,tmp1,t16);
    
    mpz_pow_ui(tmp1,params.order_r,2);
    
    //    printf("r^2 divisible %d\n",(int)mpz_divisible_p(sEFp_16,tmp1));
    mpz_tdiv_q(r16_div_r2,sEFp_16,params.order_r);
    mpz_tdiv_q(r16_div_r2,r16_div_r2,params.order_r);
    do{
        Fp16_random(&x);
        Fp16_mul(&a,&x,&x);
        Fp16_mul(&a,&a,&x);//x^3
        Fp16_mul_mpz(&tmp16,&x, kss_curve_const.a); //ax
        Fp16_add(&a, &a, &tmp16);//x^3+ax
    }while(Fp16_legendre(&a)!=1);
    Fp16_sqrt(&P.y,&a);
    Fp16_set(&P.x,&x);
    EFp16_scm_bin(ANS,&P,r16_div_r2);//R
    
    
    EFp16_scm_bin(&ans_temp, ANS, params.order_r);//T
    if (ans_temp.PoI == TRUE)
    {
        printf("Check Successful \n");
    }
    
    EFp16_clear(&P);
    Fp16_clear(&a);
    Fp16_clear(&x);
    mpz_clear(t2);
    mpz_clear(p2);
    mpz_clear(p22);
    mpz_clear(p4);
    mpz_clear(p8);
    mpz_clear(tmp1);
    mpz_clear(tmp2);
    mpz_clear(t16);
}

void EFp16_frobenius_map(struct EFp16 *ANS,struct EFp16 *A){
    struct EFp16 tmp;
    EFp16_init(&tmp);
    
    Fp16_frobenius_map(&tmp.x,&A->x);
    Fp16_frobenius_map(&tmp.y,&A->y);
    
    EFp16_set(ANS,&tmp);
    
    EFp16_clear(&tmp);
}

void EFp16_random_set_G2(struct EFp16 *ANS){
    struct EFp16 P,P_frobenius,tmp_EFp16;
    EFp16_init(&P);
    EFp16_init(&P_frobenius);
    EFp16_init(&tmp_EFp16);
    
    EFp16_random_set(&P);
    
    EFp16_frobenius_map(&P_frobenius,&P);
    Fp16_neg(&tmp_EFp16.y,&P.y);
    Fp16_set(&tmp_EFp16.x,&P.x);
    
    EFp16_eca(&tmp_EFp16,&tmp_EFp16,&P_frobenius);
    //    EFp16_printf(&tmp_EFp16);
    EFp16_set(ANS,&tmp_EFp16);
    
    EFp16_clear(&P);
    EFp16_clear(&P_frobenius);
    EFp16_clear(&tmp_EFp16);
}


void EFp16_to_EFp4_map(struct EFp4 *ANS,struct EFp16 *A){
    Fp4_set_ui(&ANS->x,0);
    Fp4_set_ui(&ANS->y,0);
    Fp4_set(&ANS->x,&A->x.x0.x1);
    Fp4_set(&ANS->y,&A->y.x1.x1);
    ANS->PoI=A->PoI;
}

void EFp4_to_EFp16_map(struct EFp16 *ANS,struct EFp4 *A){
    Fp16_set_ui(&ANS->x,0);
    Fp16_set_ui(&ANS->y,0);
    Fp4_set(&ANS->x.x0.x1,&A->x);
    Fp4_set(&ANS->y.x1.x1,&A->y);
    ANS->PoI=A->PoI;
}
