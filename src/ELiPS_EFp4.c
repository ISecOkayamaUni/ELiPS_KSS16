//
//  EFp4.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_EFp4.h>

#pragma mark EFp4 method implementations

void EFp4_init(struct EFp4 *A){
    Fp4_init(&A->x);
    Fp4_init(&A->y);
    A->PoI=FALSE;
}
void EFp4_set(struct EFp4 *A,struct EFp4 *B){
    Fp4_set(&A->x,&B->x);
    Fp4_set(&A->y,&B->y);
    A->PoI=B->PoI;
}
void EFp4_set_PoI(struct EFp4 *A){
    Fp4_set_ui(&A->x,0);
    Fp4_set_ui(&A->y,0);
    A->PoI=TRUE;
}
void EFp4_set_EFp(struct EFp4 *ANS,struct EFp *A){
    Fp4_set_ui(&ANS->x,0);
    Fp4_set_ui(&ANS->y,0);
    
    Fp_set(&ANS->x.x0.x0,&A->x);
    Fp_set(&ANS->y.x0.x0,&A->y);
    ANS->PoI=A->PoI;
}
void EFp4_clear(struct EFp4 *A){
    Fp4_clear(&A->x);
    Fp4_clear(&A->y);
}
void EFp4_printf(struct EFp4 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd)\n",A->x.x0.x0.x0,A->x.x0.x1.x0,A->x.x1.x0.x0, A->x.x1.x1.x0);
    gmp_printf("(%Zd,%Zd,%Zd,%Zd)\n",A->y.x0.x0.x0,A->y.x0.x1.x0,A->y.x1.x0.x0, A->y.x1.x1.x0);
}
void EFp4_eca(struct EFp4 *ANS, struct EFp4 *P1, struct EFp4 *P2){
    if(P2->PoI==TRUE){//if P2==inf
        EFp4_set(ANS,P1);
        return;
    }
    else if(P1->PoI==TRUE){//if P1==inf
        EFp4_set(ANS,P2);
        return;
    }
    else if(Fp4_cmp(&P1->x,&P2->x)==0&&Fp4_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp4_set_PoI(ANS);
        return;
    }
    else if(EFp4_cmp(P1,P2)==0){ // P=Q
        EFp4_ecd(ANS,P1);
        return;
    }
    
    struct Fp4 x,y,lambda,tmp;
    struct EFp4 t_ans;
    
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&lambda);
    Fp4_init(&tmp);
    EFp4_init(&t_ans);
    
    Fp4_sub(&x,&P2->x,&P1->x);
    Fp4_sub(&y,&P2->y,&P1->y);
    Fp4_div(&lambda,&y,&x);
    Fp4_mul(&tmp,&lambda,&lambda);
    Fp4_add(&x,&P1->x,&P2->x);
    Fp4_sub(&x,&tmp,&x);
    Fp4_sub(&tmp,&P1->x,&x);
    Fp4_set(&t_ans.x,&x);
    Fp4_mul(&tmp,&tmp,&lambda);
    Fp4_sub(&t_ans.y,&tmp,&P1->y);
    
    EFp4_set(ANS,&t_ans);
    
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&lambda);
    Fp4_clear(&tmp);
    EFp4_clear(&t_ans);
}
void EFp4_ecd(struct EFp4 *ANS, struct EFp4 *P){
    if(P->PoI==TRUE){
        EFp4_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp4_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp4_set_PoI(ANS);
        return;
    }
    
    struct Fp4 x,y,lambda,tmp,a_betainv;
    struct EFp4 t_ans;
    Fp4_init(&x);
    Fp4_init(&lambda);
    Fp4_init(&tmp);
    Fp4_init(&y);
    Fp4_init(&a_betainv);
    
    EFp4_init(&t_ans);
    
    
    Fp4_mul(&x,&P->x,&P->x);
    Fp4_add(&tmp,&x,&x);
    Fp4_add(&x,&tmp,&x);
    
    Fp_add_mpz(&x.x0.x0,&x.x0.x0,kss_curve_const.a);
    Fp4_add(&y,&P->y,&P->y);
    Fp4_div(&lambda,&x,&y);
    Fp4_mul(&tmp,&lambda,&lambda);
    Fp4_add(&x,&P->x,&P->x);
    Fp4_sub(&x,&tmp,&x);
    Fp4_sub(&tmp,&P->x,&x);
    Fp4_set(&t_ans.x,&x);
    Fp4_mul(&tmp,&tmp,&lambda);
    Fp4_sub(&t_ans.y,&tmp,&P->y);
    
    EFp4_set(ANS,&t_ans);
    
    Fp4_clear(&x);
    Fp4_clear(&lambda);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&t_ans);
}
void EFp4_ecd_sparse(struct EFp4 *ANS, struct EFp4 *P){
    if(P->PoI==TRUE){
        EFp4_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp4_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp4_set_PoI(ANS);
        return;
    }
    
    struct Fp4 x,y,lambda,tmp,a_betainv;
    struct EFp4 t_ans;
    Fp4_init(&x);
    Fp4_init(&lambda);
    Fp4_init(&tmp);
    Fp4_init(&y);
    Fp4_init(&a_betainv);
    
    EFp4_init(&t_ans);
    
    
    Fp4_mul(&x,&P->x,&P->x);
    Fp4_add(&tmp,&x,&x);
    Fp4_add(&x,&tmp,&x);
    
    Fp4_mul_beta_inv(&a_betainv);
    Fp4_add(&x, &x, &a_betainv);
    
    // Fp_add_mpz(&x.x0.x0,&x.x0.x0,a_x);
    Fp4_add(&y,&P->y,&P->y);
    Fp4_div(&lambda,&x,&y);
    Fp4_mul(&tmp,&lambda,&lambda);
    Fp4_add(&x,&P->x,&P->x);
    Fp4_sub(&x,&tmp,&x);
    Fp4_sub(&tmp,&P->x,&x);
    Fp4_set(&t_ans.x,&x);
    Fp4_mul(&tmp,&tmp,&lambda);
    Fp4_sub(&t_ans.y,&tmp,&P->y);
    
    EFp4_set(ANS,&t_ans);
    
    Fp4_clear(&x);
    Fp4_clear(&lambda);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&t_ans);
}

int EFp4_cmp(struct EFp4 *A,struct EFp4 *B){
    if(Fp4_cmp(&A->x,&B->x)==0 && Fp4_cmp(&A->y,&B->y)==0){
        return 0;
    }
    return 1;
}

void EFp4_scm_bin(struct EFp4 *ANS,struct EFp4 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp4 Q,R;
    EFp4_init(&Q);
    EFp4_set(&Q,P);
    EFp4_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp4_ecd(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp4_eca(&Q,&Q,P);
        }
    }
    EFp4_set(ANS,&Q);
    
    EFp4_clear(&Q);
    EFp4_clear(&R);
    return;
}
void EFp4_scm_win(struct EFp4 *ANS, struct EFp4 *P, mpz_t scalar){
    int i,length;
    length= (int)mpz_sizeinbase(scalar,2);
    char r_binary[length];
    mpz_get_str(r_binary,2,scalar);
    
    int win_size = 4;
    int rem = length % win_size;
    int mul3 = length - rem;
    
    struct EFp4 R,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15,ANS_temp;
    EFp4_init(&R);
    EFp4_init(&R2);
    EFp4_init(&R3);
    EFp4_init(&R4);
    EFp4_init(&R5);
    EFp4_init(&R6);
    EFp4_init(&R7);
    EFp4_init(&R8);
    EFp4_init(&R9);
    EFp4_init(&R10);
    EFp4_init(&R11);
    EFp4_init(&R12);
    EFp4_init(&R13);
    EFp4_init(&R14);
    EFp4_init(&R15);
    EFp4_init(&ANS_temp);
    
    EFp4_set_PoI(&ANS_temp);
    EFp4_set(&R, P);
    EFp4_eca(&R2, &R, &R);
    EFp4_eca(&R3, &R2, &R);
    EFp4_eca(&R4, &R3, &R);
    EFp4_eca(&R5, &R4, &R);
    EFp4_eca(&R6, &R5, &R);
    EFp4_eca(&R7, &R6, &R);
    EFp4_eca(&R8, &R7, &R);
    EFp4_eca(&R9, &R8, &R);
    EFp4_eca(&R10, &R9, &R);
    EFp4_eca(&R11, &R10, &R);
    EFp4_eca(&R12, &R11, &R);
    EFp4_eca(&R13, &R12, &R);
    EFp4_eca(&R14, &R13, &R);
    EFp4_eca(&R15, &R14, &R);
    
    for(i=0; i< mul3;i=i+4){
        
        EFp4_ecd(&ANS_temp,&ANS_temp);
        EFp4_ecd(&ANS_temp,&ANS_temp);
        EFp4_ecd(&ANS_temp,&ANS_temp);
        EFp4_ecd(&ANS_temp,&ANS_temp);
        
        if (i+2 < length)
        {
            if(r_binary[i] =='1' && r_binary[i+1] =='1' && r_binary[i+2] =='1' && r_binary[i+3] =='1'){
                EFp4_eca(&ANS_temp, &ANS_temp, &R15);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='1' && r_binary[i+2] =='1' && r_binary[i+3] =='0'){
                EFp4_eca(&ANS_temp, &ANS_temp, &R14);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='1' && r_binary[i+2] =='0' && r_binary[i+3] =='1'){
                EFp4_eca(&ANS_temp, &ANS_temp, &R13);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='1' && r_binary[i+2] =='0' && r_binary[i+3] =='0'){
                EFp4_eca(&ANS_temp, &ANS_temp, &R12);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='0' && r_binary[i+2] =='1' && r_binary[i+3] =='1'){
                EFp4_eca(&ANS_temp, &ANS_temp, &R11);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='0' && r_binary[i+2] =='1' && r_binary[i+3] =='0'){
                EFp4_eca(&ANS_temp, &ANS_temp, &R10);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='0' && r_binary[i+2] =='0' && r_binary[i+3] =='1'){
                EFp4_eca(&ANS_temp, &ANS_temp, &R9);
            }
            else if(r_binary[i] =='1' && r_binary[i+1] =='0' && r_binary[i+2] =='0' && r_binary[i+3] =='0'){
                EFp4_eca(&ANS_temp, &ANS_temp, &R8);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='1' && r_binary[i+2] =='1' && r_binary[i+3] =='1'){
                EFp4_eca(&ANS_temp, &ANS_temp, &R7);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='1' && r_binary[i+2] =='1' && r_binary[i+3] =='0'){
                EFp4_eca(&ANS_temp, &ANS_temp, &R6);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='1' && r_binary[i+2] =='0' && r_binary[i+3] =='1'){
                EFp4_eca(&ANS_temp, &ANS_temp, &R5);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='1' && r_binary[i+2] =='0' && r_binary[i+3] =='0'){
                EFp4_eca(&ANS_temp, &ANS_temp, &R4);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='0' && r_binary[i+2] =='1' && r_binary[i+3] =='1'){
                EFp4_eca(&ANS_temp, &ANS_temp, &R3);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='0' && r_binary[i+2] =='1' && r_binary[i+3] =='0'){
                EFp4_eca(&ANS_temp, &ANS_temp, &R2);
            }
            else if(r_binary[i] =='0' && r_binary[i+1] =='0' && r_binary[i+2] =='0' && r_binary[i+3] =='1'){
                EFp4_eca(&ANS_temp, &ANS_temp, &R);
            }
            else{
            }
        }
    }
    
    for(i = mul3;r_binary[i]!='\0';i++){
        EFp4_ecd(&ANS_temp,&ANS_temp);
        if(r_binary[i]=='1'){
            EFp4_eca(&ANS_temp,&ANS_temp,P);
        }
    }
    EFp4_set(ANS,&ANS_temp);
    EFp4_clear(&ANS_temp);
    EFp4_clear(&R);
    EFp4_clear(&R2);
    EFp4_clear(&R3);
    EFp4_clear(&R4);
    EFp4_clear(&R5);
    EFp4_clear(&R6);
    EFp4_clear(&R7);
    return;
}
void EFp4_scm_ml(struct EFp4 *RES, struct EFp4 *P,mpz_t scalar){
    int i,length;
    length= (int)mpz_sizeinbase(scalar,2);
    char r_binary[length];
    mpz_get_str(r_binary,2,scalar);
    struct EFp4 T0,T1;
    EFp4_init(&T0);
    
    EFp4_set_PoI(&T0);
    EFp4_init(&T1);
    EFp4_set(&T1, P);
    //    EFp_ecd(&T1,&T1);
    for(i=0;r_binary[i]!='\0';i++){
        if(r_binary[i]=='1'){
            EFp4_eca(&T0,&T0,&T1);
            EFp4_ecd(&T1,&T1);
        }
        else if(r_binary[i]=='0'){
            EFp4_eca(&T1,&T0,&T1);
            EFp4_ecd(&T0,&T0);
        }
    }
    EFp4_set(RES,&T0);
    EFp4_clear(&T0);
    EFp4_clear(&T1);
    return;
}

void EFp4_scm_bin_sparse(struct EFp4 *ANS,struct EFp4 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp4 Q,R;
    EFp4_init(&Q);
    EFp4_set(&Q,P);
    EFp4_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp4_ecd_sparse(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp4_eca(&Q,&Q,P);
        }
    }
    EFp4_set(ANS,&Q);
    
    EFp4_clear(&Q);
    EFp4_clear(&R);
    return;
}

void kss16_skew_frobenius_map(struct EFp4 *ANS, struct EFp4 *Qt)
{
    struct EFp4 tmp_ans;
    EFp4_init(&tmp_ans);
    
    struct Fp4 Qt_x, Qt_y;
    Fp4_init(&Qt_x);
    Fp4_init(&Qt_y);
    Fp4_set(&tmp_ans.x, &Qt->x);
    Fp4_set(&tmp_ans.y, &Qt->y);
    
    
    Fp_mul(&Qt_x.x0.x0, &tmp_ans.x.x0.x1,&m_cpm5d8);
    Fp_mul(&Qt_x.x0.x1,&tmp_ans.x.x0.x0,&pm5d8);
    Fp_mul(&Qt_x.x1.x0, &tmp_ans.x.x1.x1, &m_cpm1d4pm5d8);
    Fp_mul(&Qt_x.x1.x1, &tmp_ans.x.x1.x0, &pm1d4pm5d8);
    
    Fp_mul(&Qt_y.x0.x0, &tmp_ans.y.x1.x1,&m_ccpm1d4pm5d8p13d16);
    Fp_mul(&Qt_y.x0.x1,&tmp_ans.y.x1.x0,&cpm1d4pm5d6pm13d16);
    Fp_mul(&Qt_y.x1.x0, &tmp_ans.y.x0.x0, &cpm5d8pm13d16);
    Fp_mul(&Qt_y.x1.x1, &tmp_ans.y.x0.x1, &m_cpm5d8pm13d16);
    
    Fp4_set(&tmp_ans.x, &Qt_x);
    Fp4_set(&tmp_ans.y, &Qt_y);
    
    EFp4_set(ANS,&tmp_ans);
    EFp4_clear(&tmp_ans);
    Fp4_clear(&Qt_x);
    Fp4_clear(&Qt_y);
}


//void EFp4_ECD_Pseudo_Sparse(struct EFp4 *ANS, struct EFp4 *P){
//    if(P->PoI==TRUE){
//        EFp4_set(ANS,P);
//        return;
//    }
//    mpz_t cmp;
//    mpz_init(cmp);
//    mpz_set_ui(cmp,0);
//    if(Fp4_cmp_mpz(&P->y,cmp)==0){//P.y==0
//        EFp4_set_PoI(ANS);
//        return;
//    }
//    
//    struct Fp4 x,y,lambda,tmp,a_betainv;
//    struct EFp4 t_ans;
//    Fp4_init(&x);
//    Fp4_init(&lambda);
//    Fp4_init(&tmp);
//    Fp4_init(&y);
//    Fp4_init(&a_betainv);
//    
//    EFp4_init(&t_ans);
//    
//    
//    Fp4_mul(&x,&P->x,&P->x);
//    Fp4_add(&tmp,&x,&x);
//    Fp4_add(&x,&tmp,&x);
//    
//    Fp4_mul_beta_inv(&a_betainv);
//    Fp4_mul(&a_betainv, &a_betainv, &z_inv2);
//    Fp4_add(&x, &x, &a_betainv);
//    
//    Fp4_add(&y,&P->y,&P->y);
//    Fp4_div(&lambda,&x,&y);
//    Fp4_mul(&tmp,&lambda,&lambda);
//    Fp4_add(&x,&P->x,&P->x);
//    Fp4_sub(&x,&tmp,&x);
//    Fp4_sub(&tmp,&P->x,&x);
//    Fp4_set(&t_ans.x,&x);
//    Fp4_mul(&tmp,&tmp,&lambda);
//    Fp4_sub(&t_ans.y,&tmp,&P->y);
//    
//    EFp4_set(ANS,&t_ans);
//    
//    Fp4_clear(&x);
//    Fp4_clear(&lambda);
//    Fp4_clear(&y);
//    Fp4_clear(&tmp);
//    EFp4_clear(&t_ans);
//}

//void EFp4_random_set(struct EFp4 *ANS){
//    struct EFp4 P;
//    EFp4_init(&P);
//    
//    struct Fp4 x,a;
//    Fp4_init(&a);
//    Fp4_init(&x);
//    
//    mpz_t t2,p2,p22,tmp,r2;
//    mpz_t set_3;
//    mpz_init(set_3);
//    mpz_set_ui(set_3,3);
//    
//    mpz_init(t2);
//    mpz_init(p2);
//    mpz_init(p22);
//    mpz_init(tmp);
//    mpz_init(r2);
//    
//    mpz_add(p22,params.prime,params.prime);
//    mpz_mul(p2,params.prime,params.prime);
//    mpz_mul(t2,params.trace_t,params.trace_t);
//    mpz_sub(tmp,t2,p22);
//    mpz_add_ui(r2,p2,1);
//    mpz_sub(r2,r2,tmp);
//    
//    do{
//        Fp4_random(&x);
//        Fp4_pow(&a,&x,set_3);
//        mpz_add(a.x0.x0.x0,a.x0.x0.x0,b);
//    }while(Fp4_legendre(&a)!=1);
//    Fp4_sqrt(&P.y,&a);
//    Fp4_set(&P.x,&x);
//    
//    mpz_t r12_div_r2;
//    mpz_init(r12_div_r2);
//    mpz_div(r12_div_r2,r2,order_r);
//    mpz_div(r12_div_r2,r12_div_r2,order_r);
//    
//    EFp4_scm_bin(ANS,&P,r12_div_r2);
//    
//    EFp4_clear(&P);
//    Fp4_clear(&a);
//    Fp4_clear(&x);
//}
