//
//  bn_fp.c
//  bn_Single_File_Elips
//
//  Created by Y.Nanjo and M.A.A Khandaker on 1/25/18.
//  Copyright © 2018 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_bn_bls/bn_fp.h>

#define c1 2
#define C1_INV
void Fp_init(Fp *A){
    mpz_init(A->x0);
}

void Fp_clear(Fp *A){
    mpz_clear(A->x0);
}

void Fp_printf(Fp *A,char *str){
    gmp_printf("%s%Zd",str,A->x0);
}

void Fp_set(Fp *ANS,Fp *A){
    mpz_set(ANS->x0,A->x0);
}

void Fp_set_ui(Fp *ANS,unsigned long int UI){
    mpz_set_ui(ANS->x0,UI);
}

void Fp_set_mpz(Fp *ANS,mpz_t A){
    mpz_set(ANS->x0,A);
}

void Fp_set_neg(Fp *ANS,Fp *A){
    mpz_neg(ANS->x0,A->x0);
}

void Fp_set_random(Fp *ANS,gmp_randstate_t state){
    mpz_urandomm(ANS->x0,state,curve_parameters.prime);
}

void Fp_mul(Fp *ANS,Fp *A,Fp *B){
    mpz_mul(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,curve_parameters.prime);
}

void Fp_mul_ui(Fp *ANS,Fp *A,unsigned long int UI){
    mpz_mul_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,curve_parameters.prime);
}

void Fp_mul_mpz(Fp *ANS,Fp *A,mpz_t B){
    mpz_mul(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,curve_parameters.prime);
}

void Fp_mul_basis(Fp *ANS,Fp *A){
    mpz_sub(ANS->x0,curve_parameters.prime,A->x0);
}

void Fp_mul_basis_KSS16(Fp *ANS,Fp *A){
    mpz_mul_ui(ANS->x0,A->x0,c1);
    mpz_mod(ANS->x0,ANS->x0,curve_parameters.curve_a);
}

void Fp_add(Fp *ANS,Fp *A,Fp *B){
    mpz_add(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,curve_parameters.prime);
}

void Fp_add_ui(Fp *ANS,Fp *A,unsigned long int UI){
    mpz_add_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,curve_parameters.prime);
}

void Fp_add_mpz(Fp *ANS,Fp *A,mpz_t B){
    mpz_add(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,curve_parameters.prime);
}

void Fp_sub(Fp *ANS,Fp *A,Fp *B){
    mpz_sub(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,curve_parameters.prime);
}

void Fp_sub_ui(Fp *ANS,Fp *A,unsigned long int UI){
    mpz_sub_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,curve_parameters.prime);
}

void Fp_sub_mpz(Fp *ANS,Fp *A,mpz_t B){
    mpz_sub(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,curve_parameters.prime);
}

void Fp_inv(Fp *ANS,Fp *A){
    mpz_invert(ANS->x0,A->x0,curve_parameters.prime);
}
void Fp_div(Fp *ANS, Fp *A, Fp *B){
    mpz_invert(ANS->x0,B->x0,curve_parameters.prime);
    mpz_mul(ANS->x0,A->x0,ANS->x0);
    mpz_mod(ANS->x0,ANS->x0,curve_parameters.prime);
}

int  Fp_legendre(Fp *A){
    return mpz_legendre(A->x0,curve_parameters.prime);
}

int  Fp_isCNR(Fp *A){
    Fp tmp;
    Fp_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_sub_ui(exp,curve_parameters.prime,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp_pow(&tmp,A,exp);
    
    if(Fp_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp_clear(&tmp);
        return -1;
    }
}

void Fp_sqrt(Fp *ANS,Fp *A){
    Fp x,y,t,k,n,tmp;
    Fp_init(&x);
    Fp_init(&y);
    Fp_init(&t);
    Fp_init(&k);
    Fp_init(&n);
    Fp_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    Fp_set_random(&n,state);
    
    while(Fp_legendre(&n)!=-1){
        Fp_set_random(&n,state);
    }
    mpz_sub_ui(q,curve_parameters.prime,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp_pow(&x,A,exp);
    Fp_mul(&tmp,&x,&x);
    Fp_mul(&k,&tmp,A);
    Fp_mul(&x,&x,A);
    while(mpz_cmp_ui(k.x0,1)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp_pow(&tmp,&k,exp);
        while(mpz_cmp_ui(tmp.x0,1)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp_pow(&t,&y,result);
        Fp_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp_mul(&x,&x,&t);
        Fp_mul(&k,&k,&y);
    }
    Fp_set(ANS,&x);
    
    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
    Fp_clear(&x);
    Fp_clear(&y);
    Fp_clear(&t);
    Fp_clear(&k);
    Fp_clear(&n);
    Fp_clear(&tmp);
}

void Fp_pow(Fp *ANS,Fp *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp tmp;
    Fp_init(&tmp);
    
    Fp_set(&tmp,A);
    
    for(i=1; i<length; i++){
        Fp_mul(&tmp,&tmp,&tmp);
        if(binary[i]=='1'){
            Fp_mul(&tmp,A,&tmp);
        }
    }
    Fp_set(ANS,&tmp);
    
    Fp_clear(&tmp);
}

int  Fp_cmp(Fp *A,Fp *B){
    if(mpz_cmp(A->x0,B->x0)==0){
        return 0;
    }
    return 1;
}

int  Fp_cmp_ui(Fp *A,unsigned long int UI){
    if(mpz_cmp_ui(A->x0,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp_cmp_mpz(Fp *A,mpz_t B){
    if(mpz_cmp(A->x0,B)==0){
        return 0;
    }
    return 1;
}

int  Fp_cmp_zero(Fp *A){
    if(mpz_cmp_ui(A->x0,0)==0){
        return 0;
    }
    return 1;
}

int  Fp_cmp_one(Fp *A){
    if(mpz_cmp_ui(A->x0,1)==0){
        return 0;
    }
    return 1;
}

void Fp_neg(struct Fp *ANS,struct Fp *A){
    mpz_sub(ANS->x0,curve_parameters.prime,A->x0);
}


