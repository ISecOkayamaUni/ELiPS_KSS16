//
//  bn_fp12.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 1/31/18.

/*
 * ELiPS is an Efficient Library for Pairing-based Systems
 * Copyright (C) 2008-2018 ELiPS Authors. Please refer to the COPYRIGHT file
 * for contact information.
 *
 * This file is part of ELiPS Project. ELiPS is legal property of its
 * developers, whose names are not listed here. 
 *
 * ELiPS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * ELiPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ELiPS. If not, see <http://www.gnu.org/licenses/>.
 */

#include <ELiPS_bn_bls/bn_fp12.h>


void Fp12_init(Fp12 *A){
    Fp6_init(&A->x0);
    Fp6_init(&A->x1);
}

void Fp12_clear(Fp12 *A){
    Fp6_clear(&A->x0);
    Fp6_clear(&A->x1);
}

void Fp12_printf(Fp12 *A,char *str){
    gmp_printf("%s(",str);
    Fp6_printf(&A->x0,"");
    gmp_printf(",");
    Fp6_printf(&A->x1,"");
    gmp_printf(")");
}

void Fp12_set(Fp12 *ANS,Fp12 *A){
    Fp6_set(&ANS->x0,&A->x0);
    Fp6_set(&ANS->x1,&A->x1);
}

void Fp12_set_ui(Fp12 *ANS,unsigned long int UI){
    Fp6_set_ui(&ANS->x0,UI);
    Fp6_set_ui(&ANS->x1,UI);
}

void Fp12_set_mpz(Fp12 *ANS,mpz_t A){
    Fp6_set_mpz(&ANS->x0,A);
    Fp6_set_mpz(&ANS->x1,A);
}

void Fp12_set_neg(Fp12 *ANS,Fp12 *A){
    Fp6_set_neg(&ANS->x0,&A->x0);
    Fp6_set_neg(&ANS->x1,&A->x1);
}

void Fp12_set_random(Fp12 *ANS,gmp_randstate_t state){
    Fp6_set_random(&ANS->x0,state);
    Fp6_set_random(&ANS->x1,state);
}

void Fp12_mul(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6 tmp1,tmp2;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    
    //set
    Fp6_mul(&tmp2,&A->x1,&B->x1);//b*d
    Fp6_add(&tmp1,&A->x0,&A->x1);//a+b
    Fp6_add(&ANS->x1,&B->x0,&B->x1);//c+d
    Fp6_mul(&ANS->x1,&tmp1,&ANS->x1);//(a+b)(c+d)
    Fp6_mul(&tmp1,&A->x0,&B->x0);//a*c
    //x0
    Fp6_mul_basis(&ANS->x0,&tmp2);//b*d*v
    Fp6_add(&ANS->x0,&ANS->x0,&tmp1);//a*c+b*d*v
    //x1
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp1);
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp2);
    
    //clear
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
}

void Fp12_mul_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_mul_ui(&ANS->x0,&A->x0,UI);
    Fp6_mul_ui(&ANS->x1,&A->x1,UI);
}

void Fp12_mul_mpz(Fp12 *ANS,Fp12 *A,mpz_t B){
    Fp6_mul_mpz(&ANS->x0,&A->x0,B);
    Fp6_mul_mpz(&ANS->x1,&A->x1,B);
}

void Fp12_sqr(Fp12 *ANS,Fp12 *A){
    Fp6 tmp1,tmp2,tmp3;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    Fp6_init(&tmp3);
    
    Fp6_add(&tmp1,&A->x0,&A->x1);
    Fp6_mul_basis(&tmp2,&A->x1);
    Fp6_add(&tmp2,&tmp2,&A->x0);
    Fp6_mul(&tmp3,&A->x0,&A->x1);
    
    //x0
    Fp6_mul(&ANS->x0,&tmp1,&tmp2);
    Fp6_sub(&ANS->x0,&ANS->x0,&tmp3);
    Fp6_mul_basis(&tmp1,&tmp3);
    Fp6_sub(&ANS->x0,&ANS->x0,&tmp1);
    //x1
    Fp6_add(&ANS->x1,&tmp3,&tmp3);
    
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
    Fp6_clear(&tmp3);
}

void Fp12_add(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6_add(&ANS->x0,&A->x0,&B->x0);
    Fp6_add(&ANS->x1,&A->x1,&B->x1);
}

void Fp12_add_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_add_ui(&ANS->x0,&A->x0,UI);
    Fp6_add_ui(&ANS->x1,&A->x1,UI);
}

void Fp12_add_mpz(Fp12 *ANS,Fp12 *A,mpz_t B){
    Fp6_add_mpz(&ANS->x0,&ANS->x0,B);
    Fp6_add_mpz(&ANS->x1,&ANS->x1,B);
}

void Fp12_sub(Fp12 *ANS,Fp12 *A,Fp12 *B){
    Fp6_sub(&ANS->x0,&A->x0,&B->x0);
    Fp6_sub(&ANS->x1,&A->x1,&B->x1);
}

void Fp12_sub_ui(Fp12 *ANS,Fp12 *A,unsigned long int UI){
    Fp6_sub_ui(&ANS->x0,&ANS->x0,UI);
    Fp6_sub_ui(&ANS->x1,&ANS->x1,UI);
}

void Fp12_sub_mpz(Fp12 *ANS,Fp12 *A,mpz_t B){
    Fp6_sub_mpz(&ANS->x0,&ANS->x0,B);
    Fp6_sub_mpz(&ANS->x1,&ANS->x1,B);
}

void Fp12_inv(Fp12 *ANS,Fp12 *A){
    Fp12 frob,tmp;
    Fp12_init(&frob);
    Fp12_init(&tmp);
    
    Fp6_set(&frob.x0,&A->x0);
    Fp6_set_neg(&frob.x1,&A->x1);
    
    Fp12_mul(&tmp,A,&frob);
    Fp6_inv(&tmp.x0,&tmp.x0);
    Fp6_mul(&ANS->x0,&frob.x0,&tmp.x0);
    Fp6_mul(&ANS->x1,&frob.x1,&tmp.x0);
    
    Fp12_clear(&frob);
    Fp12_clear(&tmp);
}

int  Fp12_legendre(Fp12 *A){
    mpz_t exp;
    mpz_init(exp);
    Fp12 tmp;
    Fp12_init(&tmp);
    
    mpz_pow_ui(exp,curve_parameters.prime,12);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp12_pow(&tmp,A,exp);
    
    if(Fp12_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp12_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp12_clear(&tmp);
        return -1;
    }
}

void Fp12_sqrt(Fp12 *ANS,Fp12 *A){
    Fp12 x,y,t,k,n,tmp;
    Fp12_init(&x);
    Fp12_init(&y);
    Fp12_init(&t);
    Fp12_init(&k);
    Fp12_init(&n);
    Fp12_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    Fp12_set_random(&n,state);
    while(Fp12_legendre(&n)!=-1){
        Fp12_set_random(&n,state);
    }
    mpz_pow_ui(q,curve_parameters.prime,12);
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp12_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp12_pow(&x,A,exp);
    Fp12_mul(&tmp,&x,&x);
    Fp12_mul(&k,&tmp,A);
    Fp12_mul(&x,&x,A);
    while(Fp12_cmp_one(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp12_pow(&tmp,&k,exp);
        while(Fp12_cmp_one(&tmp)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp12_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp12_pow(&t,&y,result);
        Fp12_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp12_mul(&x,&x,&t);
        Fp12_mul(&k,&k,&y);
    }
    Fp12_set(ANS,&x);
    
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(exp);
    mpz_clear(result);
    Fp12_clear(&x);
    Fp12_clear(&y);
    Fp12_clear(&t);
    Fp12_clear(&k);
    Fp12_clear(&n);
    Fp12_clear(&tmp);
}

void Fp12_pow(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp12 tmp;
    Fp12_init(&tmp);
    Fp12_set(&tmp,A);
    
    for(i=1; i<length; i++){
        Fp12_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp12_mul(&tmp,A,&tmp);
        }
    }
    
    Fp12_set(ANS,&tmp);
    Fp12_clear(&tmp);
}

int  Fp12_cmp(Fp12 *A,Fp12 *B){
    if(Fp6_cmp(&A->x0,&B->x0)==0 && Fp6_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}

int  Fp12_cmp_ui(Fp12 *A,unsigned long int UI){
    if(Fp6_cmp_ui(&A->x0,UI)==0 && Fp6_cmp_ui(&A->x1,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp12_cmp_mpz(Fp12 *A,mpz_t B){
    if(Fp6_cmp_mpz(&A->x0,B)==0 && Fp6_cmp_mpz(&A->x1,B)==0){
        return 0;
    }
    return 1;
}

int  Fp12_cmp_zero(Fp12 *A){
    if(Fp6_cmp_zero(&A->x0)==0 && Fp6_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

int  Fp12_cmp_one(Fp12 *A){
    if(Fp6_cmp_one(&A->x0)==0 && Fp6_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

