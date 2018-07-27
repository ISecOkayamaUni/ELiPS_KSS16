//
//  fp8.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/26/18.

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

#include <ELiPS_bn_bls/fp8.h>

void Fp8_init(struct Fp8 *A){
    Fp4_init(&A->x0);
    Fp4_init(&A->x1);
}
void Fp8_set(struct Fp8 *ANS,struct Fp8 *A){
    Fp4_set(&ANS->x0,&A->x0);
    Fp4_set(&ANS->x1,&A->x1);
}
void Fp8_set_ui(struct Fp8 *A,signed long int B){
    Fp4_set_ui(&A->x0,B);
    Fp4_set_ui(&A->x1,B);
}
//void Fp8_random(struct Fp8 *A){
//    Fp4_random(&A->x0);
//    Fp4_random(&A->x1);
//}
void Fp8_set_random(Fp8 *A, gmp_randstate_t state){
    Fp4_set_random(&A->x0,state);
    Fp4_set_random(&A->x1,state);
}
void Fp8_clear(struct Fp8 *A){
    Fp4_clear(&A->x0);
    Fp4_clear(&A->x1);
}
void Fp8_printf(struct Fp8 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd\n",A->x0.x0.x0.x0,A->x0.x0.x1.x0,A->x0.x1.x0.x0,A->x0.x1.x1.x0);
    gmp_printf("(%Zd,%Zd,%Zd,%Zd\n",A->x1.x0.x0.x0,A->x1.x0.x1.x0,A->x1.x1.x0.x0,A->x1.x1.x1.x0);
}
void Fp8_add(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B){
    Fp4_add(&ANS->x0,&A->x0,&B->x0);
    Fp4_add(&ANS->x1,&A->x1,&B->x1);
}
void Fp8_add_ui(struct Fp8 *ANS,struct Fp8 *A,unsigned long int B){
    Fp4_add_ui(&ANS->x0,&A->x0,B);
    Fp4_add_ui(&ANS->x1,&A->x1,B);
}
void Fp8_sub(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B){
    Fp4_sub(&ANS->x0,&A->x0,&B->x0);
    Fp4_sub(&ANS->x1,&A->x1,&B->x1);
}

void Fp8_sqr(struct Fp8 *ANS,struct Fp8 *A){
    
    struct Fp4 tmp1,tmp2,tmp3;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    
    Fp4_add(&tmp1,&A->x0,&A->x1);
    Fp4_mul_basis(&tmp2,&A->x1);
    Fp4_add(&tmp2,&tmp2,&A->x0);
    Fp4_mul(&tmp3,&A->x0,&A->x1);
    
    //x0
    Fp4_mul(&ANS->x0,&tmp1,&tmp2);
    Fp4_sub(&ANS->x0,&ANS->x0,&tmp3);
    Fp4_mul_basis(&tmp1,&tmp3);
    Fp4_sub(&ANS->x0,&ANS->x0,&tmp1);
    //x1
    Fp4_add(&ANS->x1,&tmp3,&tmp3);
    
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    
}

void Fp8_mul(struct Fp8 *ANS,struct Fp8 *A, struct Fp8 *B){
    struct Fp4 tmp1,tmp2;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    
    //set
    Fp4_mul(&tmp2,&A->x1,&B->x1);//b*d
    Fp4_add(&tmp1,&A->x0,&A->x1);//a+b
    Fp4_add(&ANS->x1,&B->x0,&B->x1);//c+d
    Fp4_mul(&ANS->x1,&tmp1,&ANS->x1);//(a+b)(c+d)
    Fp4_mul(&tmp1,&A->x0,&B->x0);//a*c
    //x0
    Fp4_mul_basis(&ANS->x0,&tmp2);//b*d*v
    Fp4_add(&ANS->x0,&ANS->x0,&tmp1);//a*c+b*d*v
    //x1
    Fp4_sub(&ANS->x1,&ANS->x1,&tmp1);
    Fp4_sub(&ANS->x1,&ANS->x1,&tmp2);
    
    //clear
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
}
void Fp8_mul_v(struct Fp8 *ANS,struct Fp8 *A){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    Fp8_set(&tmp,A);
    
    Fp_mul_basis(&ANS->x0.x0.x0,&tmp.x1.x1.x1);
    Fp_set(&ANS->x0.x0.x1,&tmp.x1.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&tmp.x1.x0.x0);
    Fp_set(&ANS->x0.x1.x1,&tmp.x1.x0.x1);
    Fp_set(&ANS->x1.x0.x0,&tmp.x0.x0.x0);
    Fp_set(&ANS->x1.x0.x1,&tmp.x0.x0.x1);
    Fp_set(&ANS->x1.x1.x0,&tmp.x0.x1.x0);
    Fp_set(&ANS->x1.x1.x1,&tmp.x0.x1.x1);
    
    Fp8_clear(&tmp);
}
void Fp8_mul_ui(struct Fp8 *ANS,struct Fp8 *A,unsigned long int B){
    Fp4_mul_ui(&ANS->x0,&A->x0,B);
    Fp4_mul_ui(&ANS->x1,&A->x1,B);
}
void Fp8_mul_Fp(struct Fp8 *ANS,struct Fp8 *A,struct Fp *B){
    Fp4_mul_Fp(&ANS->x0,&A->x0,B);
    Fp4_mul_Fp(&ANS->x1,&A->x1,B);
}
void Fp8_mul_mpz(struct Fp8 *ANS,struct Fp8 *A,mpz_t B){
    Fp4_mul_mpz(&ANS->x0,&A->x0,B);
    Fp4_mul_mpz(&ANS->x1,&A->x1,B);
}
void Fp8_invert(struct Fp8 *ANS,struct Fp8 *A){
    struct Fp8 frob,buf;
    Fp8_init(&frob);
    Fp8_init(&buf);
    
    Fp4_set(&frob.x0,&A->x0);
    Fp4_neg(&frob.x1,&A->x1);
    Fp8_mul(&buf,A,&frob);
    Fp4_invert(&buf.x0,&buf.x0);
    Fp4_mul(&ANS->x0,&frob.x0,&buf.x0);
    Fp4_mul(&ANS->x1,&frob.x1,&buf.x0);
    
    Fp8_clear(&frob);
    Fp8_clear(&buf);
}
void Fp8_div(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    Fp8_invert(&tmp,B);
    Fp8_mul(ANS,A,&tmp);
    
    Fp8_clear(&tmp);
}
void Fp8_pow(struct Fp8 *ANS,struct Fp8 *A,mpz_t B){
    int i,length;
    length=(int)mpz_sizeinbase(B,2);
    char binary[length];
    mpz_get_str(binary,2,B);
    struct Fp8 buf;
    Fp8_init(&buf);
    Fp8_set(&buf,A);
    
    for(i=1; i<length; i++){
        //Fp2_mul(&buf,&buf,&buf);
        Fp8_sqr(&buf,&buf);
        if(binary[i]=='1'){
            Fp8_mul(&buf,A,&buf);
        }
    }
    
    Fp8_set(ANS,&buf);
    Fp8_clear(&buf);
}
void Fp8_sqrt(struct Fp8 *ANS,struct Fp8 *A){
    struct Fp8 n,y,x,b,t,tmp_Fp4;
    Fp8_init(&n);
    Fp8_init(&y);
    Fp8_init(&x);
    Fp8_init(&b);
    Fp8_init(&t);
    Fp8_init(&tmp_Fp4);
    Fp8_set(&n,A);
    
    mpz_t tmp_mpz,q,e,r,set_1,set_2;
    mpz_init(tmp_mpz);
    mpz_init(q);
    mpz_init(e);
    mpz_init(r);
    mpz_init(set_1);
    mpz_init(set_2);
    mpz_set_ui(set_1,1);
    mpz_set_ui(set_2,2);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    Fp8_set_random(&n,state);
    while(Fp8_legendre(&n)!=-1){
        Fp8_set_random(&n,state);
    }
    
    mpz_pow_ui(q,curve_parameters.prime,12);
    mpz_sub_ui(q,q,1);
    mpz_set_ui(e,0);
    while(mpz_odd_p(q)==0){
        mpz_add_ui(e,e,1);
        mpz_div_ui(q,q,2);
    }
    Fp8_pow(&y,&n,q);
    
    mpz_set(r,e);
    
    mpz_sub_ui(tmp_mpz,q,1);
    mpz_div_ui(tmp_mpz,tmp_mpz,2);
    
    Fp8_pow(&x,A,tmp_mpz);
    Fp8_pow(&tmp_Fp4,&x,set_2);
    Fp8_mul(&b,&tmp_Fp4,A);
    Fp8_mul(&x,&x,A);
    
    int m;
    
    while(Fp8_cmp_mpz(&b,set_1)==1){
        m=-1;
        Fp8_set(&tmp_Fp4,&b);
        while(Fp8_cmp_mpz(&tmp_Fp4,set_1)==1){
            m++;
            mpz_pow_ui(tmp_mpz,set_2,m);
            Fp8_pow(&tmp_Fp4,&b,tmp_mpz);
        }
        mpz_sub_ui(tmp_mpz,r,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,curve_parameters.prime);
        // gmp_printf("%Zd,%Zd,%d\n",tmp_mpz,r,m);
        Fp8_pow(&t,&y,tmp_mpz);
        Fp8_pow(&y,&t,set_2);
        // gmp_printf("%Zd,%Zd,\n",y.x0.x0.x0,y.x0.x1.x0);
        mpz_set_ui(r,m);
        Fp8_mul(&x,&x,&t);
        Fp8_mul(&b,&b,&y);
    }
    
    Fp8_set(ANS,&x);
    
    Fp8_clear(&n);
    Fp8_clear(&y);
    Fp8_clear(&x);
    Fp8_clear(&b);
    Fp8_clear(&t);
    Fp8_clear(&tmp_Fp4);
    mpz_clear(tmp_mpz);
    mpz_clear(q);
    mpz_clear(e);
    mpz_clear(r);
    mpz_clear(set_1);
}
int Fp8_legendre(struct Fp8 *a){
    mpz_t i,cmp;
    struct Fp8 tmp;
    Fp8_init(&tmp);
    mpz_init(i);
    mpz_init(cmp);
    mpz_set_ui(cmp,1);
    mpz_pow_ui(i,curve_parameters.prime,8);
    mpz_sub_ui(i,i,1);
    mpz_tdiv_q_ui(i,i,2);
    Fp8_pow(&tmp,a,i);
    
    if((Fp8_cmp_mpz(&tmp,cmp))==0){
        Fp8_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return 1;
    }else{
        Fp8_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return -1;
    }
}
int Fp8_cmp(struct Fp8 *A,struct Fp8 *B){
    if(Fp4_cmp(&A->x0,&B->x0)==0 && Fp4_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}
int Fp8_cmp_mpz(struct Fp8 *A,mpz_t B){
    if(Fp4_cmp_mpz(&A->x0,B)==0 && Fp4_cmp_mpz(&A->x1,B)==0){
        return 0;
    }
    return 1;
}
void Fp8_neg(struct Fp8 *ans,struct Fp8 *a){
    Fp4_neg(&ans->x0,&a->x0);
    Fp4_neg(&ans->x1,&a->x1);
}
