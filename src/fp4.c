//
//  fp4.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/22/18.

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

#include <ELiPS_bn_bls/fp4.h>

void Fp4_init(Fp4 *A){
    Fp2_init(&A->x0);
    Fp2_init(&A->x1);
}
void Fp4_set( Fp4 *ANS, Fp4 *A){
    Fp2_set(&ANS->x0,&A->x0);
    Fp2_set(&ANS->x1,&A->x1);
}
void Fp4_set_ui( Fp4 *A,signed long int B){
    Fp2_set_ui(&A->x0,B);
    Fp2_set_ui(&A->x1,B);
}
void Fp4_set_random(Fp4 *A, gmp_randstate_t state){
    Fp2_set_random(&A->x0, state);
    Fp2_set_random(&A->x1, state);
}
void Fp4_clear( Fp4 *A){
    Fp2_clear(&A->x0);
    Fp2_clear(&A->x1);
}
void Fp4_printf( Fp4 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd\n",A->x0.x0.x0,A->x0.x1.x0,A->x1.x0.x0,A->x1.x1.x0);
}
void Fp4_add( Fp4 *ANS, Fp4 *A, Fp4 *B){
    Fp2_add(&ANS->x0,&A->x0,&B->x0);
    Fp2_add(&ANS->x1,&A->x1,&B->x1);
}
void Fp4_add_ui( Fp4 *ANS, Fp4 *A,unsigned long int B){
    Fp2_add_ui(&ANS->x0,&A->x0,B);
    Fp2_add_ui(&ANS->x1,&A->x1,B);
}
void Fp4_sub( Fp4 *ANS, Fp4 *A, Fp4 *B){
    Fp2_sub(&ANS->x0,&A->x0,&B->x0);
    Fp2_sub(&ANS->x1,&A->x1,&B->x1);
}

void Fp4_sqr( Fp4 *ANS, Fp4 *A){
     Fp2 tmp1,tmp2,tmp3;
    Fp2_init(&tmp1);
    Fp2_init(&tmp2);
    Fp2_init(&tmp3);
    
    Fp2_add(&tmp1,&A->x0,&A->x1);
    Fp2_mul_basis_KSS16(&tmp2,&A->x1);
    Fp2_add(&tmp2,&tmp2,&A->x0);
    Fp2_mul(&tmp3,&A->x0,&A->x1);
    
    //x0
    Fp2_mul(&ANS->x0,&tmp1,&tmp2);
    Fp2_sub(&ANS->x0,&ANS->x0,&tmp3);
    Fp2_mul_basis_KSS16(&tmp1,&tmp3);
    Fp2_sub(&ANS->x0,&ANS->x0,&tmp1);
    //x1
    Fp2_add(&ANS->x1,&tmp3,&tmp3);
    
    Fp2_clear(&tmp1);
    Fp2_clear(&tmp2);
    Fp2_clear(&tmp3);
}

void Fp4_mul( Fp4 *ANS, Fp4 *A, Fp4 *B){
    struct Fp2 tmp1,tmp2;
    Fp2_init(&tmp1);
    Fp2_init(&tmp2);
    
    //set
    Fp2_mul(&tmp2,&A->x1,&B->x1);//b*d
    Fp2_add(&tmp1,&A->x0,&A->x1);//a+b
    Fp2_add(&ANS->x1,&B->x0,&B->x1);//c+d
    Fp2_mul(&ANS->x1,&tmp1,&ANS->x1);//(a+b)(c+d)
    Fp2_mul(&tmp1,&A->x0,&B->x0);//a*c
    //x0
    Fp2_mul_basis_KSS16(&ANS->x0,&tmp2);//b*d*v
    Fp2_add(&ANS->x0,&ANS->x0,&tmp1);//a*c+b*d*v
    //x1
    Fp2_sub(&ANS->x1,&ANS->x1,&tmp1);
    Fp2_sub(&ANS->x1,&ANS->x1,&tmp2);
    
    //clear
    Fp2_clear(&tmp1);
    Fp2_clear(&tmp2);
}
void Fp4_mul_basis( Fp4 *ANS, Fp4 *A){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    Fp4_set(&tmp,A);
    
    Fp_mul_basis_KSS16(&ANS->x0.x0,&tmp.x1.x1);
    Fp_set(&ANS->x0.x1,&tmp.x1.x0);
    Fp_set(&ANS->x1.x0,&tmp.x0.x0);
    Fp_set(&ANS->x1.x1,&tmp.x0.x1);
    
    Fp4_clear(&tmp);
}
void Fp4_mul_ui( Fp4 *ANS, Fp4 *A,unsigned long int B){
    Fp2_mul_ui(&ANS->x0,&A->x0,B);
    Fp2_mul_ui(&ANS->x1,&A->x1,B);
}
void Fp4_mul_mpz( Fp4 *ANS, Fp4 *A,mpz_t B){
    Fp2_mul_mpz(&ANS->x0,&A->x0,B);
    Fp2_mul_mpz(&ANS->x1,&A->x1,B);
}
void Fp4_mul_Fp( Fp4 *ANS, Fp4 *A, Fp *B){
    Fp2_mul_Fp(&ANS->x0,&A->x0,B);
    Fp2_mul_Fp(&ANS->x1,&A->x1,B);
}
void Fp4_invert( Fp4 *ANS, Fp4 *A){
    struct Fp4 frob,buf;
    Fp4_init(&frob);
    Fp4_init(&buf);
    
    Fp2_set(&frob.x0,&A->x0);
    Fp2_neg(&frob.x1,&A->x1);
    Fp4_mul(&buf,A,&frob);
    Fp2_invert_kss16(&buf.x0,&buf.x0);
    Fp2_mul(&ANS->x0,&frob.x0,&buf.x0);
    Fp2_mul(&ANS->x1,&frob.x1,&buf.x0);
    
    Fp4_clear(&frob);
    Fp4_clear(&buf);
}
void Fp4_div( Fp4 *ANS, Fp4 *A, Fp4 *B){
    struct Fp4 tmp,t_ans;
    Fp4_init(&tmp);
    Fp4_init(&t_ans);
    
    Fp4_invert(&tmp,B);
    Fp4_mul(&t_ans,A,&tmp);
    
    Fp4_set(ANS,&t_ans);
    
    Fp4_clear(&tmp);
    Fp4_clear(&t_ans);
}
void Fp4_pow( Fp4 *ANS, Fp4 *A,mpz_t B){
    int i,length;
    length=(int)mpz_sizeinbase(B,2);
    char binary[length];
    mpz_get_str(binary,2,B);
     Fp4 buf;
    Fp4_init(&buf);
    Fp4_set(&buf,A);
    
    for(i=1; i<length; i++){
        //Fp2_mul(&buf,&buf,&buf);
        Fp4_sqr(&buf,&buf);
        if(binary[i]=='1'){
            Fp4_mul(&buf,A,&buf);
        }
    }
    
    Fp4_set(ANS,&buf);
    Fp4_clear(&buf);
}
void Fp4_sqrt( Fp4 *ANS, Fp4 *A){
    struct Fp4 n,y,x,b,t,tmp_Fp4;
    Fp4_init(&n);
    Fp4_init(&y);
    Fp4_init(&x);
    Fp4_init(&b);
    Fp4_init(&t);
    Fp4_init(&tmp_Fp4);
    Fp4_set(&n,A);
    
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
    
    Fp4_set_random(&n,state);
    while(Fp4_legendre(&n)!=-1){
        Fp4_set_random(&n,state);
    }
    mpz_pow_ui(q,curve_parameters.prime,12);
    mpz_sub_ui(q,q,1);
    mpz_set_ui(e,0);
    while(mpz_odd_p(q)==0){
        mpz_add_ui(e,e,1);
        mpz_div_ui(q,q,2);
    }
    Fp4_pow(&y,&n,q);
    
    mpz_set(r,e);
    
    mpz_sub_ui(tmp_mpz,q,1);
    mpz_div_ui(tmp_mpz,tmp_mpz,2);
    
    Fp4_pow(&x,A,tmp_mpz);
    Fp4_pow(&tmp_Fp4,&x,set_2);
    Fp4_mul(&b,&tmp_Fp4,A);
    Fp4_mul(&x,&x,A);
    
    int m;
    
    while(Fp4_cmp_mpz(&b,set_1)==1){
        m=-1;
        Fp4_set(&tmp_Fp4,&b);
        while(Fp4_cmp_mpz(&tmp_Fp4,set_1)==1){
            m++;
            mpz_pow_ui(tmp_mpz,set_2,m);
            Fp4_pow(&tmp_Fp4,&b,tmp_mpz);
        }
        mpz_sub_ui(tmp_mpz,r,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,curve_parameters.prime);
        // gmp_printf("%Zd,%Zd,%d\n",tmp_mpz,r,m);
        Fp4_pow(&t,&y,tmp_mpz);
        Fp4_pow(&y,&t,set_2);
        // gmp_printf("%Zd,%Zd,\n",y.x0.x0.x0,y.x0.x1.x0);
        mpz_set_ui(r,m);
        Fp4_mul(&x,&x,&t);
        Fp4_mul(&b,&b,&y);
    }
    
    Fp4_set(ANS,&x);
    
    Fp4_clear(&n);
    Fp4_clear(&y);
    Fp4_clear(&x);
    Fp4_clear(&b);
    Fp4_clear(&t);
    Fp4_clear(&tmp_Fp4);
    mpz_clear(tmp_mpz);
    mpz_clear(q);
    mpz_clear(e);
    mpz_clear(r);
    mpz_clear(set_1);
}
int Fp4_legendre( Fp4 *a){
    mpz_t i,cmp;
     Fp4 tmp;
    Fp4_init(&tmp);
    mpz_init(i);
    mpz_init(cmp);
    mpz_set_ui(cmp,1);
    mpz_pow_ui(i,curve_parameters.prime,4);
    mpz_sub_ui(i,i,1);
    mpz_tdiv_q_ui(i,i,2);
    Fp4_pow(&tmp,a,i);
    
    if((Fp4_cmp_mpz(&tmp,cmp))==0){
        Fp4_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return 1;
    }else{
        Fp4_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return -1;
    }
}
int Fp4_cmp( Fp4 *A, Fp4 *B){
    if(Fp2_cmp(&A->x0,&B->x0)==0 && Fp2_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}
int Fp4_cmp_mpz( Fp4 *A,mpz_t B){
    if(Fp2_cmp_mpz(&A->x0,B)==0 && Fp2_cmp_mpz(&A->x1,B)==0){
        return 0;
    }
    return 1;
}
void Fp4_neg( Fp4 *ans, Fp4 *a){
    Fp2_neg(&ans->x0,&a->x0);
    Fp2_neg(&ans->x1,&a->x1);
}

void Fp4_mul_betainv( Fp4 *ANS)
{
     Fp4 tmp;
    Fp4_init(&tmp);
    Fp4_set_ui(&tmp, 0);
    mpz_set(tmp.x1.x1.x0,curve_parameters.curve_a);
    mpz_mul(tmp.x1.x1.x0,tmp.x1.x1.x0,C1_INV);
    
    Fp4_set(ANS, &tmp);
    Fp4_clear(&tmp);
}

