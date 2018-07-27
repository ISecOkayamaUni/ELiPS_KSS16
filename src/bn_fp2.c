//
//  bn_fp2.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 1/29/18.

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

#include <ELiPS_bn_bls/bn_fp2.h>

//struct Fp2 Fp2_basis_inv;

void Fp2_init(Fp2 *A){
    Fp_init(&A->x0);
    Fp_init(&A->x1);
}

void Fp2_clear(Fp2 *A){
    Fp_clear(&A->x0);
    Fp_clear(&A->x1);
}

void Fp2_printf(Fp2 *A,char *str){
    gmp_printf("%s(",str);
    Fp_printf(&A->x0,"");
    gmp_printf(",");
    Fp_printf(&A->x1,"");
    gmp_printf(")");
}

void Fp2_set(Fp2 *ANS,Fp2 *A){
    Fp_set(&ANS->x0,&A->x0);
    Fp_set(&ANS->x1,&A->x1);
}

void Fp2_set_ui(Fp2 *ANS,unsigned long int UI){
    Fp_set_ui(&ANS->x0,UI);
    Fp_set_ui(&ANS->x1,UI);
}

void Fp2_set_mpz(Fp2 *ANS,mpz_t A){
    Fp_set_mpz(&ANS->x0,A);
    Fp_set_mpz(&ANS->x1,A);
}

void Fp2_set_neg(Fp2 *ANS,Fp2 *A){
    Fp_set_neg(&ANS->x0,&A->x0);
    Fp_set_neg(&ANS->x1,&A->x1);
}

void Fp2_set_random(Fp2 *ANS,gmp_randstate_t state){
    Fp_set_random(&ANS->x0,state);
    Fp_set_random(&ANS->x1,state);
}

void Fp2_mul(Fp2 *ANS,Fp2 *A,Fp2 *B){
    Fp tmp1,tmp2,tmp3,tmp4;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&tmp3);
    Fp_init(&tmp4);
    
    //set
    Fp_mul(&tmp1,&A->x0,&B->x0);//a*c
    Fp_mul(&tmp2,&A->x1,&B->x1);//b*d
    Fp_add(&tmp3,&A->x0,&A->x1);//a+b
    Fp_add(&tmp4,&B->x0,&B->x1);//c+d
    //x0
    Fp_mul_basis(&ANS->x0,&tmp2);//b*d*v
    Fp_add(&ANS->x0,&ANS->x0,&tmp1);//a*c+b*d*v
    //x1
    Fp_mul(&ANS->x1,&tmp3,&tmp4);//(a+b)(c+d)
    Fp_sub(&ANS->x1,&ANS->x1,&tmp1);
    Fp_sub(&ANS->x1,&ANS->x1,&tmp2);
    
    //clear
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&tmp3);
    Fp_clear(&tmp4);
}

void Fp2_mul_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI){
    Fp_mul_ui(&ANS->x0,&A->x0,UI);
    Fp_mul_ui(&ANS->x1,&A->x1,UI);
}

void Fp2_mul_mpz(Fp2 *ANS,Fp2 *A,mpz_t B){
    Fp_mul_mpz(&ANS->x0,&A->x0,B);
    Fp_mul_mpz(&ANS->x1,&A->x1,B);
}

void Fp2_mul_Fp(struct Fp2 *ANS,struct Fp2 *A,struct Fp *B){
    Fp_mul(&ANS->x0,&ANS->x0,B);
    Fp_mul(&ANS->x1,&ANS->x1,B);
}
void Fp2_mul_basis(Fp2 *ANS,Fp2 *A){
    Fp tmp;
    Fp_init(&tmp);
    Fp_set(&tmp,&A->x0);
    
    Fp_sub(&ANS->x0,&tmp,&A->x1);
    Fp_add(&ANS->x1,&tmp,&A->x1);
    
    Fp_clear(&tmp);
}

void Fp2_mul_basis_KSS16(Fp2 *ANS, Fp2 *A){
    Fp_mul_basis_KSS16(&ANS->x0,&A->x1);
    Fp_set(&ANS->x1,&A->x0);
}

void Fp2_mul_KSS16(Fp2 *ANS,Fp2 *A,Fp2 *B){
    struct Fp tmp1,tmp2;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    
    //set
    Fp_mul(&tmp2,&A->x1,&B->x1);//b*d
    Fp_add(&tmp1,&A->x0,&A->x1);//a+b
    Fp_add(&ANS->x1,&B->x0,&B->x1);//c+d
    Fp_mul(&ANS->x1,&tmp1,&ANS->x1);//(a+b)(c+d)
    Fp_mul(&tmp1,&A->x0,&B->x0);//a*c
    //x0
    Fp_mul_basis_KSS16(&ANS->x0,&tmp2);//b*d*v
    Fp_add(&ANS->x0,&ANS->x0,&tmp1);//a*c+b*d*v
    //x1
    Fp_sub(&ANS->x1,&ANS->x1,&tmp1);
    Fp_sub(&ANS->x1,&ANS->x1,&tmp2);
    
    //clear
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
}

void Fp2_sqr_KSS16(Fp2 *ANS,Fp2 *A){
    struct Fp tmp1,tmp2,tmp3;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&tmp3);
    
    Fp_add(&tmp1,&A->x0,&A->x1);
    Fp_mul_basis(&tmp2,&A->x1);
    Fp_add(&tmp2,&tmp2,&A->x0);
    Fp_mul(&tmp3,&A->x0,&A->x1);
    
    //x0
    Fp_mul(&ANS->x0,&tmp1,&tmp2);
    Fp_sub(&ANS->x0,&ANS->x0,&tmp3);
    Fp_mul_basis(&tmp1,&tmp3);
    Fp_sub(&ANS->x0,&ANS->x0,&tmp1);
    //x1
    Fp_add(&ANS->x1,&tmp3,&tmp3);
    
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&tmp3);
}

void Fp2_inv_basis(Fp2 *ANS,Fp2 *A){
    //Fp2_mul(ANS,A,&Fp2_basis_inv);
    Fp2 tmp;
    Fp2_init(&tmp);
    Fp2_set(&tmp,A);
    //TODO
    Fp_add(&ANS->x0,&tmp.x0,&tmp.x1);
    Fp_mul_mpz(&ANS->x0,&ANS->x0,Fp2_basis_inv.x0.x0);
    Fp_sub(&ANS->x1,&tmp.x1,&tmp.x0);
    Fp_mul_mpz(&ANS->x1,&ANS->x1,Fp2_basis_inv.x0.x0);
    
    Fp2_clear(&tmp);
}

void Fp2_sqr(Fp2 *ANS,Fp2 *A){
    Fp tmp1,tmp2;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    
    Fp_add(&tmp1,&A->x0,&A->x1);
    Fp_sub(&tmp2,&A->x0,&A->x1);
    //x1
    Fp_mul(&ANS->x1,&A->x0,&A->x1);
    Fp_add(&ANS->x1,&ANS->x1,&ANS->x1);
    //x0
    Fp_mul(&ANS->x0,&tmp1,&tmp2);
    
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
}

void Fp2_sqr_kss16(struct Fp2 *ANS,struct Fp2 *A){
    struct Fp tmp1,tmp2,tmp3;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&tmp3);
    
    Fp_add(&tmp1,&A->x0,&A->x1);
    Fp_mul_basis_KSS16(&tmp2,&A->x1);
    Fp_add(&tmp2,&tmp2,&A->x0);
    Fp_mul(&tmp3,&A->x0,&A->x1);
    
    //x0
    Fp_mul(&ANS->x0,&tmp1,&tmp2);
    Fp_sub(&ANS->x0,&ANS->x0,&tmp3);
    Fp_mul_basis_KSS16(&tmp1,&tmp3);
    Fp_sub(&ANS->x0,&ANS->x0,&tmp1);
    //x1
    Fp_add(&ANS->x1,&tmp3,&tmp3);
    
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&tmp3);
}


void Fp2_add(Fp2 *ANS,Fp2 *A,Fp2 *B){
    Fp_add(&ANS->x0,&A->x0,&B->x0);
    Fp_add(&ANS->x1,&A->x1,&B->x1);
}

void Fp2_add_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI){
    Fp_add_ui(&ANS->x0,&A->x0,UI);
    Fp_add_ui(&ANS->x1,&A->x1,UI);
}

void Fp2_add_mpz(Fp2 *ANS,Fp2 *A,mpz_t B){
    Fp_add_mpz(&ANS->x0,&A->x0,B);
    Fp_add_mpz(&ANS->x1,&A->x1,B);
}

void Fp2_sub(Fp2 *ANS,Fp2 *A,Fp2 *B){
    Fp_sub(&ANS->x0,&A->x0,&B->x0);
    Fp_sub(&ANS->x1,&A->x1,&B->x1);
}

void Fp2_sub_ui(Fp2 *ANS,Fp2 *A,unsigned long int UI){
    Fp_sub_ui(&ANS->x0,&A->x0,UI);
    Fp_sub_ui(&ANS->x1,&A->x1,UI);
}

void Fp2_sub_mpz(Fp2 *ANS,Fp2 *A,mpz_t B){
    Fp_sub_mpz(&ANS->x0,&A->x0,B);
    Fp_sub_mpz(&ANS->x1,&A->x1,B);
}

void Fp2_inv(Fp2 *ANS,Fp2 *A){
    Fp2 frob,tmp;
    Fp2_init(&frob);
    Fp2_init(&tmp);
    
    Fp_set(&frob.x0,&A->x0);
    Fp_set_neg(&frob.x1,&A->x1);
    
    Fp2_mul(&tmp,A,&frob);
    Fp_inv(&tmp.x0,&tmp.x0);
    Fp2_mul_mpz(ANS,&frob,tmp.x0.x0);
    
    Fp2_clear(&frob);
    Fp2_clear(&tmp);
}

void Fp2_invert_kss16(Fp2 *ANS,Fp2 *A){
    struct Fp2 frob,buf;
    Fp2_init(&frob);
    Fp2_init(&buf);
    
    Fp_set(&frob.x0,&A->x0);
    Fp_neg(&frob.x1,&A->x1);
    Fp2_mul_KSS16(&buf,A,&frob);
    mpz_invert(buf.x0.x0,buf.x0.x0,curve_parameters.prime);
    Fp2_mul_mpz(ANS,&frob,buf.x0.x0);
    
    Fp2_clear(&frob);
    Fp2_clear(&buf);
}

int  Fp2_legendre(Fp2 *A){
    Fp2 tmp;
    Fp2_init(&tmp);
    
    mpz_t exp;
    mpz_init(exp);
    mpz_pow_ui(exp,curve_parameters.prime,2);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp,A,exp);
    
    if(Fp2_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp2_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp2_clear(&tmp);
        return -1;
    }
}

int  Fp2_isCNR(Fp2 *A){
    Fp2 tmp;
    Fp2_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_pow_ui(exp,curve_parameters.prime,2);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&tmp,A,exp);
    
    if(Fp2_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp2_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp2_clear(&tmp);
        return -1;
    }
    
}

void Fp2_sqrt(Fp2 *ANS,Fp2 *A){
    Fp2 x,y,t,k,n,tmp;
    Fp2_init(&x);
    Fp2_init(&y);
    Fp2_init(&t);
    Fp2_init(&k);
    Fp2_init(&n);
    Fp2_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    Fp2_set_random(&n,state);
    while(Fp2_legendre(&n)!=-1){
        Fp2_set_random(&n,state);
    }
    mpz_pow_ui(q,curve_parameters.prime,2);
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp2_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&x,A,exp);
    Fp2_mul(&tmp,&x,&x);
    Fp2_mul(&k,&tmp,A);
    Fp2_mul(&x,&x,A);
    while(Fp2_cmp_one(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp2_pow(&tmp,&k,exp);
        while(Fp2_cmp_one(&tmp)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp2_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp2_pow(&t,&y,result);
        Fp2_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp2_mul(&x,&x,&t);
        Fp2_mul(&k,&k,&y);
    }
    Fp2_set(ANS,&x);
    
    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
    Fp2_clear(&x);
    Fp2_clear(&y);
    Fp2_clear(&t);
    Fp2_clear(&k);
    Fp2_clear(&n);
    Fp2_clear(&tmp);
}

void Fp2_pow(Fp2 *ANS,Fp2 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp2 tmp;
    Fp2_init(&tmp);
    Fp2_set(&tmp,A);
    
    for(i=1; i<length; i++){
        Fp2_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp2_mul(&tmp,A,&tmp);
        }
    }
    
    Fp2_set(ANS,&tmp);
    Fp2_clear(&tmp);
}

int  Fp2_cmp(Fp2 *A,Fp2 *B){
    if(Fp_cmp(&A->x0,&B->x0)==0 && Fp_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}

int  Fp2_cmp_ui(Fp2 *A,unsigned long int UI){
    if(Fp_cmp_ui(&A->x0,UI)==0 && Fp_cmp_ui(&A->x1,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp2_cmp_mpz(Fp2 *A,mpz_t B){
    if(Fp_cmp_mpz(&A->x0,B)==0 && Fp_cmp_mpz(&A->x1,B)==0){
        return 0;
    }
    return 1;
}

int  Fp2_cmp_zero(Fp2 *A){
    if(Fp_cmp_zero(&A->x0)==0 && Fp_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

int  Fp2_cmp_one(Fp2 *A){
    if(Fp_cmp_one(&A->x0)==0 && Fp_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}
void Fp2_neg(Fp2 *ans,Fp2 *a){
    Fp_neg(&ans->x0,&a->x0);
    Fp_neg(&ans->x1,&a->x1);
}

