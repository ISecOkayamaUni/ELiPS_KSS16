//
//  bn_fp6.c
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

#include <ELiPS_bn_bls/bn_fp6.h>
#include <ELiPS_bn_bls/bn_bls12_precoms.h>


void Fp6_init(Fp6 *A){
    Fp2_init(&A->x0);
    Fp2_init(&A->x1);
    Fp2_init(&A->x2);
}

void Fp6_clear(Fp6 *A){
    Fp2_clear(&A->x0);
    Fp2_clear(&A->x1);
    Fp2_clear(&A->x2);
}

void Fp6_printf(Fp6 *A,char *str){
    gmp_printf("%s(",str);
    Fp2_printf(&A->x0,"");
    gmp_printf(",");
    Fp2_printf(&A->x1,"");
    gmp_printf(",");
    Fp2_printf(&A->x2,"");
    gmp_printf(")");
}

void Fp6_set(Fp6 *ANS,Fp6 *A){
    Fp2_set(&ANS->x0,&A->x0);
    Fp2_set(&ANS->x1,&A->x1);
    Fp2_set(&ANS->x2,&A->x2);
}

void Fp6_set_ui(Fp6 *ANS,unsigned long int UI){
    Fp2_set_ui(&ANS->x0,UI);
    Fp2_set_ui(&ANS->x1,UI);
    Fp2_set_ui(&ANS->x2,UI);
}

void Fp6_set_mpz(Fp6 *ANS,mpz_t A){
    Fp2_set_mpz(&ANS->x0,A);
    Fp2_set_mpz(&ANS->x1,A);
    Fp2_set_mpz(&ANS->x2,A);
}

void Fp6_set_neg(Fp6 *ANS,Fp6 *A){
    Fp2_set_neg(&ANS->x0,&A->x0);
    Fp2_set_neg(&ANS->x1,&A->x1);
    Fp2_set_neg(&ANS->x2,&A->x2);
}

void Fp6_set_random(Fp6 *ANS,gmp_randstate_t state){
    Fp2_set_random(&ANS->x0,state);
    Fp2_set_random(&ANS->x1,state);
    Fp2_set_random(&ANS->x2,state);
}

void Fp6_mul(Fp6 *ANS,Fp6 *A,Fp6 *B){
    Fp2 tmp00,tmp11,tmp22,tmp,t0,t1,t2;
    Fp2_init(&tmp00);
    Fp2_init(&tmp11);
    Fp2_init(&tmp22);
    Fp2_init(&tmp);
    Fp2_init(&t0);
    Fp2_init(&t1);
    Fp2_init(&t2);
    
    //set
    Fp2_mul(&tmp00,&A->x0,&B->x0);//x0*y0
    Fp2_mul(&tmp11,&A->x1,&B->x1);//x1*y1
    Fp2_mul(&tmp22,&A->x2,&B->x2);//x2*y2
    
    Fp2_add(&t0,&A->x0,&A->x1);//x0+x1
    Fp2_add(&tmp,&B->x0,&B->x1);//y0+y1
    Fp2_mul(&t0,&t0,&tmp);//(x0+x1)(y0+y1)
    
    Fp2_add(&t1,&A->x1,&A->x2);//x1+x2
    Fp2_add(&tmp,&B->x1,&B->x2);//y1+y2
    Fp2_mul(&t1,&t1,&tmp);//(x1+x2)(y1+y2)
    
    Fp2_add(&t2,&B->x0,&B->x2);//y2+y0
    Fp2_add(&tmp,&A->x0,&A->x2);//x2+x0
    Fp2_mul(&t2,&t2,&tmp);//(x2+x0)(y2+y0)
    //x0
    Fp2_sub(&t1,&t1,&tmp11);
    Fp2_sub(&t1,&t1,&tmp22);//(x1+x2)(y1+y2)-x1y1-x2y2
    Fp2_mul_basis(&tmp,&t1);
    Fp2_add(&ANS->x0,&tmp00,&tmp);
    //x1
    Fp2_sub(&t0,&t0,&tmp00);
    Fp2_sub(&t0,&t0,&tmp11);
    Fp2_mul_basis(&tmp,&tmp22);
    Fp2_add(&ANS->x1,&tmp,&t0);
    //x2
    Fp2_sub(&t2,&t2,&tmp00);
    Fp2_sub(&t2,&t2,&tmp22);
    Fp2_add(&ANS->x2,&tmp11,&t2);
    
    //clear
    Fp2_clear(&tmp00);
    Fp2_clear(&tmp11);
    Fp2_clear(&tmp22);
    Fp2_clear(&tmp);
    Fp2_clear(&t0);
    Fp2_clear(&t1);
    Fp2_clear(&t2);
}

void Fp6_mul_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI){
    Fp2_mul_ui(&ANS->x0,&A->x0,UI);
    Fp2_mul_ui(&ANS->x1,&A->x1,UI);
    Fp2_mul_ui(&ANS->x2,&A->x2,UI);
}

void Fp6_mul_mpz(Fp6 *ANS,Fp6 *A,mpz_t B){
    Fp2_mul_mpz(&ANS->x0,&A->x0,B);
    Fp2_mul_mpz(&ANS->x1,&A->x1,B);
    Fp2_mul_mpz(&ANS->x2,&A->x2,B);
}

void Fp6_mul_basis(Fp6 *ANS,Fp6 *A){
    Fp6 tmp;
    Fp6_init(&tmp);
    Fp6_set(&tmp,A);
    
    Fp_sub(&ANS->x0.x0,&tmp.x2.x0,&tmp.x2.x1);
    Fp_add(&ANS->x0.x1,&tmp.x2.x0,&tmp.x2.x1);
    Fp_set(&ANS->x1.x0,&tmp.x0.x0);
    Fp_set(&ANS->x1.x1,&tmp.x0.x1);
    Fp_set(&ANS->x2.x0,&tmp.x1.x0);
    Fp_set(&ANS->x2.x1,&tmp.x1.x1);
    
    Fp6_clear(&tmp);
}

void Fp6_sqr(Fp6 *ANS,Fp6 *A){
    Fp2 tmp00,tmp12_2,tmp01_2,tmp22,tmp;
    Fp2_init(&tmp00);
    Fp2_init(&tmp22);
    Fp2_init(&tmp12_2);
    Fp2_init(&tmp01_2);
    Fp2_init(&tmp);
    
    Fp2_sqr(&tmp00,&A->x0);        //x0^2
    Fp2_sqr(&tmp22,&A->x2);        //x2^2
    Fp2_add(&tmp,&A->x1,&A->x1);        //2x1
    Fp2_mul(&tmp12_2,&tmp,&A->x2);  //2x1x2
    Fp2_mul(&tmp01_2,&A->x0,&tmp);  //2x0x1
    Fp2_add(&tmp,&A->x0,&A->x1);        //x0+x1+x2
    Fp2_add(&tmp,&tmp,&A->x2);
    
    //x0
    Fp2_mul_basis(&ANS->x0,&tmp12_2);
    Fp2_add(&ANS->x0,&ANS->x0,&tmp00);
    //x1
    Fp2_mul_basis(&ANS->x1,&tmp22);
    Fp2_add(&ANS->x1,&ANS->x1,&tmp01_2);
    //x2
    Fp2_sqr(&ANS->x2,&tmp);
    Fp2_add(&tmp,&tmp00,&tmp22);
    Fp2_add(&tmp,&tmp,&tmp12_2);
    Fp2_add(&tmp,&tmp,&tmp01_2);
    Fp2_sub(&ANS->x2,&ANS->x2,&tmp);
    
    Fp2_clear(&tmp00);
    Fp2_clear(&tmp22);
    Fp2_clear(&tmp12_2);
    Fp2_clear(&tmp01_2);
    Fp2_clear(&tmp);
}

void Fp6_add(Fp6 *ANS,Fp6 *A,Fp6 *B){
    Fp2_add(&ANS->x0,&A->x0,&B->x0);
    Fp2_add(&ANS->x1,&A->x1,&B->x1);
    Fp2_add(&ANS->x2,&A->x2,&B->x2);
}

void Fp6_add_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI){
    Fp2_add_ui(&ANS->x0,&A->x0,UI);
    Fp2_add_ui(&ANS->x1,&A->x1,UI);
    Fp2_add_ui(&ANS->x2,&A->x2,UI);
}

void Fp6_add_mpz(Fp6 *ANS,Fp6 *A,mpz_t B){
    Fp2_add_mpz(&ANS->x0,&A->x0,B);
    Fp2_add_mpz(&ANS->x1,&A->x1,B);
    Fp2_add_mpz(&ANS->x2,&A->x2,B);
}

void Fp6_sub(Fp6 *ANS,Fp6 *A,Fp6 *B){
    Fp2_sub(&ANS->x0,&A->x0,&B->x0);
    Fp2_sub(&ANS->x1,&A->x1,&B->x1);
    Fp2_sub(&ANS->x2,&A->x2,&B->x2);
}

void Fp6_sub_ui(Fp6 *ANS,Fp6 *A,unsigned long int UI){
    Fp2_sub_ui(&ANS->x0,&A->x0,UI);
    Fp2_sub_ui(&ANS->x1,&A->x1,UI);
    Fp2_sub_ui(&ANS->x2,&A->x2,UI);
}

void Fp6_sub_mpz(Fp6 *ANS,Fp6 *A,mpz_t B){
    Fp2_sub_mpz(&ANS->x0,&A->x0,B);
    Fp2_sub_mpz(&ANS->x1,&A->x1,B);
    Fp2_sub_mpz(&ANS->x2,&A->x2,B);
}

void Fp6_inv(Fp6 *ANS,Fp6 *A){
    Fp6 frob1,frob2,tmp1,tmp2;
    Fp6_init(&frob1);
    Fp6_init(&frob2);
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    
    Fp2_set(&frob1.x0,&A->x0);
    Fp2_mul_mpz(&frob1.x1,&A->x1,epsilon1);
    Fp2_mul_mpz(&frob1.x2,&A->x2,epsilon2);
    
    Fp2_set(&frob2.x0,&A->x0);
    Fp2_mul_mpz(&frob2.x1,&A->x1,epsilon2);
    Fp2_mul_mpz(&frob2.x2,&A->x2,epsilon1);
    
    Fp6_mul(&tmp1,&frob1,&frob2);
    Fp6_mul(&tmp2,&tmp1,A);
    Fp2_inv(&tmp2.x0,&tmp2.x0);
    Fp2_mul(&ANS->x0,&tmp1.x0,&tmp2.x0);
    Fp2_mul(&ANS->x1,&tmp1.x1,&tmp2.x0);
    Fp2_mul(&ANS->x2,&tmp1.x2,&tmp2.x0);
    
    Fp6_clear(&frob1);
    Fp6_clear(&frob2);
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
}

int  Fp6_legendre(Fp6 *A){
    mpz_t exp;
    mpz_init(exp);
    Fp6 tmp;
    Fp6_init(&tmp);
    
    mpz_pow_ui(exp,curve_parameters.prime,6);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp6_pow(&tmp,A,exp);
    
    if(Fp6_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp6_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp6_clear(&tmp);
        return -1;
    }
}

int  Fp6_isCNR(Fp6 *A){
    Fp6 tmp;
    Fp6_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_pow_ui(exp,curve_parameters.prime,6);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp6_pow(&tmp,A,exp);
    
    if(Fp6_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp6_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp6_clear(&tmp);
        return -1;
    }
}

void Fp6_sqrt(Fp6 *ANS,Fp6 *A){
    Fp6 tmp1,tmp2;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    mpz_t exp,buf;
    mpz_init(exp);
    mpz_init(buf);
    
    Fp2_set(&tmp1.x0,&A->x0);
//    Fp2_mul_mpz(&tmp1.x1,&A->x1,d12_frobenius_constant[f_p4][1].x0.x0);
//    Fp2_mul_mpz(&tmp1.x2,&A->x2,d12_frobenius_constant[f_p4][2].x0.x0);
    
    Fp2_set(&tmp2.x0,&A->x0);
//    Fp2_mul_mpz(&tmp2.x1,&A->x1,d12_frobenius_constant[f_p2][1].x0.x0);
//    Fp2_mul_mpz(&tmp2.x2,&A->x2,d12_frobenius_constant[f_p2][2].x0.x0);
    
    Fp6_mul(&tmp1,&tmp1,&tmp2);
    Fp6_mul(&tmp1,&tmp1,A);
    Fp6_set_ui(&tmp2,0);
    Fp2_sqrt(&tmp2.x0,&tmp1.x0);
    Fp2_inv(&tmp2.x0,&tmp2.x0);
    Fp2_set(&tmp2.x0,&tmp2.x0);
    mpz_pow_ui(exp,curve_parameters.prime,8);
    mpz_pow_ui(buf,curve_parameters.prime,4);
    mpz_add(exp,exp,buf);
    mpz_add_ui(exp,exp,2);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp6_pow(&tmp1,A,exp);
    Fp6_mul(&tmp1,&tmp1,&tmp2);
    Fp6_set(ANS,&tmp1);
    
    mpz_clear(exp);
    mpz_clear(buf);
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
}

void Fp6_pow(Fp6 *ANS,Fp6 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp6 tmp;
    Fp6_init(&tmp);
    Fp6_set(&tmp,A);
    
    for(i=1; i<length; i++){
        Fp6_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp6_mul(&tmp,A,&tmp);
        }
    }
    
    Fp6_set(ANS,&tmp);
    Fp6_clear(&tmp);
}

int  Fp6_cmp(Fp6 *A,Fp6 *B){
    if(Fp2_cmp(&A->x0,&B->x0)==0 && Fp2_cmp(&A->x1,&B->x1)==0 && Fp2_cmp(&A->x2,&B->x2)==0){
        return 0;
    }
    return 1;
}

int  Fp6_cmp_ui(Fp6 *A,unsigned long int UI){
    if(Fp2_cmp_ui(&A->x0,UI)==0 && Fp2_cmp_ui(&A->x1,UI)==0 && Fp2_cmp_ui(&A->x2,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp6_cmp_mpz(Fp6 *A,mpz_t B){
    if(Fp2_cmp_mpz(&A->x0,B)==0 && Fp2_cmp_mpz(&A->x1,B)==0 && Fp2_cmp_mpz(&A->x2,B)==0){
        return 0;
    }
    return 1;
}

int  Fp6_cmp_zero(Fp6 *A){
    if(Fp2_cmp_zero(&A->x0)==0 && Fp2_cmp_zero(&A->x1)==0 && Fp2_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}

int  Fp6_cmp_one(Fp6 *A){
    if(Fp2_cmp_one(&A->x0)==0 && Fp2_cmp_zero(&A->x1)==0 && Fp2_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}

