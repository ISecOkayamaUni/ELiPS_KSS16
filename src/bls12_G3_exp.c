//
//  bls12_G3_exp.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/6/18.

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

#include <ELiPS_bn_bls/bls12_G3_exp.h>

void bls12_plain_G3_exp(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    gettimeofday(&t0,NULL);
    
    Fp12_pow(ANS,A,scalar);
    
    gettimeofday(&t1,NULL);
    bls12_G3EXP_PLAIN=timedifference_msec(t0,t1);
}

void bls12_2split_G3_exp(Fp12 *ANS,Fp12 *A,mpz_t scalar){
    gettimeofday(&t0,NULL);
    
    int i,length_s[2],loop_length;
    Fp12 Buf;
    Fp12_init(&Buf);
    Fp12 next_f,f,frobenius_f_2x;
    Fp12_init(&next_f);
    Fp12_init(&f);
    Fp12_init(&frobenius_f_2x);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    Fp12 table[4];
    for(i=0; i<4; i++){
        Fp12_init(&table[i]);
    }
    
    //set
    Fp12_set(&f, A);
    bls12_Fp12_frobenius_map_p2(&frobenius_f_2x, A);
    
    //set table
    Fp_set_ui(&table[0].x0.x0.x0,1);    //00
    Fp12_set(&table[1],&f);                //01
    Fp12_set(&table[2],&frobenius_f_2x);    //10
    Fp12_mul(&table[3],&table[1],&table[2]);    //11
    
    //s0,s1
    mpz_neg(buf,bls12_X);
    mpz_pow_ui(buf,buf,2);
    mpz_tdiv_qr(s[1],s[0],scalar,buf);
    //binary
    loop_length=0;
    for(i=0; i<2; i++){
        length_s[i]=(int)mpz_sizeinbase(s[i],2);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }
    }
    //set binary
    char binary_s[2][loop_length+1];
    char str[5],*e;
    int binary[loop_length+1];
    for(i=0; i<2; i++){
        if(length_s[i]==loop_length){
            mpz_get_str(binary_s[i],2,s[i]);
        }else{
            char binary_buf[loop_length+1];
            mpz_get_str(binary_buf,2,s[i]);
            memset(binary_s[i],'0',sizeof(binary_s[i]));
            memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
        }
    }
    for(i=0; i<loop_length; i++){
        sprintf(str,"%c%c",binary_s[1][i],binary_s[0][i]);
        binary[i]=(int)strtol(str,&e,2);
    }
    Fp12_set(&next_f,&table[binary[0]]);
    
    //SCM
    for(i=1; i<loop_length; i++){
        Fp12_mul(&next_f,&next_f,&next_f);
        Fp12_mul(&next_f,&next_f,&table[binary[i]]);
    }
    
    Fp12_set(ANS,&next_f);
    
    mpz_clear(buf);
    Fp12_clear(&next_f);
    Fp12_clear(&f);
    Fp12_clear(&frobenius_f_2x);
    Fp12_clear(&Buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
    for(i=0; i<4; i++){
        Fp12_clear(&table[i]);
    }
    
    gettimeofday(&t1,NULL);
    bls12_G3EXP_2SPLIT=timedifference_msec(t0,t1);
}

void bls12_4split_G3_exp( Fp12 *ANS,Fp12 *A,mpz_t scalar){
    gettimeofday(&t0,NULL);
    
    int i,length_s[4],loop_length;
    Fp12 Buf;
    Fp12_init(&Buf);
    Fp12 next_f,f,frobenius_f_x,frobenius_f_2x,frobenius_f_3x;
    Fp12_init(&next_f);
    Fp12_init(&f);
    Fp12_init(&frobenius_f_x);
    Fp12_init(&frobenius_f_2x);
    Fp12_init(&frobenius_f_3x);
    mpz_t C,D,s[4],x_2,x_1;
    mpz_init(C);
    mpz_init(D);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    mpz_init(x_1);
    mpz_init(x_2);
    //table
    Fp12 table[16];
    for(i=0; i<16; i++){
        Fp12_init(&table[i]);
    }
    
    //set
    Fp12_set(&f, A);
    bls12_Fp12_frobenius_map_p1(&frobenius_f_x, A);
    Fp12_inv(&frobenius_f_x,&frobenius_f_x);
    bls12_Fp12_frobenius_map_p2(&frobenius_f_2x, A);
    bls12_Fp12_frobenius_map_p3(&frobenius_f_3x, A);
    Fp12_inv(&frobenius_f_3x,&frobenius_f_3x);
    //set table
    Fp_set_ui(&table[0].x0.x0.x0,1);            //0000
    Fp12_set(&table[1],&f);                    //0001
    Fp12_set(&table[2],&frobenius_f_x);        //0010
    Fp12_mul(&table[3],&table[1],&table[2]);    //0011
    Fp12_set(&table[4],&frobenius_f_2x);        //0100
    Fp12_mul(&table[5],&table[4],&table[1]);    //0101
    Fp12_mul(&table[6],&table[2],&table[4]);    //0110
    Fp12_mul(&table[7],&table[6],&table[1]);    //0111
    Fp12_set(&table[8],&frobenius_f_3x);        //1000
    Fp12_mul(&table[9],&table[8],&table[1]);    //1001
    Fp12_mul(&table[10],&table[8],&table[2]);    //1010
    Fp12_mul(&table[11],&table[10],&table[1]);    //1011
    Fp12_mul(&table[12],&table[8],&table[4]);    //1100
    Fp12_mul(&table[13],&table[12],&table[1]);    //1101
    Fp12_mul(&table[14],&table[12],&table[2]);    //1110
    Fp12_mul(&table[15],&table[14],&table[1]);    //1111
    //set
    //s0,s1,s2,s3
    mpz_neg(x_1,bls12_X);
    mpz_mul(x_2,x_1,x_1);
    mpz_tdiv_qr(D,C,scalar,x_2);
    mpz_tdiv_qr(s[1],s[0],C,x_1);
    mpz_tdiv_qr(s[3],s[2],D,x_1);
    
    //binary
    loop_length=0;
    for(i=0; i<4; i++){
        length_s[i]=(int)mpz_sizeinbase(s[i],2);
        if(loop_length<length_s[i]){
            loop_length=length_s[i];
        }
    }
    //set binary
    char binary_s[4][loop_length+1];
    char str[5],*e;
    int binary[loop_length+1];
    for(i=0; i<4; i++){
        if(length_s[i]==loop_length){
            mpz_get_str(binary_s[i],2,s[i]);
        }else{
            char binary_buf[loop_length+1];
            mpz_get_str(binary_buf,2,s[i]);
            memset(binary_s[i],'0',sizeof(binary_s[i]));
            memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
        }
    }
    for(i=0; i<loop_length; i++){
        sprintf(str,"%c%c%c%c",binary_s[3][i],binary_s[2][i],binary_s[1][i],binary_s[0][i]);
        binary[i]=(int)strtol(str,&e,2);
    }
    
    Fp12_set(&next_f,&table[binary[0]]);
    
    //SCM
    for(i=1; i<loop_length; i++){
        Fp12_mul(&next_f,&next_f,&next_f);
        Fp12_mul(&next_f,&next_f,&table[binary[i]]);
    }
    
    Fp12_set(ANS,&next_f);
    
    Fp12_clear(&Buf);
    Fp12_clear(&next_f);
    Fp12_clear(&f);
    Fp12_clear(&frobenius_f_x);
    Fp12_clear(&frobenius_f_2x);
    Fp12_clear(&frobenius_f_3x);
    mpz_clear(x_1);
    mpz_clear(x_2);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
    for(i=0; i<16; i++){
        Fp12_clear(&table[i]);
    }
    
    gettimeofday(&t1,NULL);
    bls12_G3EXP_4SPLIT=timedifference_msec(t0,t1);
}
