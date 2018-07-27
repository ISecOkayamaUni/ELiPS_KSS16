//
//  bls12_scm.c
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

#include <ELiPS_bn_bls/bls12_scm.h>

void bls12_plain_G1_scm(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    gettimeofday(&t0,NULL);
    
    EFp Tmp_P;
    EFp_init(&Tmp_P);
    
    bls12_EFp12_to_EFp(&Tmp_P,P);
    EFp_SCM(&Tmp_P,&Tmp_P,scalar);
    bls12_EFp_to_EFp12(ANS,&Tmp_P);
    
    EFp_clear(&Tmp_P);
    
    gettimeofday(&t1,NULL);
    bls12_G1SCM_PLAIN=timedifference_msec(t0,t1);
}

void bls12_2split_G1_scm(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    gettimeofday(&t0,NULL);
    
    int i,length_s[2],loop_length;
    EFp12 Buf;
    EFp12_init(&Buf);
    EFp next_P,tmp_P,skew_P_2;
    EFp_init(&next_P);
    EFp_init(&tmp_P);
    EFp_init(&skew_P_2);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    EFp table[4];
    for(i=0; i<4; i++){
        EFp_init(&table[i]);
    }
    
    //set
    bls12_EFp12_to_EFp(&tmp_P,P);
    bls12_EFp_skew_frobenius_map_p2(&skew_P_2,&tmp_P);
    //set table
    table[0].infinity=1;    //00
    EFp_set(&table[1],&tmp_P);    //01
    EFp_set(&table[2],&skew_P_2);    //10
    EFp_ECA(&table[3],&tmp_P,&skew_P_2);    //11
    
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
    EFp_set(&next_P,&table[binary[0]]);
    
    //SCM
    for(i=1; i<loop_length; i++){
        EFp_ECD(&next_P,&next_P);
        EFp_ECA(&next_P,&next_P,&table[binary[i]]);
    }
    
    bls12_EFp_to_EFp12(ANS,&next_P);
    
    
    EFp12_clear(&Buf);
    EFp_clear(&next_P);
    EFp_clear(&tmp_P);
    EFp_clear(&skew_P_2);
    mpz_clear(buf);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
    //table
    for(i=0; i<4; i++){
        EFp_clear(&table[i]);
    }
    
    gettimeofday(&t1,NULL);
    bls12_G1SCM_2SPLIT=timedifference_msec(t0,t1);
}

void bls12_plain_G2_scm(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
    gettimeofday(&t0,NULL);
    
    EFp2 twisted_Q;
    EFp2_init(&twisted_Q);
    
    bls12_EFp12_to_EFp2(&twisted_Q,Q);
    EFp2_SCM(&twisted_Q,&twisted_Q,scalar);
    bls12_EFp2_to_EFp12(ANS,&twisted_Q);
    
    EFp2_clear(&twisted_Q);
    
    gettimeofday(&t1,NULL);
    bls12_G2SCM_PLAIN=timedifference_msec(t0,t1);
}

void bls12_2split_G2_scm(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
    gettimeofday(&t0,NULL);
    
    int i,length_s[2],loop_length;
    EFp2 next_twisted_Q,twisted_Q,twisted_Q_2x;
    EFp2_init(&next_twisted_Q);
    EFp2_init(&twisted_Q);
    EFp2_init(&twisted_Q_2x);
    mpz_t s[2],buf;
    mpz_init(buf);
    for(i=0; i<2; i++){
        mpz_init(s[i]);
    }
    //table
    EFp2 table[4];
    for(i=0; i<4; i++){
        EFp2_init(&table[i]);
    }
    
    //set
    bls12_EFp12_to_EFp2(&twisted_Q,Q);                //twisted_Q
    bls12_EFp2_skew_frobenius_map_p2(&twisted_Q_2x,&twisted_Q);//twisted_Q_2x
    //set table
    table[0].infinity=1;                        //00
    EFp2_set(&table[1],&twisted_Q);            //01
    EFp2_set(&table[2],&twisted_Q_2x);            //10
    EFp2_ECA(&table[3],&twisted_Q,&twisted_Q_2x);    //11
    
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
    EFp2_set(&next_twisted_Q,&table[binary[0]]);
    
    //SCM
    for(i=1; i<loop_length; i++){
        EFp2_ECD(&next_twisted_Q,&next_twisted_Q);
        EFp2_ECA(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
    }
    
    bls12_EFp2_to_EFp12(ANS,&next_twisted_Q);
    
    mpz_clear(buf);
    EFp2_clear(&next_twisted_Q);
    EFp2_clear(&twisted_Q);
    EFp2_clear(&twisted_Q_2x);
    for(i=0; i<2; i++){
        mpz_clear(s[i]);
    }
    for(i=0; i<4; i++){
        EFp2_clear(&table[i]);
    }
    
    gettimeofday(&t1,NULL);
    bls12_G2SCM_2SPLIT=timedifference_msec(t0,t1);
}

void bls12_4split_G2_scm(EFp12 *ANS,EFp12 *Q,mpz_t scalar){
    gettimeofday(&t0,NULL);
    
    int i,length_s[4],loop_length;
    EFp2 next_twisted_Q,twisted_Q,twisted_Q_x,twisted_Q_2x,twisted_Q_3x;
    EFp2_init(&next_twisted_Q);
    EFp2_init(&twisted_Q);
    EFp2_init(&twisted_Q_x);
    EFp2_init(&twisted_Q_2x);
    EFp2_init(&twisted_Q_3x);
    mpz_t A,B,s[4],x_2,x_1;
    mpz_init(A);
    mpz_init(B);
    mpz_init(x_1);
    mpz_init(x_2);
    for(i=0; i<4; i++){
        mpz_init(s[i]);
    }
    //table
    EFp2 table[16];
    for(i=0; i<16; i++){
        EFp2_init(&table[i]);
    }
    
    //set twisted_Q
    bls12_EFp12_to_EFp2(&twisted_Q,Q);                    //twisted_Q
    bls12_EFp2_skew_frobenius_map_p1(&twisted_Q_x,&twisted_Q);        //twisted_Q_x
    EFp2_set_neg(&twisted_Q_x,&twisted_Q_x);
    bls12_EFp2_skew_frobenius_map_p2(&twisted_Q_2x,&twisted_Q);    //twisted_Q_2x
    bls12_EFp2_skew_frobenius_map_p3(&twisted_Q_3x,&twisted_Q);    //twisted_Q_3x
    EFp2_set_neg(&twisted_Q_3x,&twisted_Q_3x);
    
    //set table
    table[0].infinity=1;                            //0000
    EFp2_set(&table[1],&twisted_Q);                //0001
    EFp2_set(&table[2],&twisted_Q_x);                //0010
    EFp2_ECA(&table[3],&twisted_Q_x,&twisted_Q);        //0011
    EFp2_set(&table[4],&twisted_Q_2x);                //0100
    EFp2_ECA(&table[5],&twisted_Q_2x,&twisted_Q);        //0101
    EFp2_ECA(&table[6],&twisted_Q_2x,&twisted_Q_x);        //0110
    EFp2_ECA(&table[7],&table[6],&twisted_Q);            //0111
    EFp2_set(&table[8],&twisted_Q_3x);                //1000
    EFp2_ECA(&table[9],&twisted_Q_3x,&twisted_Q);        //1001
    EFp2_ECA(&table[10],&twisted_Q_3x,&twisted_Q_x);    //1010
    EFp2_ECA(&table[11],&twisted_Q_3x,&table[3]);        //1011
    EFp2_ECA(&table[12],&twisted_Q_3x,&twisted_Q_2x);    //1100
    EFp2_ECA(&table[13],&table[12],&twisted_Q);        //1101
    EFp2_ECA(&table[14],&table[12],&twisted_Q_x);        //1110
    EFp2_ECA(&table[15],&table[14],&twisted_Q);        //1111
    
    //set
    //s0,s1,s2,s3
    mpz_neg(x_1,bls12_X);
    mpz_mul(x_2,x_1,x_1);
    mpz_tdiv_qr(B,A,scalar,x_2);
    mpz_tdiv_qr(s[1],s[0],A,x_1);
    mpz_tdiv_qr(s[3],s[2],B,x_1);
    
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
    
    EFp2_set(&next_twisted_Q,&table[binary[0]]);
    
    //SCM
    for(i=1; i<loop_length; i++){
        EFp2_ECD(&next_twisted_Q,&next_twisted_Q);
        EFp2_ECA(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
    }
    
    bls12_EFp2_to_EFp12(ANS,&next_twisted_Q);
    
    EFp2_clear(&next_twisted_Q);
    EFp2_clear(&twisted_Q);
    EFp2_clear(&twisted_Q_x);
    EFp2_clear(&twisted_Q_2x);
    EFp2_init(&twisted_Q_3x);
    mpz_clear(x_1);
    mpz_clear(x_2);
    for(i=0; i<4; i++){
        mpz_clear(s[i]);
    }
    for(i=0; i<16; i++){
        EFp2_clear(&table[i]);
    }
    
    gettimeofday(&t1,NULL);
    bls12_G2SCM_4SPLIT=timedifference_msec(t0,t1);
}


