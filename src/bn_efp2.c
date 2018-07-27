//
//  bn_efp2.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/1/18.

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

#include <ELiPS_bn_bls/bn_efp2.h>


void EFp2_init(EFp2 *P){
    Fp2_init(&P->x);
    Fp2_init(&P->y);
    P->infinity=0;
}

void EFp2_clear(EFp2 *P){
    Fp2_clear(&P->x);
    Fp2_clear(&P->y);
}

void EFp2_printf(EFp2 *P,char *str){
    gmp_printf("%s",str);
    if(P->infinity==0){
        gmp_printf("(");
        Fp2_printf(&P->x,"");
        gmp_printf(",");
        Fp2_printf(&P->y,"");
        gmp_printf(")");
    }else{
        gmp_printf("infinity");
    }
}

void EFp2_set(EFp2 *ANS,EFp2 *P){
    Fp2_set(&ANS->x,&P->x);
    Fp2_set(&ANS->y,&P->y);
    ANS->infinity=P->infinity;
}

void EFp2_set_ui(EFp2 *ANS,unsigned long int UI){
    Fp2_set_ui(&ANS->x,UI);
    Fp2_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp2_set_mpz(EFp2 *ANS,mpz_t A){
    Fp2_set_mpz(&ANS->x,A);
    Fp2_set_mpz(&ANS->y,A);
    ANS->infinity=0;
}

void EFp2_set_neg(EFp2 *ANS,EFp2 *P){
    Fp2_set(&ANS->x,&P->x);
    Fp2_set_neg(&ANS->y,&P->y);
    ANS->infinity=P->infinity;
}

void EFp2_rational_point(EFp2 *P){
    Fp2 tmp1,tmp2;
    Fp2_init(&tmp1);
    Fp2_init(&tmp2);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp2_set_random(&P->x,state);
        Fp2_sqr(&tmp1,&P->x);
        Fp2_mul(&tmp2,&tmp1,&P->x);
        mpz_sub(tmp2.x0.x0,tmp2.x0.x0,curve_parameters.curve_b);
        if(Fp2_legendre(&tmp2)==1){
            Fp2_sqrt(&P->y,&tmp2);
            break;
        }
    }
    
    Fp2_clear(&tmp1);
    Fp2_clear(&tmp2);
}

void EFp2_ECD(EFp2 *ANS,EFp2 *P){
    if(Fp2_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp2 Tmp_P;
    EFp2_init(&Tmp_P);
    EFp2_set(&Tmp_P,P);
    Fp2 tmp1,tmp2,lambda;
    Fp2_init(&tmp1);
    Fp2_init(&tmp2);
    Fp2_init(&lambda);
    
    Fp2_mul_ui(&tmp1,&Tmp_P.y,2);
    
    Fp2_inv(&tmp1,&tmp1);
    Fp2_mul(&tmp2,&Tmp_P.x,&Tmp_P.x);
    Fp2_mul_ui(&tmp2,&tmp2,3);
    Fp2_mul(&lambda,&tmp1,&tmp2);
    
    Fp2_sqr(&tmp1,&lambda);
    Fp2_mul_ui(&tmp2,&Tmp_P.x,2);
    Fp2_sub(&ANS->x,&tmp1,&tmp2);
    
    Fp2_sub(&tmp1,&Tmp_P.x,&ANS->x);
    Fp2_mul(&tmp2,&lambda,&tmp1);
    Fp2_sub(&ANS->y,&tmp2,&Tmp_P.y);
    
    Fp2_clear(&tmp1);
    Fp2_clear(&tmp2);
    Fp2_clear(&lambda);
    EFp2_clear(&Tmp_P);
}

void EFp2_ECA(EFp2 *ANS,EFp2 *P1,EFp2 *P2){
    if(P1->infinity==1){
        EFp2_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp2_set(ANS,P1);
        return;
    }else if(Fp2_cmp(&P1->x,&P2->x)==0){
        if(Fp2_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp2_ECD(ANS,P1);
            return;
        }
    }
    
    EFp2 Tmp_P1,Tmp_P2;
    EFp2_init(&Tmp_P1);
    EFp2_set(&Tmp_P1,P1);
    EFp2_init(&Tmp_P2);
    EFp2_set(&Tmp_P2,P2);
    Fp2 tmp1,tmp2,lambda;
    Fp2_init(&tmp1);
    Fp2_init(&tmp2);
    Fp2_init(&lambda);
    
    Fp2_sub(&tmp1,&Tmp_P2.x,&Tmp_P1.x);
    Fp2_inv(&tmp1,&tmp1);
    Fp2_sub(&tmp2,&Tmp_P2.y,&Tmp_P1.y);
    Fp2_mul(&lambda,&tmp1,&tmp2);
    Fp2_sqr(&tmp1,&lambda);
    Fp2_sub(&tmp2,&tmp1,&Tmp_P1.x);
    Fp2_sub(&ANS->x,&tmp2,&Tmp_P2.x);
    Fp2_sub(&tmp1,&Tmp_P1.x,&ANS->x);
    Fp2_mul(&tmp2,&lambda,&tmp1);
    Fp2_sub(&ANS->y,&tmp2,&Tmp_P1.y);
    
    //clear
    Fp2_clear(&tmp1);
    Fp2_clear(&tmp2);
    Fp2_clear(&lambda);
    EFp2_clear(&Tmp_P1);
    EFp2_clear(&Tmp_P2);
}

void EFp2_SCM(EFp2 *ANS,EFp2 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp2_set(ANS,P);
        return;
    }
    
    EFp2 Tmp_P,Next_P;
    EFp2_init(&Tmp_P);
    EFp2_set(&Tmp_P,P);
    EFp2_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp2_set(&Next_P,&Tmp_P);
    for(i=1; i<length; i++){
        EFp2_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp2_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp2_set(ANS,&Next_P);
    
    EFp2_clear(&Next_P);
    EFp2_clear(&Tmp_P);
}
