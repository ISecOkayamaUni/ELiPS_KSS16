//
//  bn_efp6.c
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
#include <ELiPS_bn_bls/bn_efp6.h>

void EFp6_init(EFp6 *P){
    Fp6_init(&P->x);
    Fp6_init(&P->y);
    P->infinity=0;
}

void EFp6_clear(EFp6 *P){
    Fp6_clear(&P->x);
    Fp6_clear(&P->y);
}

void EFp6_printf(EFp6 *P,char *str){
    gmp_printf("%s",str);
    if(P->infinity==0){
        gmp_printf("(");
        Fp6_printf(&P->x,"");
        gmp_printf(",");
        Fp6_printf(&P->y,"");
        gmp_printf(")");
    }else{
        gmp_printf("infinity");
    }
}

void EFp6_set(EFp6 *ANS,EFp6 *P){
    Fp6_set(&ANS->x,&P->x);
    Fp6_set(&ANS->y,&P->y);
    ANS->infinity=P->infinity;
}

void EFp6_set_ui(EFp6 *ANS,unsigned long int UI){
    Fp6_set_ui(&ANS->x,UI);
    Fp6_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp6_set_mpz(EFp6 *ANS,mpz_t A){
    Fp6_set_mpz(&ANS->x,A);
    Fp6_set_mpz(&ANS->y,A);
    ANS->infinity=0;
}

void EFp6_set_neg(EFp6 *ANS,EFp6 *P){
    Fp6_set(&ANS->x,&P->x);
    Fp6_set_neg(&ANS->y,&P->y);
    ANS->infinity=P->infinity;
}

void EFp6_rational_point(EFp6 *P){
    Fp6 tmp1,tmp2;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp6_set_random(&P->x,state);
        Fp6_sqr(&tmp1,&P->x);
        Fp6_mul(&tmp2,&tmp1,&P->x);
        mpz_sub(tmp2.x0.x0.x0,tmp2.x0.x0.x0,curve_parameters.curve_b);
        if(Fp6_legendre(&tmp2)==1){
            Fp6_sqrt(&P->y,&tmp2);
            break;
        }
    }
    
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
}

void EFp6_ECD(EFp6 *ANS,EFp6 *P){
    if(Fp6_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp6 Tmp_P;
    EFp6_init(&Tmp_P);
    EFp6_set(&Tmp_P,P);
    Fp6 tmp1,tmp2,lambda;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    Fp6_init(&lambda);
    
    Fp6_mul_ui(&tmp1,&Tmp_P.y,2);
    
    Fp6_inv(&tmp1,&tmp1);
    Fp6_mul(&tmp2,&Tmp_P.x,&Tmp_P.x);
    Fp6_mul_ui(&tmp2,&tmp2,3);
    Fp6_mul(&lambda,&tmp1,&tmp2);
    Fp6_sqr(&tmp1,&lambda);
    Fp6_mul_ui(&tmp2,&Tmp_P.x,2);
    Fp6_sub(&ANS->x,&tmp1,&tmp2);
    Fp6_sub(&tmp1,&Tmp_P.x,&ANS->x);
    Fp6_mul(&tmp2,&lambda,&tmp1);
    Fp6_sub(&ANS->y,&tmp2,&Tmp_P.y);
    
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
    Fp6_clear(&lambda);
    EFp6_clear(&Tmp_P);
}

void EFp6_ECA(EFp6 *ANS,EFp6 *P1,EFp6 *P2){
    if(P1->infinity==1){
        EFp6_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp6_set(ANS,P1);
        return;
    }else if(Fp6_cmp(&P1->x,&P2->x)==0){
        if(Fp6_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp6_ECD(ANS,P1);
            return;
        }
    }
    
    EFp6 Tmp_P1,Tmp_P2;
    EFp6_init(&Tmp_P1);
    EFp6_set(&Tmp_P1,P1);
    EFp6_init(&Tmp_P2);
    EFp6_set(&Tmp_P2,P2);
    Fp6 tmp1,tmp2,lambda;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    Fp6_init(&lambda);
    
    Fp6_sub(&tmp1,&Tmp_P2.x,&Tmp_P1.x);
    Fp6_inv(&tmp1,&tmp1);
    Fp6_sub(&tmp2,&Tmp_P2.y,&Tmp_P1.y);
    Fp6_mul(&lambda,&tmp1,&tmp2);
    Fp6_sqr(&tmp1,&lambda);
    Fp6_sub(&tmp2,&tmp1,&Tmp_P1.x);
    Fp6_sub(&ANS->x,&tmp2,&Tmp_P2.x);
    Fp6_sub(&tmp1,&Tmp_P1.x,&ANS->x);
    Fp6_mul(&tmp2,&lambda,&tmp1);
    Fp6_sub(&ANS->y,&tmp2,&Tmp_P1.y);
    
    //clear
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
    Fp6_clear(&lambda);
    EFp6_clear(&Tmp_P1);
    EFp6_clear(&Tmp_P2);
}

void EFp6_SCM(EFp6 *ANS,EFp6 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp6_set(ANS,P);
        return;
    }
    
    EFp6 Tmp_P,Next_P;
    EFp6_init(&Tmp_P);
    EFp6_set(&Tmp_P,P);
    EFp6_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp6_set(&Next_P,&Tmp_P);
    for(i=1; i<length; i++){
        EFp6_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp6_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    
    EFp6_set(ANS,&Next_P);
    
    EFp6_clear(&Next_P);
    EFp6_clear(&Tmp_P);
}
