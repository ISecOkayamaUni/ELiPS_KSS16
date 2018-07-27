//
//  bn_efp12.c
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
#include <ELiPS_bn_bls/bn_efp12.h>

void EFp12_init(EFp12 *P){
    Fp12_init(&P->x);
    Fp12_init(&P->y);
    P->infinity=0;
}

void EFp12_set(EFp12 *ANS,EFp12 *A){
    Fp12_set(&ANS->x,&A->x);
    Fp12_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp12_set_ui(EFp12 *ANS,unsigned long int UI){
    Fp12_set_ui(&ANS->x,UI);
    Fp12_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp12_set_mpz(EFp12 *ANS,mpz_t A){
    Fp12_set_mpz(&ANS->x,A);
    Fp12_set_mpz(&ANS->y,A);
    ANS->infinity=0;
}

void EFp12_set_neg(EFp12 *ANS,EFp12 *A){
    Fp12_set(&ANS->x,&A->x);
    Fp12_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp12_clear(EFp12 *P){
    Fp12_clear(&P->x);
    Fp12_clear(&P->y);
}

void EFp12_printf(EFp12 *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp12_printf(&P->x,"");
        printf(",");
        Fp12_printf(&P->y,"");
        printf(")");
    }else{
        printf("0");
    }
}

void EFp12_rational_point_bn(EFp12 *P){
    Fp12 tmp1,tmp2;
    Fp12_init(&tmp1);
    Fp12_init(&tmp2);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp12_set_random(&P->x,state);
        Fp12_sqr(&tmp1,&P->x);
        Fp12_mul(&tmp2,&tmp1,&P->x);
        mpz_sub(tmp2.x0.x0.x0.x0,tmp2.x0.x0.x0.x0,curve_parameters.curve_b);
        if(Fp12_legendre(&tmp2)==1){
            Fp12_sqrt(&P->y,&tmp2);
            break;
        }
    }
    
    Fp12_clear(&tmp1);
    Fp12_clear(&tmp2);
}

void EFp12_rational_point_bls12(EFp12 *P){
    Fp12 tmp1,tmp2;
    Fp12_init(&tmp1);
    Fp12_init(&tmp2);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp12_set_random(&P->x,state);
        Fp12_sqr(&tmp1,&P->x);
        Fp12_mul(&tmp2,&tmp1,&P->x);
        mpz_add(tmp2.x0.x0.x0.x0,tmp2.x0.x0.x0.x0, curve_parameters.curve_b);
        if(Fp12_legendre(&tmp2)==1){
            Fp12_sqrt(&P->y,&tmp2);
            break;
        }
    }
    
    Fp12_clear(&tmp1);
    Fp12_clear(&tmp2);
}

void EFp12_ECD(EFp12 *ANS,EFp12 *P){
    if(Fp12_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp12 Tmp_P;
    EFp12_init(&Tmp_P);
    EFp12_set(&Tmp_P,P);
    Fp12 tmp1,tmp2,lambda;
    Fp12_init(&tmp1);
    Fp12_init(&tmp2);
    Fp12_init(&lambda);
    
    Fp12_add(&tmp1,&Tmp_P.y,&Tmp_P.y);
    
    Fp12_inv(&tmp1,&tmp1);
    Fp12_sqr(&tmp2,&Tmp_P.x);
    Fp12_add(&lambda,&tmp2,&tmp2);
    Fp12_add(&tmp2,&tmp2,&lambda);
    Fp12_mul(&lambda,&tmp1,&tmp2);
    
    Fp12_sqr(&tmp1,&lambda);
    Fp12_add(&tmp2,&Tmp_P.x,&Tmp_P.x);
    Fp12_sub(&ANS->x,&tmp1,&tmp2);
    
    Fp12_sub(&tmp1,&Tmp_P.x,&ANS->x);
    Fp12_mul(&tmp2,&lambda,&tmp1);
    Fp12_sub(&ANS->y,&tmp2,&Tmp_P.y);
    
    Fp12_clear(&tmp1);
    Fp12_clear(&tmp2);
    Fp12_clear(&lambda);
    EFp12_clear(&Tmp_P);
}

void EFp12_ECA(EFp12 *ANS,EFp12 *P1,EFp12 *P2){
    if(P1->infinity==1){
        EFp12_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp12_set(ANS,P1);
        return;
    }else if(Fp12_cmp(&P1->x,&P2->x)==0){
        if(Fp12_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp12_ECD(ANS,P1);
            return;
        }
    }
    
    EFp12 Tmp_P1,Tmp_P2;
    EFp12_init(&Tmp_P1);
    EFp12_set(&Tmp_P1,P1);
    EFp12_init(&Tmp_P2);
    EFp12_set(&Tmp_P2,P2);
    Fp12 tmp1,tmp2,lambda;
    Fp12_init(&tmp1);
    Fp12_init(&tmp2);
    Fp12_init(&lambda);
    
    Fp12_sub(&tmp1,&Tmp_P2.x,&Tmp_P1.x);
    Fp12_inv(&tmp1,&tmp1);
    Fp12_sub(&tmp2,&Tmp_P2.y,&Tmp_P1.y);
    Fp12_mul(&lambda,&tmp1,&tmp2);
    Fp12_sqr(&tmp1,&lambda);
    Fp12_sub(&tmp2,&tmp1,&Tmp_P1.x);
    Fp12_sub(&ANS->x,&tmp2,&Tmp_P2.x);
    Fp12_sub(&tmp1,&Tmp_P1.x,&ANS->x);
    Fp12_mul(&tmp2,&lambda,&tmp1);
    Fp12_sub(&ANS->y,&tmp2,&Tmp_P1.y);
        
    //clear 
    Fp12_clear(&tmp1);
    Fp12_clear(&tmp2);
    Fp12_clear(&lambda);
    EFp12_clear(&Tmp_P1);
    EFp12_clear(&Tmp_P2);
}

void EFp12_SCM(EFp12 *ANS,EFp12 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp12_set(ANS,P);
        return;
    }
    
    EFp12 Tmp_P,Next_P;
    EFp12_init(&Tmp_P);
    EFp12_set(&Tmp_P,P);
    EFp12_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp12_set(&Next_P,&Tmp_P);
    for(i=1; i<length; i++){
        EFp12_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp12_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp12_set(ANS,&Next_P);
    
    EFp12_clear(&Next_P);
    EFp12_clear(&Tmp_P);
}
