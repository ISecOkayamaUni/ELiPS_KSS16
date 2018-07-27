//
//  bls12_line_tate.c
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

#include <ELiPS_bn_bls/bls12_line_tate.h>

static void bls12_EFp_ECD_return_lambda(EFp *ANS,Fp *lambda,EFp *P){
    if(Fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp Tmp_P;
    EFp_init(&Tmp_P);
    EFp_set(&Tmp_P,P);
    Fp tmp1,tmp2;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    
    Fp_mul_ui(&tmp1,&Tmp_P.y,2);
    Fp_inv(&tmp1,&tmp1);
    Fp_mul(&tmp2,&Tmp_P.x,&Tmp_P.x);
    Fp_mul_ui(&tmp2,&tmp2,3);
    Fp_mul(lambda,&tmp1,&tmp2);
    Fp_mul(&tmp1,lambda,lambda);
    Fp_mul_ui(&tmp2,&Tmp_P.x,2);
    Fp_sub(&ANS->x,&tmp1,&tmp2);
    Fp_sub(&tmp1,&Tmp_P.x,&ANS->x);
    Fp_mul(&tmp2,lambda,&tmp1);
    Fp_sub(&ANS->y,&tmp2,&Tmp_P.y);
    
    //clear
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    EFp_clear(&Tmp_P);
}

static void bls12_EFp_ECA_return_lambda(EFp *ANS,Fp *lambda,EFp *P1, EFp *P2){
    if(P1->infinity==1){
        EFp_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp_set(ANS,P1);
        return;
    }else if(Fp_cmp(&P1->x,&P2->x)==0){
        if(Fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp_ECD(ANS,P1);
            return;
        }
    }
    
    EFp Tmp_P1,Tmp_P2;
    EFp_init(&Tmp_P1);
    EFp_set(&Tmp_P1,P1);
    EFp_init(&Tmp_P2);
    EFp_set(&Tmp_P2,P2);
    Fp tmp1,tmp2;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    
    Fp_sub(&tmp1,&Tmp_P2.x,&Tmp_P1.x);
    Fp_inv(&tmp1,&tmp1);
    Fp_sub(&tmp2,&Tmp_P2.y,&Tmp_P1.y);
    Fp_mul(lambda,&tmp1,&tmp2);
    Fp_mul(&tmp1,lambda,lambda);
    Fp_sub(&tmp2,&tmp1,&Tmp_P1.x);
    Fp_sub(&ANS->x,&tmp2,&Tmp_P2.x);
    Fp_sub(&tmp1,&Tmp_P1.x,&ANS->x);
    Fp_mul(&tmp2,lambda,&tmp1);
    Fp_sub(&ANS->y,&tmp2,&Tmp_P1.y);
    
    //clear
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    EFp_clear(&Tmp_P1);
    EFp_clear(&Tmp_P2);
}

void bls12_ff_ltt_vtt_for_tate(Fp12 *f,EFp *T,EFp12 *Q){
    EFp Next_T;
    EFp_init(&Next_T);
    Fp12 tmp1,tmp2,tmp3;
    Fp12_init(&tmp1);
    Fp12_init(&tmp2);
    Fp12_init(&tmp3);
    Fp lambda;
    Fp_init(&lambda);
    
    bls12_EFp_ECD_return_lambda(&Next_T,&lambda,T);
    switch(Next_T.infinity){
        case 0:
            Fp12_sqr(&tmp1,f);
            Fp12_set(&tmp3,&Q->x);
            Fp_sub(&tmp3.x0.x0.x0,&Q->x.x0.x0.x0,&T->x);
            Fp12_mul_mpz(&tmp3,&tmp3,lambda.x0);
            Fp12_set(&tmp2,&Q->y);
            Fp_sub(&tmp2.x0.x0.x0,&Q->y.x0.x0.x0,&T->y);
            Fp12_sub(&tmp2,&tmp2,&tmp3);
            
            Fp12_set(&tmp3,&Q->x);
            Fp_sub(&tmp3.x0.x0.x0,&Q->x.x0.x0.x0,&Next_T.x);
            Fp12_inv(&tmp3,&tmp3);
            Fp12_mul(&tmp1,&tmp1,&tmp2);
            Fp12_mul(f,&tmp1,&tmp3);
            break;
        case 1:
            Fp12_sqr(&tmp1,f);
            Fp12_set(&tmp2,&Q->x);
            Fp_sub(&tmp2.x0.x0.x0,&Q->x.x0.x0.x0,&T->x);
            Fp12_mul(f,&tmp1,&tmp2);
            break;
        default:
            break;
    }
    EFp_set(T,&Next_T);
    
    EFp_clear(&Next_T);
    Fp12_clear(&tmp1);
    Fp12_clear(&tmp2);
    Fp12_clear(&tmp3);
    Fp_clear(&lambda);
}

void bls12_f_ltp_vtp_for_tate(Fp12 *f,EFp *T,EFp *P,EFp12 *Q){
    EFp Next_T;
    EFp_init(&Next_T);
    Fp12 tmp1,tmp2,tmp3;
    Fp12_init(&tmp1);
    Fp12_init(&tmp2);
    Fp12_init(&tmp3);
    Fp lambda;
    Fp_init(&lambda);
    
    bls12_EFp_ECA_return_lambda(&Next_T,&lambda,T,P);
    
    switch(Next_T.infinity){
        case 0:
            Fp12_set(&tmp3,&Q->x);
            Fp_sub(&tmp3.x0.x0.x0,&Q->x.x0.x0.x0,&P->x);
            Fp12_mul_mpz(&tmp3,&tmp3,lambda.x0);
            Fp12_set(&tmp1,&Q->y);
            Fp_sub(&tmp1.x0.x0.x0,&Q->y.x0.x0.x0,&P->y);
            Fp12_sub(&tmp1,&tmp1,&tmp3);
            
            Fp12_set(&tmp2,&Q->x);
            Fp_sub(&tmp2.x0.x0.x0,&Q->x.x0.x0.x0,&Next_T.x);
            
            Fp12_inv(&tmp3,&tmp2);
            Fp12_mul(&tmp2,f,&tmp1);
            Fp12_mul(f,&tmp2,&tmp3);
            break;
        case 1:
            Fp12_set(&tmp1,&Q->x);
            Fp_sub(&tmp1.x0.x0.x0,&Q->x.x0.x0.x0,&T->x);
            Fp12_mul(f,f,&tmp1);
            break;
        default:
            break;
    }
    EFp_set(T,&Next_T);
    
    EFp_clear(&Next_T);
    Fp12_clear(&tmp1);
    Fp12_clear(&tmp2);
    Fp12_clear(&tmp3);
    Fp_clear(&lambda);
}

