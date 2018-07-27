//
//  bn_efp.c
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

#include <ELiPS_bn_bls/bn_efp.h>

void EFp_init(EFp *P){
    Fp_init(&P->x);
    Fp_init(&P->y);
    P->infinity=0;
}

void EFp_clear(EFp *P){
    Fp_clear(&P->x);
    Fp_clear(&P->y);
}

void EFp_printf(EFp *P,char *str){
    gmp_printf("%s",str);
    if(P->infinity==0){
        gmp_printf("(");
        Fp_printf(&P->x,"");
        gmp_printf(",");
        Fp_printf(&P->y,"");
        gmp_printf(")");
    }else{
        gmp_printf("infinity");
    }
}

void EFp_set(EFp *ANS,EFp *P){
    Fp_set(&ANS->x,&P->x);
    Fp_set(&ANS->y,&P->y);
    ANS->infinity=P->infinity;
}

void EFp_set_ui(EFp *ANS,unsigned long int UI){
    Fp_set_ui(&ANS->x,UI);
    Fp_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp_set_mpz(EFp *ANS,mpz_t A){
    Fp_set_mpz(&ANS->x,A);
    Fp_set_mpz(&ANS->y,A);
    ANS->infinity=0;
}

void EFp_set_neg(EFp *ANS,EFp *P){
    Fp_set(&ANS->x,&P->x);
    Fp_set_neg(&ANS->y,&P->y);
    ANS->infinity=P->infinity;
}

void EFp_rational_point_bn(EFp *P){
    Fp tmp1,tmp2,tmp_x;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&tmp_x);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp_set_random(&P->x,state);
        Fp_mul(&tmp1,&P->x,&P->x);
        Fp_mul(&tmp2,&tmp1,&P->x);
        Fp_sub_mpz(&tmp_x,&tmp2,curve_parameters.curve_b);
        if(Fp_legendre(&tmp_x)==1){
            Fp_sqrt(&P->y,&tmp_x);
            break;
        }
    }
    
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&tmp_x);
}

void EFp_rational_point_bls12(EFp *P){
    Fp tmp1,tmp2,tmp_x;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&tmp_x);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp_set_random(&P->x,state);
        Fp_mul(&tmp1,&P->x,&P->x);
        Fp_mul(&tmp2,&tmp1,&P->x);
        Fp_add_mpz(&tmp_x,&tmp2,curve_parameters.curve_b);
        if(Fp_legendre(&tmp_x)==1){
            Fp_sqrt(&P->y,&tmp_x);
            break;
        }
    }
    
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&tmp_x);
}

void EFp_ECD(EFp *ANS,EFp *P){
    if(Fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp Tmp_P;
    EFp_init(&Tmp_P);
    EFp_set(&Tmp_P,P);
    Fp tmp1,tmp2,lambda;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&lambda);
    
    Fp_mul_ui(&tmp1,&Tmp_P.y,2);
    Fp_inv(&tmp1,&tmp1);
    Fp_mul(&tmp2,&Tmp_P.x,&Tmp_P.x);
    Fp_mul_ui(&tmp2,&tmp2,3);
    Fp_mul(&lambda,&tmp1,&tmp2);
    Fp_mul(&tmp1,&lambda,&lambda);
    Fp_mul_ui(&tmp2,&Tmp_P.x,2);
    Fp_sub(&ANS->x,&tmp1,&tmp2);
    Fp_sub(&tmp1,&Tmp_P.x,&ANS->x);
    Fp_mul(&tmp2,&lambda,&tmp1);
    Fp_sub(&ANS->y,&tmp2,&Tmp_P.y);
    
    //clear
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&lambda);
    EFp_clear(&Tmp_P);
}

void EFp_ECA(EFp *ANS,EFp *P1,EFp *P2){
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
    Fp tmp1,tmp2,lambda;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&lambda);
    
    Fp_sub(&tmp1,&Tmp_P2.x,&Tmp_P1.x);
    Fp_inv(&tmp1,&tmp1);
    Fp_sub(&tmp2,&Tmp_P2.y,&Tmp_P1.y);
    Fp_mul(&lambda,&tmp1,&tmp2);
    Fp_mul(&tmp1,&lambda,&lambda);
    Fp_sub(&tmp2,&tmp1,&Tmp_P1.x);
    Fp_sub(&ANS->x,&tmp2,&Tmp_P2.x);
    Fp_sub(&tmp1,&Tmp_P1.x,&ANS->x);
    Fp_mul(&tmp2,&lambda,&tmp1);
    Fp_sub(&ANS->y,&tmp2,&Tmp_P1.y);
    
    //clear
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&lambda);
    EFp_clear(&Tmp_P1);
    EFp_clear(&Tmp_P2);
}

void EFp_SCM(EFp *ANS,EFp *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp_set(ANS,P);
        return;
    }
    
    EFp Tmp_P,Next_P;
    EFp_init(&Tmp_P);
    EFp_set(&Tmp_P,P);
    EFp_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp_set(&Next_P,&Tmp_P);
    for(i=1; i<length; i++){
        EFp_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    
    EFp_set(ANS,&Next_P);
    
    EFp_clear(&Next_P);
    EFp_clear(&Tmp_P);
}

