//
//  bn_frobenius.c
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
#include <ELiPS_bn_bls/bn_frobenius.h>

void Fp12_frobenius_map_p1(Fp12 *ANS,Fp12 *A){
    Fp tmp;
    Fp_init(&tmp);
    
    //x0
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    Fp_set(&tmp,&A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1,&tmp);
    Fp2_mul_mpz(&ANS->x0.x1,&ANS->x0.x1,d12_frobenius_constant[f_p1][1].x1.x0);
    Fp_set(&ANS->x0.x2.x0,&A->x0.x2.x0);
    Fp_set_neg(&ANS->x0.x2.x1,&A->x0.x2.x1);
    Fp2_mul_mpz(&ANS->x0.x2,&ANS->x0.x2,d12_frobenius_constant[f_p1][2].x0.x0);
    //x1
    Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    Fp2_mul(&ANS->x1.x0,&ANS->x1.x0,&d12_frobenius_constant[f_p1][3]);
    Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    Fp2_mul(&ANS->x1.x1,&ANS->x1.x1,&d12_frobenius_constant[f_p1][4]);
    Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    Fp2_mul(&ANS->x1.x2,&ANS->x1.x2,&d12_frobenius_constant[f_p1][5]);
    
    Fp_clear(&tmp);
}

void Fp12_frobenius_map_p2(Fp12 *ANS,Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpz(&ANS->x0.x1,&A->x0.x1,d12_frobenius_constant[f_p2][1].x0.x0);
    Fp2_mul_mpz(&ANS->x0.x2,&A->x0.x2,d12_frobenius_constant[f_p2][2].x0.x0);
    //x1
    Fp2_mul_mpz(&ANS->x1.x0,&A->x1.x0,d12_frobenius_constant[f_p2][3].x0.x0);
    Fp2_set_neg(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpz(&ANS->x1.x2,&A->x1.x2,d12_frobenius_constant[f_p2][5].x0.x0);
}

void Fp12_frobenius_map_p3(Fp12 *ANS,Fp12 *A){
    Fp tmp;
    Fp_init(&tmp);
    
    //x0
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    Fp_set(&tmp,&A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1,&tmp);
    Fp_set_neg(&ANS->x0.x2.x0,&A->x0.x2.x0);
    Fp_set(&ANS->x0.x2.x1,&A->x0.x2.x1);
    //x1
    Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    Fp2_mul(&ANS->x1.x0,&ANS->x1.x0,&d12_frobenius_constant[f_p3][3]);
    Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    Fp2_mul(&ANS->x1.x1,&ANS->x1.x1,&d12_frobenius_constant[f_p3][4]);
    Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    Fp2_mul(&ANS->x1.x2,&ANS->x1.x2,&d12_frobenius_constant[f_p3][5]);
    
    Fp_clear(&tmp);
}

void Fp12_frobenius_map_p4(Fp12 *ANS,Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpz(&ANS->x0.x1,&A->x0.x1,d12_frobenius_constant[f_p4][1].x0.x0);
    Fp2_mul_mpz(&ANS->x0.x2,&A->x0.x2,d12_frobenius_constant[f_p4][2].x0.x0);
    //x1
    Fp2_mul_mpz(&ANS->x1.x0,&A->x1.x0,d12_frobenius_constant[f_p4][3].x0.x0);
    Fp2_set(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpz(&ANS->x1.x2,&A->x1.x2,d12_frobenius_constant[f_p4][5].x0.x0);
}

void Fp12_frobenius_map_p6(Fp12 *ANS,Fp12 *A){
    //x0
    Fp6_set(&ANS->x0,&A->x0);
    //x1
    Fp6_set_neg(&ANS->x1,&A->x1);
}

void Fp12_frobenius_map_p8(Fp12 *ANS,Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpz(&ANS->x0.x1,&A->x0.x1,d12_frobenius_constant[f_p8][1].x0.x0);
    Fp2_mul_mpz(&ANS->x0.x2,&A->x0.x2,d12_frobenius_constant[f_p8][2].x0.x0);
    //x1
    Fp2_mul_mpz(&ANS->x1.x0,&A->x1.x0,d12_frobenius_constant[f_p8][3].x0.x0);
    Fp2_set(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpz(&ANS->x1.x2,&A->x1.x2,d12_frobenius_constant[f_p8][5].x0.x0);
}

void Fp12_frobenius_map_p10(Fp12 *ANS,Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpz(&ANS->x0.x1,&A->x0.x1,d12_frobenius_constant[f_p10][1].x0.x0);
    Fp2_mul_mpz(&ANS->x0.x2,&A->x0.x2,d12_frobenius_constant[f_p10][2].x0.x0);
    //x1
    Fp2_mul_mpz(&ANS->x1.x0,&A->x1.x0,d12_frobenius_constant[f_p10][3].x0.x0);
    Fp2_set_neg(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpz(&ANS->x1.x2,&A->x1.x2,d12_frobenius_constant[f_p10][5].x0.x0);
}
