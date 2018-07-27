//
//  bn_skew_frobenius.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/2/18.

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

#include <ELiPS_bn_bls/bn_skew_frobenius.h>

void EFp_skew_frobenius_map_p2(EFp *ANS,EFp *A){
    Fp_mul_mpz(&ANS->x,&A->x,epsilon1);
    Fp_set_neg(&ANS->y,&A->y);
}

void EFp2_skew_frobenius_map_p1(EFp2 *ANS,EFp2 *A){
    //x
    Fp_set(&ANS->x.x0,&A->x.x0);
    Fp_set_neg(&ANS->x.x1,&A->x.x1);
    Fp2_mul(&ANS->x,&ANS->x,&d12_skew_frobenius_constant[f_p1][0]);
    //y
    Fp_set(&ANS->y.x0,&A->y.x0);
    Fp_set_neg(&ANS->y.x1,&A->y.x1);
    Fp2_mul(&ANS->y,&ANS->y,&d12_skew_frobenius_constant[f_p1][1]);
}

void EFp2_skew_frobenius_map_p2(EFp2 *ANS,EFp2 *A){
    //x
    Fp2_mul(&ANS->x,&A->x,&d12_skew_frobenius_constant[f_p2][0]);
    //y
    Fp2_mul(&ANS->y,&A->y,&d12_skew_frobenius_constant[f_p2][1]);
}

void EFp2_skew_frobenius_map_p3(EFp2 *ANS,EFp2 *A){
    //x
    Fp_set(&ANS->x.x0,&A->x.x0);
    Fp_set_neg(&ANS->x.x1,&A->x.x1);
    Fp2_mul(&ANS->x,&ANS->x,&d12_skew_frobenius_constant[f_p3][0]);
    //y
    Fp_set(&ANS->y.x0,&A->y.x0);
    Fp_set_neg(&ANS->y.x1,&A->y.x1);
    Fp2_mul(&ANS->y,&ANS->y,&d12_skew_frobenius_constant[f_p3][1]);
}

void EFp2_skew_frobenius_map_p10(EFp2 *ANS,EFp2 *A){
    //x
    Fp2_mul(&ANS->x,&A->x,&d12_skew_frobenius_constant[f_p10][0]);
    //y
    Fp2_mul(&ANS->y,&A->y,&d12_skew_frobenius_constant[f_p10][1]);
}
