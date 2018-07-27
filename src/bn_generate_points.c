//
//  bn_generate_points.c
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

#include <ELiPS_bn_bls/bn_generate_points.h>

EFp12 G1_P;
EFp12 G2_Q;
EFp12 Random_R;

void bn12_generate_G1_point(EFp12 *P){
    EFp12_init(P);
    EFp Tmp;
    EFp_init(&Tmp);
    
    EFp_rational_point_bn(&Tmp);
    EFp12_set_ui(P,0);
    Fp_set(&P->x.x0.x0.x0,&Tmp.x);
    Fp_set(&P->y.x0.x0.x0,&Tmp.y);
    P->infinity=Tmp.infinity;
    
    EFp_clear(&Tmp);
}

//TODO
void bn12_generate_G2_point(EFp12 *Q){
    EFp12_init(Q);
    EFp12 random_P,P,frobenius_P;
    EFp12_init(&random_P);
    EFp12_init(&P);
    EFp12_init(&frobenius_P);
    mpz_t exp;
    mpz_init(exp);
    
    EFp12_rational_point_bn(&random_P);
    mpz_pow_ui(exp,curve_parameters.order,2);
    mpz_tdiv_q(exp,curve_parameters.EFpd_total,exp);
    EFp12_SCM(&P,&random_P,exp);
    Fp12_frobenius_map_p1(&frobenius_P.x,&P.x);
    Fp12_frobenius_map_p1(&frobenius_P.y,&P.y);
    EFp12_set_neg(&P,&P);
    EFp12_ECA(Q,&P,&frobenius_P);
    
    mpz_clear(exp);
    EFp12_clear(&random_P);
    EFp12_clear(&P);
    EFp12_clear(&frobenius_P);
}

void bn12_generate_random_point(EFp12 *R){
    EFp12_init(R);
    EFp12_rational_point_bn(R);
}

