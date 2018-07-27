//
//  bls12_generate_points.c
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

#include <ELiPS_bn_bls/bls12_generate_points.h>
#include <ELiPS_bn_bls/curve_settings.h>

//rational points
void bls12_generate_G1_point(EFp12 *P){
    EFp Tmp_P;
    EFp_init(&Tmp_P);
    mpz_t scalar;
    mpz_init(scalar);
    
    EFp_rational_point_bls12(&Tmp_P);
    bls12_EFp_to_EFp12(P,&Tmp_P);
    mpz_tdiv_q(scalar,curve_parameters.EFp_total,curve_parameters.order);
    EFp12_SCM(P,P,scalar);
    P->infinity=Tmp_P.infinity;
    
    EFp_clear(&Tmp_P);
    mpz_clear(scalar);
}

void bls12_generate_G2_point(EFp12 *Q){
    EFp12 random_P,P,frobenius_P;
    EFp12_init(&random_P);
    EFp12_init(&P);
    EFp12_init(&frobenius_P);
    mpz_t exp;
    mpz_init(exp);
    
    EFp12_rational_point_bls12(&random_P);
    mpz_pow_ui(exp,curve_parameters.order,2);
    mpz_tdiv_q(exp,curve_parameters.EFpd_total,exp);
    EFp12_SCM(&P,&random_P,exp);
    bls12_Fp12_frobenius_map_p1(&frobenius_P.x,&P.x);
    bls12_Fp12_frobenius_map_p1(&frobenius_P.y,&P.y);
    frobenius_P.infinity=P.infinity;
    EFp12_set_neg(&P,&P);
    EFp12_ECA(Q,&P,&frobenius_P);
    
    mpz_clear(exp);
    EFp12_clear(&random_P);
    EFp12_clear(&P);
    EFp12_clear(&frobenius_P);
}

void bls12_generate_random_point(EFp12 *R){
    EFp12_init(R);
    EFp12_rational_point_bls12(R);
}
