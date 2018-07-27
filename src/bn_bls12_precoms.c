//
//  bn_precoms.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 1/31/18.

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

#include <ELiPS_bn_bls/bn_bls12_precoms.h>

gmp_randstate_t state;

struct Fp Fp_basis;
struct Fp2 Fp2_basis;
struct Fp2 Fp2_basis_inv;
struct Fp6 Fp6_basis;

mpz_t epsilon1,epsilon2;

Fp2 d12_frobenius_constant[d12][6];
Fp2 d12_skew_frobenius_constant[d12][2];

static mpz_t prime_p;

void init_precoms(int curvetype){
    
    mpz_init(prime_p);
    if (curvetype == 3) {
        mpz_set(prime_p,curve_parameters.prime);
    }
    else if (curvetype == 2){
        mpz_set(prime_p,curve_parameters.prime);
    }
    
    Fp_init(&Fp_basis);
    Fp2_init(&Fp2_basis);
    Fp2_init(&Fp2_basis_inv);
    Fp6_init(&Fp6_basis);
    mpz_init(epsilon1);
    mpz_init(epsilon2);
    
    int i,j;
    for(i=0; i<d12; i++){
        for(j=0; j<6; j++){
            Fp2_init(&d12_frobenius_constant[i][j]);
        }
        for(j=0; j<2; j++){
            Fp2_init(&d12_skew_frobenius_constant[i][j]);
        }
    }
    
    
    get_epsilon();
    set_basis();
    set_frobenius_constant();
}

void get_epsilon(){
    Fp inv,buf,result1,result2;
    Fp_init(&inv);
    Fp_init(&buf);
    Fp_init(&result1);
    Fp_init(&result2);
    
    Fp_set_ui(&buf,2);
    Fp_inv(&inv,&buf);
    mpz_sub_ui(buf.x0,prime_p,3);
    
    Fp_sqrt(&buf,&buf);
    Fp_sub_ui(&buf,&buf,1);
    Fp_mul(&result1,&buf,&inv);
    Fp_mul(&result2,&result1,&result1);
    
    mpz_set(epsilon1,result1.x0);
    mpz_set(epsilon2,result2.x0);
    
    Fp_clear(&inv);
    Fp_clear(&buf);
    Fp_clear(&result1);
    Fp_clear(&result2);
}

void set_basis(){
    Fp_set_ui(&Fp_basis,1);
    Fp2_set_ui(&Fp2_basis,1);
    Fp6_set_ui(&Fp6_basis,0);
    Fp_set_ui(&Fp6_basis.x1.x0,1);
    Fp2_inv(&Fp2_basis_inv,&Fp2_basis);
}

void set_frobenius_constant(){
    Fp2 tmp1,tmp2,tmp3;
    Fp2_init(&tmp1);
    Fp2_init(&tmp2);
    Fp2_init(&tmp3);
    
    mpz_t exp,buf,p2,p3,p4,p6,p8,p10;
    mpz_init(exp);
    mpz_init(buf);
    mpz_init(p2);
    mpz_init(p3);
    mpz_init(p4);
    mpz_init(p6);
    mpz_init(p8);
    mpz_init(p10);
    
    mpz_mul(p2,prime_p,prime_p);
    mpz_mul(p3,p2,prime_p);
    mpz_mul(p4,p3,prime_p);
    mpz_mul(p6,p4,p2);
    mpz_mul(p8,p6,p2);
    mpz_mul(p10,p8,p2);
    
    //frobenius_1
    mpz_sub_ui(exp,prime_p,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&tmp1,&Fp2_basis,exp);
    Fp2_mul(&tmp2,&tmp1,&tmp1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp3,&Fp2_basis,exp);
    //set f_p1
    Fp_set_ui(&d12_frobenius_constant[f_p1][0].x0,1);
    Fp2_set(&d12_frobenius_constant[f_p1][1],&tmp1);
    Fp2_set(&d12_frobenius_constant[f_p1][2],&tmp2);
    Fp2_set(&d12_frobenius_constant[f_p1][3],&tmp3);
    Fp2_mul(&d12_frobenius_constant[f_p1][4],&tmp1,&tmp3);
    Fp2_mul(&d12_frobenius_constant[f_p1][5],&tmp2,&tmp3);
    //set skew_f_p1
    Fp2_inv(&tmp1,&tmp1);
    mpz_sub_ui(exp,prime_p,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp2,&Fp2_basis,exp);
    Fp2_inv(&tmp2,&tmp2);
    Fp2_set(&d12_skew_frobenius_constant[f_p1][0],&tmp1);
    Fp2_set(&d12_skew_frobenius_constant[f_p1][1],&tmp2);
    
    //frobenius_2
    mpz_sub_ui(exp,p2,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&tmp1,&Fp2_basis,exp);
    Fp2_mul(&tmp2,&tmp1,&tmp1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp3,&Fp2_basis,exp);
    //set f_p2
    Fp_set_ui(&d12_frobenius_constant[f_p2][0].x0,1);
    Fp2_set(&d12_frobenius_constant[f_p2][1],&tmp1);
    Fp2_set(&d12_frobenius_constant[f_p2][2],&tmp2);
    Fp2_set(&d12_frobenius_constant[f_p2][3],&tmp3);
    Fp2_mul(&d12_frobenius_constant[f_p2][4],&tmp1,&tmp3);
    Fp2_mul(&d12_frobenius_constant[f_p2][5],&tmp2,&tmp3);
    //set skew_f_p2
    Fp2_inv(&tmp1,&tmp1);
    mpz_sub_ui(exp,p2,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp2,&Fp2_basis,exp);
    Fp2_inv(&tmp2,&tmp2);
    Fp2_set(&d12_skew_frobenius_constant[f_p2][0],&tmp1);
    Fp2_set(&d12_skew_frobenius_constant[f_p2][1],&tmp2);
    
    //frobenius_3
    mpz_sub_ui(exp,p3,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&tmp1,&Fp2_basis,exp);
    Fp2_mul(&tmp2,&tmp1,&tmp1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp3,&Fp2_basis,exp);
    //set f_p3
    Fp_set_ui(&d12_frobenius_constant[f_p3][0].x0,1);
    Fp2_set(&d12_frobenius_constant[f_p3][1],&tmp1);
    Fp2_set(&d12_frobenius_constant[f_p3][2],&tmp2);
    Fp2_set(&d12_frobenius_constant[f_p3][3],&tmp3);
    Fp2_mul(&d12_frobenius_constant[f_p3][4],&tmp1,&tmp3);
    Fp2_mul(&d12_frobenius_constant[f_p3][5],&tmp2,&tmp3);
    //set skew_f_p3
    Fp2_inv(&tmp1,&tmp1);
    mpz_sub_ui(exp,p3,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp2,&Fp2_basis,exp);
    Fp2_inv(&tmp2,&tmp2);
    Fp2_set(&d12_skew_frobenius_constant[f_p3][0],&tmp1);
    Fp2_set(&d12_skew_frobenius_constant[f_p3][1],&tmp2);
    
    //d12_frobenius_constant[f_p4]
    mpz_sub_ui(exp,p4,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&tmp1,&Fp2_basis,exp);
    Fp2_mul(&tmp2,&tmp1,&tmp1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp3,&Fp2_basis,exp);
    //set d12_frobenius_constant[f_p4]
    Fp_set_ui(&d12_frobenius_constant[f_p4][0].x0,1);
    Fp2_set(&d12_frobenius_constant[f_p4][1],&tmp1);
    Fp2_set(&d12_frobenius_constant[f_p4][2],&tmp2);
    Fp2_set(&d12_frobenius_constant[f_p4][3],&tmp3);
    Fp2_mul(&d12_frobenius_constant[f_p4][4],&tmp1,&tmp3);
    Fp2_mul(&d12_frobenius_constant[f_p4][5],&tmp2,&tmp3);
    
    //d12_frobenius_constant[f_p8]
    mpz_sub_ui(exp,p8,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&tmp1,&Fp2_basis,exp);
    Fp2_mul(&tmp2,&tmp1,&tmp1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp3,&Fp2_basis,exp);
    //set d12_frobenius_constant[f_p8]
    Fp_set_ui(&d12_frobenius_constant[f_p8][0].x0,1);
    Fp2_set(&d12_frobenius_constant[f_p8][1],&tmp1);
    Fp2_set(&d12_frobenius_constant[f_p8][2],&tmp2);
    Fp2_set(&d12_frobenius_constant[f_p8][3],&tmp3);
    Fp2_mul(&d12_frobenius_constant[f_p8][4],&tmp1,&tmp3);
    Fp2_mul(&d12_frobenius_constant[f_p8][5],&tmp2,&tmp3);
    
    //frobenius_10
    mpz_sub_ui(exp,p10,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&tmp1,&Fp2_basis,exp);
    Fp2_mul(&tmp2,&tmp1,&tmp1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp3,&Fp2_basis,exp);
    //set frobenius_10
    Fp_set_ui(&d12_frobenius_constant[f_p10][0].x0,1);
    Fp2_set(&d12_frobenius_constant[f_p10][1],&tmp1);
    Fp2_set(&d12_frobenius_constant[f_p10][2],&tmp2);
    Fp2_set(&d12_frobenius_constant[f_p10][3],&tmp3);
    Fp2_mul(&d12_frobenius_constant[f_p10][4],&tmp1,&tmp3);
    Fp2_mul(&d12_frobenius_constant[f_p10][5],&tmp2,&tmp3);
    //set skew_f_10
    Fp2_inv(&tmp1,&tmp1);
    mpz_sub_ui(exp,p10,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&tmp2,&Fp2_basis,exp);
    Fp2_inv(&tmp2,&tmp2);
    Fp2_set(&d12_skew_frobenius_constant[f_p10][0],&tmp1);
    Fp2_set(&d12_skew_frobenius_constant[f_p10][1],&tmp2);
    
    Fp2_clear(&tmp1);
    Fp2_clear(&tmp2);
    Fp2_clear(&tmp3);
    mpz_clear(exp);
    mpz_clear(buf);
    mpz_clear(p2);
    mpz_clear(p3);
    mpz_clear(p4);
    mpz_clear(p6);
    mpz_clear(p8);
    mpz_clear(p10);
}

