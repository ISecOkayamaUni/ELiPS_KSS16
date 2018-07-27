//
//  bn_final_exp.c
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

#include <ELiPS_bn_bls/bn_final_exp.h>

void bn_final_exp_plain(Fp12 *ANS,Fp12 *A){
    Fp12 t0,t1;
    Fp12_init(&t0);
    Fp12_init(&t1);
    mpz_t exp,buf;
    mpz_init(exp);
    mpz_init(buf);
    
    Fp12_frobenius_map_p6(&t0,A);
    Fp12_inv(&t1,A);
    Fp12_mul(A,&t0,&t1);
    
    Fp12_frobenius_map_p2(&t0,A);
    Fp12_mul(A,&t0,A);
    
    mpz_pow_ui(exp,curve_parameters.prime,4);
    mpz_pow_ui(buf,curve_parameters.prime,2);
    mpz_sub(exp,exp,buf);
    mpz_add_ui(exp,exp,1);
    mpz_tdiv_q(exp,exp,curve_parameters.order);
    Fp12_pow(ANS,A,exp);
    
    mpz_clear(exp);
    mpz_clear(buf);
    Fp12_clear(&t0);
    Fp12_clear(&t1);
}

void bn_final_exp_optimal(Fp12 *ANS,Fp12 *A){
    Fp12 t0,t1,t2,t3;
    Fp12_init(&t0);
    Fp12_init(&t1);
    Fp12_init(&t2);
    Fp12_init(&t3);
    
    
    //f←f^(p^6)*f^-1
    Fp12_frobenius_map_p6(&t0,A);//f^(p^6)
    Fp12_inv(&t1,A);//f^-1
    Fp12_mul(A,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    Fp12_frobenius_map_p2(&t0,A);//f^(p^2)
    Fp12_mul(A,&t0,A);//f^(p^2)*f
    
    //f←f^((p^4-p^2+1)/r)
    
    bn_fp12_power_motherparam(&t0,A);        //t0←f^(u)
    Fp12_frobenius_map_p6(&t0,&t0);            //t0←f^(-u)
    
    Fp12_sqr(&t1,&t0);                //t1←t0^2
    Fp12_sqr(&t0,&t1);                //t0←t1^2
    Fp12_mul(&t0,&t1,&t0);                //t0←t1*t0
    bn_fp12_power_motherparam(&t2,&t0);    //t2←t0^(u)
    Fp12_frobenius_map_p6(&t2,&t2);            //t2←t0^(-u)
    
    X_binary[0]=0;
    bn_fp12_power_motherparam(&t3,&t2);    //t3←t2^(u+1)
    X_binary[0]=-1;
    Fp12_sqr(&t3,&t3);                //t3←t3^2
    Fp12_frobenius_map_p6(&t0,&t0);            //t0←t0^(-1)
    Fp12_mul(&t3,&t3,&t0);                //t3←t3*t0
    Fp12_frobenius_map_p6(&t2,&t2);            //t2←t2^(-1)
    Fp12_mul(&t2,&t3,&t2);                //t2←t3*t2
    Fp12_mul(&t3,&t3,A);                //t3←t3*f
    
    Fp12_mul(&t1,&t1,&t2);                //t1←t1*t2
    Fp12_frobenius_map_p6(&t0,A);            //t0←f^(-1)
    Fp12_mul(&t0,&t1,&t0);                //t0←t1*t0
    
    Fp12_frobenius_map_p3(&t0,&t0);            //t0←t0^(p^3)
    Fp12_mul(&t0,&t0,&t3);                //t0←t0*t3
    Fp12_frobenius_map_p1(&t1,&t1);            //t1←t1^p
    Fp12_mul(&t0,&t0,&t1);                //t0←t0*t1
    Fp12_frobenius_map_p2(&t2,&t2);            //t2←t2^(p^2)
    Fp12_mul(&t0,&t0,&t2);                //t0←t0*t2
    
    Fp12_set(ANS,&t0);
    
    Fp12_clear(&t0);
    Fp12_clear(&t1);
    Fp12_clear(&t2);
    Fp12_clear(&t3);
}

void bn_fp12_power_motherparam(Fp12 *ANS,Fp12 *A){
    int i;
    Fp12 tmp,A_inv;
    Fp12_init(&tmp);
    Fp12_init(&A_inv);
    Fp12_frobenius_map_p6(&A_inv,A);
    
    Fp12_set(&tmp,A);
    for(i=bn_X_length-1; i>=0; i--){
        switch(X_binary[i]){
            case 0:
                Fp12_sqr(&tmp,&tmp);
                break;
            case 1:
                Fp12_sqr(&tmp,&tmp);
                Fp12_mul(&tmp,&tmp,A);
                break;
            case -1:
                Fp12_sqr(&tmp,&tmp);
                Fp12_mul(&tmp,&tmp,&A_inv);
                break;
            default:
                break;
        }
    }
    Fp12_set(ANS,&tmp);
    
    Fp12_clear(&tmp);
    Fp12_clear(&A_inv);
}
