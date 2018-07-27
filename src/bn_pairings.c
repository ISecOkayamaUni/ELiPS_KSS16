//
//  bn_pairings.c
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

#include <ELiPS_bn_bls/bn_pairings.h>

void bn12_tate(Fp12 *ANS,EFp12 *P, EFp12 *Q){
    //miller
    gettimeofday(&t0,NULL);
    Miller_algo_for_tate(ANS,P,Q);
    gettimeofday(&t1,NULL);
    MILLER_TATE=timedifference_msec(t0,t1);
    
    //final exp
    gettimeofday(&t0,NULL);
    bn_final_exp_optimal(ANS,ANS);
    gettimeofday(&t1,NULL);
    FINALEXP_OPT=timedifference_msec(t0,t1);
}

void bn12_plain_ate(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    //miller
    gettimeofday(&t0,NULL);
    Miller_algo_for_plain_ate(ANS,P,Q);
    gettimeofday(&t1,NULL);
    MILLER_PLAINATE=timedifference_msec(t0,t1);
    
    //final exp
    gettimeofday(&t0,NULL);
    bn_final_exp_optimal(ANS,ANS);
    gettimeofday(&t1,NULL);
    FINALEXP_OPT=timedifference_msec(t0,t1);
}

void bn12_opt_ate(Fp12 *ANS,EFp12 *P, EFp12 *Q){
    //miller
    gettimeofday(&t0,NULL);
    Miller_algo_for_opt_ate(ANS,P,Q);
    gettimeofday(&t1,NULL);
    MILLER_OPTATE=timedifference_msec(t0,t1);
    
    //final exp
    gettimeofday(&t0,NULL);
    bn_final_exp_optimal(ANS,ANS);
    gettimeofday(&t1,NULL);
    FINALEXP_OPT=timedifference_msec(t0,t1);
}

void bn12_x_ate(Fp12 *ANS,EFp12 *P, EFp12 *Q){
    //miller
    gettimeofday(&t0,NULL);
    Miller_algo_for_x_ate(ANS,P,Q);
    gettimeofday(&t1,NULL);
    MILLER_XATE=timedifference_msec(t0,t1);
    
    //final exp
    gettimeofday(&t0,NULL);
    bn_final_exp_optimal(ANS,ANS);
    gettimeofday(&t1,NULL);
    FINALEXP_OPT=timedifference_msec(t0,t1);
}
