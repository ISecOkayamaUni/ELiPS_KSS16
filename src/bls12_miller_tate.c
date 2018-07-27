//
//  bls12_miller_tate.c
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

#include <ELiPS_bn_bls/bls12_miller_tate.h>


void bls12_Miller_algo_for_tate(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    EFp Tmp_P,T;
    EFp_init(&Tmp_P);
    EFp_init(&T);
    Fp12 f;
    Fp12_init(&f);
    int i,length;
    length=(int)mpz_sizeinbase(curve_parameters.order,2);
    char binary[length];
    mpz_get_str(binary,2,curve_parameters.order);
    
    //set
    bls12_EFp12_to_EFp(&Tmp_P,P);
    EFp_set(&T,&Tmp_P);
    Fp12_set_ui(&f,0);
    Fp_set_ui(&f.x0.x0.x0,1);
    
    //miller
    for(i=1; i<length; i++){
        bls12_ff_ltt_vtt_for_tate(&f,&T,Q);
        if(binary[i]=='1'){
            bls12_f_ltp_vtp_for_tate(&f,&T,&Tmp_P,Q);
        }
    }
    Fp12_set(ANS,&f);
    
    Fp12_clear(&f);
    EFp_clear(&T);
    EFp_clear(&Tmp_P);
}

