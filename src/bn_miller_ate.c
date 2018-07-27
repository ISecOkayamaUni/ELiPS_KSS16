//
//  bn_miller_ate.c
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

#include <ELiPS_bn_bls/bn_miller_ate.h>

void Miller_algo_for_plain_ate(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    EFp2 T;
    EFp2_init(&T);
    EFp2 mapped_Q;
    EFp2_init(&mapped_Q);
    EFp mapped_P;
    EFp_init(&mapped_P);
    Fp12 f;
    Fp12_init(&f);
    Fp L;
    Fp_init(&L);
    mpz_t loop;
    mpz_init(loop);
    mpz_sub_ui(loop,curve_parameters.trace_t,1);
    int i,length;
    length=(int)mpz_sizeinbase(loop,2);
    char binary[length];
    mpz_get_str(binary,2,loop);
    
    Fp_set(&mapped_P.x,&P->x.x0.x0.x0);
    Fp_set(&mapped_P.y,&P->y.x0.x0.x0);
    mapped_P.infinity=P->infinity;
    EFp12_to_EFp2(&mapped_Q,Q);
    mapped_Q.infinity=Q->infinity;
    
    Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
    
    EFp2_set(&T,&mapped_Q);            //set T
    Fp12_set_ui(&f,0);                //set f
    Fp_set_ui(&f.x0.x0.x0,1);
    
    //miller
    for(i=1; binary[i]!='\0'; i++){
        ff_ltt(&f,&T,&mapped_P,&L);
        if(binary[i]=='1'){
            f_ltq(&f,&T,&mapped_Q,&mapped_P,&L);
        }
    }
    
    Fp12_set(ANS,&f);
    
    Fp12_clear(&f);
    EFp2_clear(&T);
    EFp2_clear(&mapped_Q);
    EFp_clear(&mapped_P);
    Fp_clear(&L);
    mpz_clear(loop);
}

