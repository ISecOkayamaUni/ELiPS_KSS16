//
//  bn_clears.c
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

#include <ELiPS_bn_bls/bn_clears.h>

int isCleared = 0;

void clear_parameters(){
    
    //    if (isCleared != 1){
    int i,j;
//    unsigned int base = 10;
    
    mpz_clear(curve_parameters.X);
    mpz_clear(curve_parameters.prime);
    mpz_clear(curve_parameters.order);
    mpz_clear(curve_parameters.trace_t);
    
    mpz_clear(curve_parameters.EFp_total);
    
//    printf(" has exact length %zu in base %d\n", strlen(mpz_get_str(NULL, base, curve_parameters.EFpd_total)), base);
//    printf("Size in base %d\n ",(int)mpz_sizeinbase(curve_parameters.EFpd_total,10));
//    i = (int)strlen(mpz_get_str(NULL, base, curve_parameters.EFpd_total));
//    if (i > 10) {
//    mpz_clear(curve_parameters.EFpd_total);
//    }
//    printf("Size in base after %d\n ",(int)mpz_sizeinbase(curve_parameters.EFpd_total,2));
    mpz_clear(curve_parameters.curve_b);
    mpz_clear(curve_parameters.curve_a);
    
    Fp_clear(&Fp_basis);
    Fp2_clear(&Fp2_basis);
    Fp6_clear(&Fp6_basis);
    mpz_clear(epsilon1);
    mpz_clear(epsilon2);
    
    for(i=0; i<d12; i++){
        for(j=0; j<6; j++){
            Fp2_clear(&d12_frobenius_constant[i][j]);
        }
        for(j=0; j<2; j++){
            Fp2_clear(&d12_skew_frobenius_constant[i][j]);
        }
    }
    //    }
    //    isCleared = 1;
}
