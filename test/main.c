//
//  main.c
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 1/25/18.

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

#include <ELiPS_bn_bls/bn_inits.h>
#include <ELiPS_bn_bls/bn_clears.h>
#include <ELiPS_bn_bls/bn_pairing_test.h>

#include <ELiPS_bn_bls/bls12_inits.h>
#include <ELiPS_bn_bls/bls12_pairings.h>
#include <ELiPS_bn_bls/bls12_test_pairings.h>


int main(int argc, const char * argv[]) {
    
    bls12_inits();
    clear_parameters();
    // bls12_print_parameters();
    // bls12_test_opt_ate_pairing();
    // clear_parameters();
    
   // init_bn();
    bls12_inits();
    mpz_t scalar1, scalar2;
    mpz_init(scalar1);
    mpz_init(scalar2);
    mpz_set_str(scalar1,"59049",10);
    mpz_set_str(scalar2,"46656",10);

    struct EFp12 P,Q, ANS;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&ANS);

    bls12_generate_G2_point(&P);
    printf("\n");
    printf("\n");
    EFp12_printf(&P,"G2 point=");
    EFp12_SCM(&ANS,&P,scalar1);
     printf("\n");
    printf("\n");
    EFp12_printf(&ANS,"SCL 1=");
    
    
//    EFp12_ECD(&Q,&P);
    EFp12_SCM(&Q,&P,scalar2);
    printf("\n");
    printf("\n");
    EFp12_printf(&P,"G2 point=");
    printf("\n");
    printf("\n");
    EFp12_printf(&Q,"SCL 2=");
     printf("\n");
    
  // print_curve_parameters();
    // bn12_test_x_ate_pairing();
    // bn12_test_plain_ate_pairing();
    // bn12_test_opt_ate_pairing();
    clear_parameters();
    
    return 0;
}
