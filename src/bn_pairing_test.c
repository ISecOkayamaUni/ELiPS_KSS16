//
//  bn_pairing_test.c
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

#include <ELiPS_bn_bls/bn_pairing_test.h>

void bn12_test_tate_pairing(){
    printf("====================================================================================\n");
    printf("tate pairing\n\n");
    EFp12 P,Q,s1_P,s2_P,s1_Q,s2_Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&s1_P);
    EFp12_init(&s2_P);
    EFp12_init(&s1_Q);
    EFp12_init(&s2_Q);
    Fp12 Z,Test1,Test2,Test3;
    Fp12_init(&Z);
    Fp12_init(&Test1);
    Fp12_init(&Test2);
    Fp12_init(&Test3);
    mpz_t s1,s2,s12;
    mpz_init(s1);
    mpz_init(s2);
    mpz_init(s12);
    
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(s1,state,curve_parameters.order);
    mpz_urandomm(s2,state,curve_parameters.order);
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,curve_parameters.order);
    
    printf("input\n");
    bn12_generate_G1_point(&P);
    EFp12_printf(&P,"P : G1 rational point\n");
    printf("\n\n");
    EFp12_rational_point_bn(&Q);
    EFp12_printf(&Q,"Q : random rational point\n");
    printf("\n\n");
    
    EFp12_SCM(&s1_P,&P,s1);
    EFp12_SCM(&s2_P,&P,s2);
    EFp12_SCM(&s1_Q,&Q,s1);
    EFp12_SCM(&s2_Q,&Q,s2);
    
    printf("bilinearity test\n");
    //test1
    printf("tate(P,Q)^s12\n");
    bn12_tate(&Z,&P,&Q);
    Fp12_pow(&Test1,&Z,s12);
    bn12_print_tate_time();
    bn12_print_final_exp_optimal_time();
    Fp12_printf(&Test1,"");
    printf("\n\n");
    //test2
    printf("tate([s1]P,[s2]Q)\n");
    bn12_tate(&Test2,&s1_P,&s2_Q);
    bn12_print_tate_time();
    bn12_print_final_exp_optimal_time();
    Fp12_printf(&Test2,"");
    printf("\n\n");
    //test3
    printf("tate([s2]P,[s1]Q)\n");
    bn12_tate(&Test3,&s2_P,&s1_Q);
    bn12_print_tate_time();
    bn12_print_final_exp_optimal_time();
    Fp12_printf(&Test3,"");
    printf("\n\n");
    
    if(Fp12_cmp_zero(&Test1)!=0 && Fp12_cmp_one(&Test1)!=0 && Fp12_cmp(&Test1,&Test2)==0 && Fp12_cmp(&Test2,&Test3)==0 && Fp12_cmp(&Test3,&Test1)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(s1);
    mpz_clear(s2);
    mpz_clear(s12);
    EFp12_clear(&P);
    EFp12_clear(&Q);
    EFp12_clear(&s1_P);
    EFp12_clear(&s2_P);
    EFp12_clear(&s1_Q);
    EFp12_clear(&s2_Q);
    Fp12_clear(&Z);
    Fp12_clear(&Test1);
    Fp12_clear(&Test2);
    Fp12_clear(&Test3);
}

void bn12_test_plain_ate_pairing(){
    printf("====================================================================================\n");
    printf("plain-ate pairing\n\n");
    EFp12 P,Q,s1_P,s2_P,s1_Q,s2_Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&s1_P);
    EFp12_init(&s2_P);
    EFp12_init(&s1_Q);
    EFp12_init(&s2_Q);
    Fp12 Z,Test1,Test2,Test3;
    Fp12_init(&Z);
    Fp12_init(&Test1);
    Fp12_init(&Test2);
    Fp12_init(&Test3);
    mpz_t s1,s2,s12;
    mpz_init(s1);
    mpz_init(s2);
    mpz_init(s12);
    
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(s1,state,curve_parameters.order);
    mpz_urandomm(s2,state,curve_parameters.order);
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,curve_parameters.order);
    
    printf("input\n");
    bn12_generate_G1_point(&P);
    EFp12_printf(&P,"P : G1 rational point\n");
    printf("\n");
    bn12_generate_G2_point(&Q);
    EFp12_printf(&Q,"Q : G2 rational point\n");
    printf("\n\n");
    
    EFp12_SCM(&s1_P,&P,s1);
    EFp12_SCM(&s2_P,&P,s2);
    EFp12_SCM(&s1_Q,&Q,s1);
    EFp12_SCM(&s2_Q,&Q,s2);
    
    printf("bilinearity test\n");
    //test1
    printf("plain-ate(Q,P)^s12\n");
    bn12_plain_ate(&Z,&P,&Q);
    Fp12_pow(&Test1,&Z,s12);
    bn12_print_plain_ate_time();
    bn12_print_final_exp_optimal_time();
    Fp12_printf(&Test1,""); printf("\n\n");
    //test2
    printf("plain-ate([s2]Q,[s1]P)\n");
    bn12_plain_ate(&Test2,&s2_P,&s1_Q);
    bn12_print_plain_ate_time();
    bn12_print_final_exp_optimal_time();
    Fp12_printf(&Test2,""); printf("\n\n");
    //test3
    printf("plain-ate([s1]Q,[s2]P)\n");
    bn12_plain_ate(&Test3,&s1_P,&s2_Q);
    bn12_print_plain_ate_time();
    bn12_print_final_exp_optimal_time();
    Fp12_printf(&Test3,""); printf("\n\n");
    
    if(Fp12_cmp_zero(&Test1)!=0 && Fp12_cmp_one(&Test1)!=0 && Fp12_cmp(&Test1,&Test2)==0 && Fp12_cmp(&Test2,&Test3)==0 && Fp12_cmp(&Test3,&Test1)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    
    mpz_clear(s1);
    mpz_clear(s2);
    mpz_clear(s12);
    EFp12_clear(&P);
    EFp12_clear(&Q);
    EFp12_clear(&s1_P);
    EFp12_clear(&s2_P);
    EFp12_clear(&s1_Q);
    EFp12_clear(&s2_Q);
    Fp12_clear(&Z);
    Fp12_clear(&Test1);
    Fp12_clear(&Test2);
    Fp12_clear(&Test3);
}

void bn12_test_opt_ate_pairing(){
    printf("====================================================================================\n");
    printf("opt-ate pairing\n\n");
    EFp12 P,Q,s1_P,s2_P,s1_Q,s2_Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&s1_P);
    EFp12_init(&s2_P);
    EFp12_init(&s1_Q);
    EFp12_init(&s2_Q);
    Fp12 Z,Test1,Test2,Test3;
    Fp12_init(&Z);
    Fp12_init(&Test1);
    Fp12_init(&Test2);
    Fp12_init(&Test3);
    mpz_t s1,s2,s12;
    mpz_init(s1);
    mpz_init(s2);
    mpz_init(s12);
    
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(s1,state,curve_parameters.order);
    mpz_urandomm(s2,state,curve_parameters.order);
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,curve_parameters.order);
    
    printf("input\n");
    bn12_generate_G1_point(&P);
    EFp12_printf(&P,"P : G1 rational point\n");
    printf("\n");
    bn12_generate_G2_point(&Q);
    EFp12_printf(&Q,"Q : G2 rational point\n");
    printf("\n\n");
    
    EFp12_SCM(&s1_P,&P,s1);
    EFp12_SCM(&s2_P,&P,s2);
    EFp12_SCM(&s1_Q,&Q,s1);
    EFp12_SCM(&s2_Q,&Q,s2);
    
    printf("bilinearity test\n");
    //test1
    printf("opt-ate(Q,P)^s12\n");
    bn12_opt_ate(&Z,&P,&Q);
    Fp12_pow(&Test1,&Z,s12);
    bn12_print_opt_ate_time();
    bn12_print_final_exp_optimal_time();
    Fp12_printf(&Test1,""); printf("\n\n");
    //test2
    printf("opt-ate([s2]Q,[s1]P)\n");
    bn12_opt_ate(&Test2,&s2_P,&s1_Q);
    bn12_print_opt_ate_time();
    bn12_print_final_exp_optimal_time();
    Fp12_printf(&Test2,""); printf("\n\n");
    //test3
    printf("opt-ate([s1]Q,[s2]P)\n");
    bn12_opt_ate(&Test3,&s1_P,&s2_Q);
    bn12_print_opt_ate_time();
    bn12_print_final_exp_optimal_time();
    Fp12_printf(&Test3,""); printf("\n\n");
    
    if(Fp12_cmp_zero(&Test1)!=0 && Fp12_cmp_one(&Test1)!=0 && Fp12_cmp(&Test1,&Test2)==0 && Fp12_cmp(&Test2,&Test3)==0 && Fp12_cmp(&Test3,&Test1)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    
    mpz_clear(s1);
    mpz_clear(s2);
    mpz_clear(s12);
    EFp12_clear(&P);
    EFp12_clear(&Q);
    EFp12_clear(&s1_P);
    EFp12_clear(&s2_P);
    EFp12_clear(&s1_Q);
    EFp12_clear(&s2_Q);
    Fp12_clear(&Z);
    Fp12_clear(&Test1);
    Fp12_clear(&Test2);
    Fp12_clear(&Test3);
}

void bn12_test_x_ate_pairing(){
    printf("====================================================================================\n");
    printf("x-ate pairing\n\n");
    EFp12 P,Q,s1_P,s2_P,s1_Q,s2_Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&s1_P);
    EFp12_init(&s2_P);
    EFp12_init(&s1_Q);
    EFp12_init(&s2_Q);
    Fp12 Z,Test1,Test2,Test3;
    Fp12_init(&Z);
    Fp12_init(&Test1);
    Fp12_init(&Test2);
    Fp12_init(&Test3);
    mpz_t s1,s2,s12;
    mpz_init(s1);
    mpz_init(s2);
    mpz_init(s12);
    
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(s1,state,curve_parameters.order);
    mpz_urandomm(s2,state,curve_parameters.order);
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,curve_parameters.order);
    
    printf("input\n");
    bn12_generate_G1_point(&P);
    EFp12_printf(&P,"P : G1 rational point\n");
    printf("\n");
    bn12_generate_G2_point(&Q);
    EFp12_printf(&Q,"Q : G2 rational point\n");
    printf("\n\n");
    
    EFp12_SCM(&s1_P,&P,s1);
    EFp12_SCM(&s2_P,&P,s2);
    EFp12_SCM(&s1_Q,&Q,s1);
    EFp12_SCM(&s2_Q,&Q,s2);
    
    printf("bilinearity test\n");
    //test1
    printf("x-ate(Q,P)^s12\n");
    bn12_x_ate(&Z,&P,&Q);
    Fp12_pow(&Test1,&Z,s12);
    bn12_print_x_ate_time();
    bn12_print_final_exp_optimal_time();
    Fp12_printf(&Test1,""); printf("\n\n");
    //test2
    printf("x-ate([s2]Q,[s1]P)\n");
    bn12_x_ate(&Test2,&s2_P,&s1_Q);
    bn12_print_x_ate_time();
    bn12_print_final_exp_optimal_time();
    Fp12_printf(&Test2,""); printf("\n\n");
    //test3
    printf("x-ate([s1]Q,[s2]P)\n");
    bn12_x_ate(&Test3,&s1_P,&s2_Q);
    bn12_print_x_ate_time();
    bn12_print_final_exp_optimal_time();
    Fp12_printf(&Test3,""); printf("\n\n");
    
    if(Fp12_cmp_zero(&Test1)!=0 && Fp12_cmp_one(&Test1)!=0 && Fp12_cmp(&Test1,&Test2)==0 && Fp12_cmp(&Test2,&Test3)==0 && Fp12_cmp(&Test3,&Test1)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    
    mpz_clear(s1);
    mpz_clear(s2);
    mpz_clear(s12);
    EFp12_clear(&P);
    EFp12_clear(&Q);
    EFp12_clear(&s1_P);
    EFp12_clear(&s2_P);
    EFp12_clear(&s1_Q);
    EFp12_clear(&s2_Q);
    Fp12_clear(&Z);
    Fp12_clear(&Test1);
    Fp12_clear(&Test2);
    Fp12_clear(&Test3);
}
