//
//  bn_utils.c
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

#include <ELiPS_bn_bls/bn_utils.h>

struct timeval t0,t1;

double MILLER_TATE,MILLER_PLAINATE,MILLER_OPTATE,MILLER_XATE;
double FINALEXP_PLAIN,FINALEXP_OPT;
double G1SCM_PLAIN,G1SCM_2SPLIT;
double G2SCM_PLAIN,G2SCM_2SPLIT,G2SCM_4SPLIT;
double G3SCM_PLAIN,G3SCM_2SPLIT,G3SCM_4SPLIT;

double timedifference_msec(struct timeval t0, struct timeval t1){
    return (double)(t1.tv_sec - t0.tv_sec) * 1000.0f + (double)(t1.tv_usec - t0.tv_usec) / 1000.0f;
}

double timedifference_usec(struct timeval t0, struct timeval t1){
    return (double)(t1.tv_sec - t0.tv_sec) * 1000000.0f + (t1.tv_usec - t0.tv_usec);
}

void bn12_print_tate_time(){
    printf("miller tate time:%.2f[ms]\n",MILLER_TATE);
}

void bn12_print_plain_ate_time(){
    printf("miller plain-ate time:%.2f[ms]\n",MILLER_PLAINATE);
}

void bn12_print_opt_ate_time(){
    printf("miller opt-ate time:%.2f[ms]\n",MILLER_OPTATE);
}

void bn12_print_x_ate_time(){
    printf("miller x-ate time:%.2f[ms]\n",MILLER_XATE);
}

void bn12_print_final_exp_plain_time(){
    printf("plain-final exp time:%.2f[ms]\n",FINALEXP_PLAIN);
}

void bn12_print_final_exp_optimal_time(){
    printf("optimal-final exp time:%.2f[ms]\n",FINALEXP_OPT);
}

void bn12_print_G1_plain_time(){
    printf("plain-G1-SCM time:%.2f[ms]\n",G1SCM_PLAIN);
}

void bn12_print_G1_2split_time(){
    printf("2split-G1-SCM time:%.2f[ms]\n",G1SCM_2SPLIT);
}

void bn12_print_G2_plain_time(){
    printf("plain-G2-SCM time:%.2f[ms]\n",G2SCM_PLAIN);
}

void bn12_print_G2_2split_time(){
    printf("2split-G2-SCM time:%.2f[ms]\n",G2SCM_2SPLIT);
}

void bn12_print_G2_4split_time(){
    printf("4split-G2-SCM time:%.2f[ms]\n",G2SCM_4SPLIT);
}

void bn12_print_G3_plain_time(){
    printf("plain-G3-SCM time:%.2f[ms]\n",G3SCM_PLAIN);
}

void bn12_print_G3_2split_time(){
    printf("2split-G3-SCM time:%.2f[ms]\n",G3SCM_2SPLIT);
}

void bn12_print_G3_4split_time(){
    printf("4split-G3-SCM time:%.2f[ms]\n",G3SCM_4SPLIT);
}


