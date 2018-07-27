//
//  bn_utils.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/1/18.

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

#ifndef bn_utils_h
#define bn_utils_h

#include <stdio.h>
#include <sys/time.h>

extern struct timeval t0,t1;

extern double MILLER_TATE,MILLER_PLAINATE,MILLER_OPTATE,MILLER_XATE;
extern double FINALEXP_PLAIN,FINALEXP_OPT;
extern double G1SCM_PLAIN,G1SCM_2SPLIT;
extern double G2SCM_PLAIN,G2SCM_2SPLIT,G2SCM_4SPLIT;
extern double G3SCM_PLAIN,G3SCM_2SPLIT,G3SCM_4SPLIT;

extern double timedifference_msec(struct timeval t0, struct timeval t1);
extern double timedifference_usec(struct timeval t0, struct timeval t1);

//print pairing
extern void bn12_print_parameters(void);
extern void bn12_print_G1_point(void);
extern void bn12_print_G2_point(void);
extern void bn12_print_tate_time(void);
extern void bn12_print_plain_ate_time(void);
extern void bn12_print_opt_ate_time(void);
extern void bn12_print_x_ate_time(void);
extern void bn12_print_final_exp_plain_time(void);
extern void bn12_print_final_exp_optimal_time(void);
extern void bn12_print_G1_plain_time(void);
extern void bn12_print_G1_2split_time(void);
extern void bn12_print_G2_plain_time(void);
extern void bn12_print_G2_2split_time(void);
extern void bn12_print_G2_4split_time(void);
extern void bn12_print_G3_plain_time(void);
extern void bn12_print_G3_2split_time(void);
extern void bn12_print_G3_4split_time(void);

#endif /* bn_utils_h */
