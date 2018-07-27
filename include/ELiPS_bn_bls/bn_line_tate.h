//
//  bn_line_tate.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/2/18.

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

#ifndef bn_line_tate_h
#define bn_line_tate_h

#include <ELiPS_bn_bls/bn_efp12.h>


extern void EFp_ECD_return_lambda(EFp *ANS,Fp *lambda,EFp *P);
extern void EFp_ECA_return_lambda(EFp *ANS,Fp *lambda,EFp *P1, EFp *P2);
extern void ff_ltt_vtt_for_tate(Fp12 *f,EFp *T,EFp12 *Q);
extern void f_ltp_vtp_for_tate(Fp12 *f,EFp *T,EFp *P,EFp12 *Q);

#endif /* bn_line_tate_h */
