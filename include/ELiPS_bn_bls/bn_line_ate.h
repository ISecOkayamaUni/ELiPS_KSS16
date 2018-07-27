//
//  bn_line_ate.h
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

#ifndef bn_line_ate_h
#define bn_line_ate_h

#include <ELiPS_bn_bls/bn_p8sparse.h>

extern void ff_ltt(Fp12 *f,EFp2 *T,EFp *P,Fp *L);
extern void f_ltq(Fp12 *f,EFp2 *T,EFp2 *Q,EFp *P,Fp *L);

#endif /* bn_line_ate_h */
