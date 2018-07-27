//
//  bn_skew_frobenius.h
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

#ifndef bn_skew_frobenius_h
#define bn_skew_frobenius_h

#include <ELiPS_bn_bls/bn_bls12_precoms.h>
#include <ELiPS_bn_bls/bn_efp12.h>

extern void EFp_skew_frobenius_map_p2(EFp *ANS,EFp *A);
extern void EFp2_skew_frobenius_map_p1(EFp2 *ANS,EFp2 *A);
extern void EFp2_skew_frobenius_map_p2(EFp2 *ANS,EFp2 *A);
extern void EFp2_skew_frobenius_map_p3(EFp2 *ANS,EFp2 *A);
extern void EFp2_skew_frobenius_map_p10(EFp2 *ANS,EFp2 *A);

#endif /* bn_skew_frobenius_h */
