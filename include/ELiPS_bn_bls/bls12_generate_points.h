//
//  bls12_generate_points.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/6/18.

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

/**
 * @file
 *
 * Header of generating rational point for BLS12 over Fp12 extension field.
 *
 * @ingroup bls12
 */

#ifndef bls12_generate_points_h
#define bls12_generate_points_h

#include <ELiPS_bn_bls/curve_settings.h>
#include <ELiPS_bn_bls/bls12_twist.h>
#include <ELiPS_bn_bls/bls12_frobenius.h>

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/
/**
 * Generates rational point P in G1 group.
 * Usally generated in prime field but maps it to Fp12 for using in pairing.
 * where G1*G2->G3 is the bilinear pairing.
 * 
 * @param[out] P				- the generated point.
 * @param[in] P				    - the input vecotr initialized as 'EFp12 P' point by 'EFp12_init' function.
 */
extern void bls12_generate_G1_point(EFp12 *P);

/**
 * Generates  special rational point Q in G2 group.
 * where G1*G2->G3 is the bilinear pairing.
 * 
 * @param[out] Q				- the generated point.
 * @param[in] Q			        - the input vecotr initialized as 'EFp12 Q' point by 'EFp12_init' function.
 */
extern void bls12_generate_G2_point(EFp12 *Q);

/**
 * Generates random rational point R in the BLS12 curve over Fp12.
 * 
 * @param[out] R				- the generated point.
 * @param[in]  R				- the input vecotr initialized as 'EFp12 P' point by 'EFp12_init'.
 */
extern void bls12_generate_random_point(EFp12 *R);

#endif /* bls12_generate_points_h */
