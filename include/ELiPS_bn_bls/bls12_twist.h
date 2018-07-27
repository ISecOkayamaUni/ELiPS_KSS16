//
//  bls12_twist.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/5/18.

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
 * Interaface for BLS12's rational point mapping functionalities
 *
 * @ingroup bls12
 */
#ifndef bls12_twist_h
#define bls12_twist_h

#include <ELiPS_bn_bls/bn_efp12.h>

/**
 * Calculate sextic twisted map G2 rational point
 *
 * @param[out] ANS			    -  output in EFp2
 * @param[in] P			        -  input P in EFp12.
 */
extern void bls12_EFp12_to_EFp2(EFp2 *ANS,EFp12 *P);

/**
 * Calculate sextic twisted map G2 rational point
 *
 * @param[out] ANS			    -  output in EFp12
 * @param[in] P			        -  input P in EFp2.
 */
extern void bls12_EFp2_to_EFp12(EFp12 *ANS,EFp2 *P);

/**
 * Calculate primary twist of  G1 rational point
 *
 * @param[out] ANS			    -  output in EFp
 * @param[in] P			        -  input P in EFp12.
 */
extern void bls12_EFp12_to_EFp(EFp *ANS,EFp12 *P);

/**
 * Calculate primary twist of  G1 rational point
 *
 * @param[out] ANS			    -  output in EFp12
 * @param[in] P			        -  input P in EFp.
 */
extern void bls12_EFp_to_EFp12(EFp12 *ANS,EFp *P);

#endif /* bls12_twist_h */
