//
//  bls12_inits.h
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
 * @defgroup bls12 Paring over BARRETO-LYNN-SCOTT curve of embedding degree 12 related files.
 */

/**
 * @file
 *
 * Interaface for pairing using BLS12 curve. It needs to included for initializing BLS12 curve in the applications.
 *
 * @ingroup bls12
 */


#ifndef bls12_inits_h
#define bls12_inits_h

#include <ELiPS_bn_bls/bn_bls12_precoms.h>

/*============================================================================*/
/* Function prototypes                                                        */
/*============================================================================*/
/**
 * This methods needs to be called before using the BLS12 curve for pairing.
 * It initialized the curve as y^2=x^3+4 with parameters suggested in https://eprint.iacr.org/2017/334. 
 */
extern void bls12_inits(void);

#endif /* bls12_inits_h */
