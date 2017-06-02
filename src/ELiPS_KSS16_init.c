//
//  KSS16_init.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_KSS16_Init.h>

void kss16_init (void) {
    set_kss16_params ();
    set_kss16_curve_const ();
    generate_mother_parameter ();
    pre_calculate ();
}
