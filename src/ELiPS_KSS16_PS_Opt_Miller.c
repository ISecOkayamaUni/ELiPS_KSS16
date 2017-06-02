//
//  KSS16_PS_Opt_Miller.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_KSS16_PS_Opt_Miller.h>

void ps_opt_miller_kss16(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop){//Q:G2,P:G1
    struct Fp16 l_sum;
    Fp16_init(&l_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    
    struct EFp4 T,P_map,Q_map,EFp4_tmp;
    EFp4_init(&T);
    EFp4_init(&P_map);
    EFp4_init(&Q_map);
    EFp4_init(&EFp4_tmp);
    
    struct Fp4 L,xy,xy_2,y_inv,tmp,y_tmp;
    Fp4_init(&L);
    Fp4_init(&xy);
    Fp4_init(&xy_2);
    Fp4_init(&y_inv);
    Fp4_init(&tmp);
    Fp4_init(&y_tmp);
    Fp4_init(&z_inv2);
    
    Fp4_invert(&y_inv,&P->y);//yp^-1
    Fp4_mul(&xy,&P->x,&y_inv);//xp.yp^-1
    
    Fp4_mul(&xy_2,&xy,&xy);//xy_2 = xp^2.yp^-2
    Fp4_mul(&P_map.x,&xy_2,&P->x);//P.x= xp^3.yp^-2
    Fp4_set(&P_map.y,&P_map.x);
    
    Fp4_mul(&y_tmp,&xy_2,&xy);// xp^2.yp^-2 * xp.yp^-1 = xp^3.yp^-3
    Fp4_mul(&Q_map.y,&y_tmp,&Q->y); //Q_map.y = yQ'.xp^3.yp^-3
    Fp4_mul(&Q_map.x,&xy_2,&Q->x); //Q_map.x = xQ'.xp^2.yp^-2
    
    Fp4_invert(&L,&P_map.y); // L =yp_bar^-1
    Fp4_set(&z_inv2,&xy_2);
    Fp4_mul(&z_inv2, &z_inv2, &z_inv2);
    
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    int i;
    
    struct EFp4 Q_neg;
    EFp4_init(&Q_neg);
    Fp4_neg(&Q_neg.y,&Q_map.y);
    Fp4_set(&Q_neg.x,&Q_map.x);
    
    
    
    if(x_signed_binary[x_bit]==-1){
        EFp4_set(&T,&Q_neg);
    }else{
        EFp4_set(&T,&Q_map);
    }
    
    for(i=x_bit-1;i>=0;i--){
        switch (x_signed_binary[i]){
            case 0:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                ps_dbl_line_kss16(&ltt,&T,&T,&P_map,&L);
                ps_mul_line_kss16(&l_sum,&l_sum,&ltt);
                // Fp16_printf(&l_sum);
                break;
            case 1:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                ps_dbl_line_kss16(&ltt,&T,&T,&P_map,&L);
                ps_add_line_kss16(&ltp,&T,&T,&Q_map,&P_map,&L);
                ps_mul_line_kss16(&l_sum,&l_sum,&ltt);
                ps_mul_line_kss16(&l_sum,&l_sum,&ltp);
                break;
            case -1:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                ps_dbl_line_kss16(&ltt,&T,&T,&P_map,&L);
                ps_add_line_kss16(&ltp,&T,&T,&Q_neg,&P_map,&L);
                ps_mul_line_kss16(&l_sum,&l_sum,&ltt);
                ps_mul_line_kss16(&l_sum,&l_sum,&ltp);
                break;
        }
    }
    
    kss16_skew_frobenius_map(&EFp4_tmp, &Q_map);
    ps_add_line_kss16(&ltp,&T,&T,&EFp4_tmp,&P_map,&L);
    ps_mul_line_kss16(&l_sum,&l_sum,&ltp);
    //  Fp16_mul(&l_sum,&l_sum,&ltp);
    
    struct Fp16 tmp_f;
    Fp16_init(&tmp_f);
    
    Fp16_frobenius_map(&tmp_f, &l_sum);
    Fp16_frobenius_map(&l_sum, &tmp_f);
    Fp16_frobenius_map(&tmp_f, &l_sum);
    
    ps_dbl_line_kss16(&ltt,&T,&Q_map,&P_map,&L);
    ps_mul_line_kss16(&l_sum,&tmp_f,&ltt);
    Fp16_set(ANS,&l_sum);
    
    EFp4_clear(&Q_neg);
    Fp16_clear(&l_sum);
    EFp4_clear(&T);
    EFp4_clear(&P_map);
    EFp4_clear(&Q_map);
    EFp4_clear(&EFp4_tmp);
    Fp4_clear(&L);
    Fp4_clear(&xy);
    Fp4_clear(&xy_2);
    Fp4_clear(&y_inv);
    Fp4_clear(&tmp);
    Fp4_clear(&y_tmp);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
    Fp16_clear(&tmp_f);
}
