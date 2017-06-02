//
//  KSS16_PS_Miller.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_KSS16_PS_Miller.h>

void ps_miller_kss16(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop){//Q:G2,P:G1
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
    
    EFp4_set(&T,&Q_map);
    Fp4_invert(&L,&P_map.y); // L =yp_bar^-1
    Fp4_set(&z_inv2,&xy_2);
    Fp4_mul(&z_inv2, &z_inv2, &z_inv2);
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    
    //    rational_point_check(&Q_map);
    //    ps_dbl_line_kss16(&ltt,&T,&T,&P_map,&L);
    //    rational_point_check(&T);
    //    ps_add_line_kss16(&ltp,&T,&T,&Q_map,&P_map,&L);
    //    Fp16_printf(&ltp);
    //    rational_point_check(&T);
    
    int i;
    int r_bit;
    r_bit= (int)mpz_sizeinbase(loop,2);
    
    for(i=r_bit-2;i>=0;i--){
        if(mpz_tstbit(loop,i)==1){
            //            printf("\n%d",i);
            Fp16_mul(&l_sum,&l_sum,&l_sum);
            ps_dbl_line_kss16(&ltt,&T,&T,&P_map,&L);
            ps_add_line_kss16(&ltp,&T,&T,&Q_map,&P_map,&L);
            //            rational_point_check(&T);
            ps_mul_line_kss16(&l_sum,&l_sum,&ltt);
            ps_mul_line_kss16(&l_sum,&l_sum,&ltp);
            //            Fp16_mul(&l_sum,&l_sum,&ltt);
            //            Fp16_mul(&l_sum,&l_sum,&ltp);
        }else{
            Fp16_mul(&l_sum,&l_sum,&l_sum);
            ps_dbl_line_kss16(&ltt,&T,&T,&P_map,&L);
            //            rational_point_check(&T);
            ps_mul_line_kss16(&l_sum,&l_sum,&ltt);
            //            Fp16_mul(&l_sum,&l_sum,&ltt);
        }
    }
    // EFp4_printf(&T);
    Fp16_set(ANS,&l_sum);
    
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
}
