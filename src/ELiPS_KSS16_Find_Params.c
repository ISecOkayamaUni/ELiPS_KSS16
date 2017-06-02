//
//  KSS16_Find_Params.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/27/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_KSS16_Params.h>
#include <stdlib.h>

void generate_kss16_params(void){
    
    mpz_t tmp1,tmp2,two;
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init(two);
    
    
    //set p,r
    mpz_t p_tmp,r_tmp,t_tmp;
    mpz_t xpow2,xpow4,xpow5,xpow6,xpow8,xpow9,xpow10;
    //    mpz_t tmp1,tmp2;
    
    mpz_init(p_tmp);
    mpz_init(r_tmp);
    mpz_init(t_tmp);
    mpz_init(xpow2);
    mpz_init(xpow4);
    mpz_init(xpow5);
    mpz_init(xpow6);
    mpz_init(xpow8);
    mpz_init(xpow9);
    mpz_init(xpow10);
    mpz_init(tmp1);
    mpz_init(tmp2);
    
    mpz_mul(xpow2,params.X,params.X);
    mpz_mul(xpow4,xpow2,xpow2);
    mpz_mul(xpow5,xpow4,params.X);
    mpz_mul(xpow6,xpow5,params.X);
    mpz_mul(xpow8,xpow6,xpow2);
    mpz_mul(xpow9,xpow8,params.X);
    mpz_mul(xpow10,xpow9,params.X);
    
    //t=1/35(2x^5+41x+35)
    mpz_mul_ui(tmp1,params.X,41);
    mpz_add_ui(tmp1,tmp1,35);
    mpz_mul_ui(tmp2,xpow5,2);
    mpz_add(t_tmp,tmp1,tmp2);
    
    mpz_div_ui(params.trace_t,t_tmp,35);
    
    gmp_printf ("trace = %Zd\n",params.trace_t);
    //r=x^8+48x^4+625
    mpz_mul_ui(tmp1,xpow4,48);
    mpz_add_ui(r_tmp,xpow8,625);
    mpz_add(tmp2,tmp1,r_tmp);
    mpz_tdiv_q_ui(params.order_r,tmp2,61250);
    //     mpz_tdiv_q_ui(order_r,tmp2,49);
    //     mpz_set(order_r,tmp2);
    gmp_printf ("order = %Zd\n",params.order_r);
    // mpz_set(r,r_tmp);
    
    //p=1/980(x^10+2x^9+5x^8+48x^6+152x^5+240x^4+625x^2+2398x+3125)
    mpz_mul_ui(tmp1,xpow9,2);
    mpz_add(p_tmp,tmp1,xpow10);
    mpz_mul_ui(tmp1,xpow8,5);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_mul_ui(tmp1,xpow6,48);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_mul_ui(tmp1,xpow5,152);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_mul_ui(tmp1,xpow4,240);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_mul_ui(tmp1,xpow2,625);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_mul_ui(tmp1,params.X,2398);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_add_ui(p_tmp,p_tmp,3125);
    
    mpz_div_ui(params.prime,p_tmp,980);
    // mpz_set(p,p_tmp);
    
    //        mpz_t mod, p1;
    //        mpz_init(mod);
    //        mpz_init(p1);
    //
    //        mpz_pow_ui(p1,prime,8);
    //        mpz_add_ui(p1,p1,1);
    //
    //        int l = (int)mpz_divisible_p(p1,order_r);
    //        int k = (int)mpz_divisible_ui_p(p1,16);
    //        gmp_printf("\n\n l , k  %d == %d\n",l,k);
    //
    //        mpz_set_ui(mod,8);
    //
    //        mpz_mod(mod,p1,order_r);
    //        gmp_printf("\n == %Zd\n",mod);
    //
    //    gmp_printf("p:%Zd\n",prime);
    
    mpz_add_ui(params.order_EFp,params.prime,1);
    mpz_sub(params.order_EFp,params.order_EFp,params.trace_t);
    
    if(mpz_probab_prime_p(params.prime,25)==0){
        gmp_printf("p:%Zd\n",params.prime);
        printf("not  prime number!\n");
        exit(0);
    }
    
//    if (mpz_divisible_p(params.order_EFp,params.order_r) == 0) {
//        printf("Not Divisible \n");
//    }
//    else{
//        printf("Divisible \n");
//    }
    
    struct EFp P,ANS;
    int legendre;
    struct Fp rhs,tmp_ax,x;
    mpz_init(kss_curve_const.tmp_a);
    Fp_init(&rhs);
    Fp_init(&tmp_ax);
    EFp_init(&P);
    EFp_init(&ANS);
    Fp_init(&x);
    mpz_init(kss_curve_const.tmp_a);
    mpz_set_ui(kss_curve_const.tmp_a,0);
    
    for(;;){
        mpz_add_ui(kss_curve_const.tmp_a,kss_curve_const.tmp_a,1);
        Fp_set_ui(&x,1);
        legendre=0;
        while(legendre !=1){
            mpz_powm_ui(rhs.x0,x.x0,3,params.prime);
//            gmp_printf("tmp %Zd\n",kss_curve_const.tmp_a);
            mpz_mul(tmp_ax.x0,x.x0,kss_curve_const.tmp_a);
            Fp_add(&rhs, &rhs, &tmp_ax);
            if((legendre = mpz_legendre(rhs.x0,params.prime))==1){
                //gmp_printf("a in while = %Zd\n",rhs.x0);
                Fp_printf(&rhs);
                Fp_sqrt(&P.y,&rhs);
                Fp_set(&P.x,&x);
                EFp_scm_bin(&ANS,&P,params.order_EFp);
//                printf("SCM\n");
//                EFp_printf(&ANS);
                if(ANS.PoI == TRUE){
                    mpz_set(kss_curve_const.a,kss_curve_const.tmp_a);
                    // mpz_clear(tmp_a);
                    Fp_clear(&rhs);
                    Fp_clear(&x);
                    EFp_clear(&P);
                    EFp_clear(&ANS);
                    return;
                }
            }
            Fp_add_ui(&x,&x,1);
        }
    }
    return;
}
