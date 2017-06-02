//
//  Settings.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/26/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_KSS16_Settings.h>

int c1 = 2; //Quardratic not residue over Prime field Fp
int is_params_set = 0;
int TRUE = TRUE_1;
int FALSE = FALSE_0;
char x_signed_binary[x_bit+1];

struct kss16_params params;
struct KSS16_constants kss_curve_const;

void set_kss16_params(){
    
    if (!is_params_set) {
        mpz_init(params.prime);
        mpz_init(params.order_r);
        mpz_init(params.trace_t);
        mpz_init(params.order_EFp);
        mpz_init(params.X);
        
        mpz_set_str(params.prime,PRIME,10);
        mpz_set_str(params.order_r,ORDER_R,10);
        mpz_set_str(params.trace_t,TRACE_T,10);
        mpz_set_str(params.order_EFp,ORDER_EFp,10);
        
        is_params_set = 1;
    }
}

void set_kss16_curve_const()
{
    mpz_init(kss_curve_const.a);
    mpz_init(kss_curve_const.tmp_a);
    
    mpz_set_ui(kss_curve_const.a,1);
    mpz_set_ui(kss_curve_const.tmp_a,1);
}


void generate_mother_parameter(){
    //c1 = 2
    // 2^35-2^32-2^18+2^8+1
    mpz_init(params.X);
    
    x_signed_binary[35]=1;
    x_signed_binary[32]=-1;
    x_signed_binary[18]=-1;
    x_signed_binary[8]=1;
    x_signed_binary[0]=1;
    
    mpz_t tmp,set_2;
    mpz_init(tmp);
    mpz_init(set_2);
    mpz_set_ui(set_2,2);
    
    int i;
    for(i=x_bit;i>=0;i--){
        printf("%d",x_signed_binary[i]);
        if(x_signed_binary[i]==1){
            mpz_pow_ui(tmp,set_2,i);
            mpz_add(params.X,params.X,tmp);
        }else if(x_signed_binary[i]==-1){
            mpz_pow_ui(tmp,set_2,i);
            mpz_sub(params.X,params.X,tmp);
        }
    }
    printf("\n");
    mpz_out_str(stdout,2,params.X);
    printf("\n");
    return;
}

float timedifference_msec(struct timeval t0, struct timeval t1){
    return (t1.tv_sec - t0.tv_sec) * 1000.0f + (t1.tv_usec - t0.tv_usec) / 1000.0f;
}
