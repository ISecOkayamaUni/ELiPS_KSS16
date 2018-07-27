//
//  bn_settings.c
//  bn_Single_File_Elips
//
//  Created by Y.Nanjo and M.A.A Khandaker on 1/25/18.
//  Copyright © 2018 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_bn_bls/curve_settings.h>
//#include "bn_bls12_precoms.h"


struct curve_params curve_parameters;
char X_binary[bn_X_length+1];
char X_binary_opt[bn_X_length+3];
char X_bit_binary_kss16[KSS16_X_length+1];
//---------------------------------------------------------------------
int bls12_X_length;
mpz_t bls12_X;
int bls12_X_binary[78];

mpz_t C1_INV;

void init_bn_settings(){
    init_bn_parameters();
    
    generate_bn_mother_parameter();
    generate_bn_prime();
    generate_bn_order();
    generate_bn_trace();
    
    weil();
    set_bn_curve_parameter();
}

void init_bls12_settings(){
    init_bls12_parameters();
    
    bls12_generate_X();
    bls12_generate_prime();
    bls12_generate_order();
    bls12_generate_trace();
    
    bls12_weil();
    bls12_set_curve_parameter();
}

void init_kss16_settings(void){
    init_kss16_parameters();
    
    generate_kss16_motherparam();
    generate_kss16_parameters();
}


/*============================================================================*/
/*  init                                                                      */
/*============================================================================*/
void init_bn_parameters(){
    //parameters
    mpz_init(C1_INV);
    mpz_set_str(C1_INV,"307811691015337575251033109375721719032053783174105059253970117417628354711317451266514326962619782791",10);
    
    mpz_init(curve_parameters.prime);
    mpz_init(curve_parameters.X);
    mpz_init(curve_parameters.trace_t);
    mpz_init(curve_parameters.order);
    mpz_init(curve_parameters.EFp_total);
    mpz_init(curve_parameters.EFp2_total);
    mpz_init(curve_parameters.EFp6_total);
    mpz_init(curve_parameters.EFp12_total);
//    mpz_init(curve_parameters.curve_a);
    mpz_init(curve_parameters.curve_b);
    
    int i;
    for(i=0; i<bn_X_length+1; i++){
        X_binary[i]=0;
    }
    for(i=0; i<bn_X_length+3; i++){
        X_binary_opt[i]=0;
    }
}

void init_bls12_parameters(){
    //parameters
    mpz_init(curve_parameters.prime);
    mpz_init(curve_parameters.X);
    mpz_init(curve_parameters.trace_t);
    mpz_init(curve_parameters.order);
    mpz_init(curve_parameters.EFp_total);
    mpz_init(curve_parameters.EFp2_total);
    mpz_init(curve_parameters.EFp6_total);
    mpz_init(curve_parameters.EFp12_total);
    //    mpz_init(curve_parameters.curve_a);
    mpz_init(curve_parameters.curve_b);
    
    int i;
    for(i=0; i<bn_X_length+1; i++){
        X_binary[i]=0;
    }
    for(i=0; i<bn_X_length+3; i++){
        X_binary_opt[i]=0;
    }
    
    bls12_X_length=77;
    for(i=0; i<bls12_X_length+1; i++){
        bls12_X_binary[i]=0;
    }
}

void init_kss16_parameters(void){
    //parameters
    mpz_init(curve_parameters.prime);
    mpz_init(curve_parameters.X);
    mpz_init(curve_parameters.trace_t);
    mpz_init(curve_parameters.order);
    mpz_init(curve_parameters.EFp_total);
    mpz_init(curve_parameters.EFp2_total);
    mpz_init(curve_parameters.EFp6_total);
    mpz_init(curve_parameters.EFp12_total);
    mpz_init(curve_parameters.curve_a);
//    mpz_init(curve_parameters.curve_b);
    
    int i;
    for(i=0; i<KSS16_X_length+1; i++){
        X_bit_binary_kss16[i]=0;
    }
}




void generate_bn_mother_parameter(){
    int i;
    mpz_t buf;
    mpz_init(buf);
    
    //X_binary
    X_binary[114]=1;
    X_binary[101]=1;
    X_binary[14]=-1;
    X_binary[0]=-1;
    
    //X_binary_opt
    X_binary_opt[116]=1;
    X_binary_opt[115]=1;
    X_binary_opt[103]=1;
    X_binary_opt[102]=1;
    X_binary_opt[16]=-1;
    X_binary_opt[15]=-1;
    X_binary_opt[2]=-1;
    
    //X
    mpz_set_ui(curve_parameters.X,0);
    for(i=bn_X_length; i>=0; i--){
        if(X_binary[i]==1){
            mpz_ui_pow_ui(buf,2,i);
            mpz_add(curve_parameters.X,curve_parameters.X,buf);
        }else if(X_binary[i]==-1){
            mpz_ui_pow_ui(buf,2,i);
            mpz_sub(curve_parameters.X,curve_parameters.X,buf);
        }
    }
    
    mpz_clear(buf);
}

int generate_bn_prime(){
    mpz_t buf,result;
    mpz_init(buf);
    mpz_init(result);
    
    //prime
    mpz_pow_ui(buf,curve_parameters.X,4);
    mpz_mul_ui(buf,buf,36);
    mpz_set(result,buf);
    mpz_pow_ui(buf,curve_parameters.X,3);
    mpz_mul_ui(buf,buf,36);
    mpz_add(result,result,buf);
    mpz_pow_ui(buf,curve_parameters.X,2);
    mpz_mul_ui(buf,buf,24);
    mpz_add(result,result,buf);
    mpz_mul_ui(buf,curve_parameters.X,6);
    mpz_add(result,result,buf);
    mpz_add_ui(result,result,1);
    
    //isprime
    if(mpz_probab_prime_p(result,25)==0){
        mpz_clear(buf);
        mpz_clear(result);
        return 0;
    }else{
        mpz_set(curve_parameters.prime,result);
        mpz_clear(buf);
        mpz_clear(result);
        return 1;
    }
}

int  generate_bn_order(){
    mpz_t buf,result;
    mpz_init(buf);
    mpz_init(result);
    
    //prime
    mpz_pow_ui(buf,curve_parameters.X,4);
    mpz_mul_ui(buf,buf,36);
    mpz_set(result,buf);
    mpz_pow_ui(buf,curve_parameters.X,3);
    mpz_mul_ui(buf,buf,36);
    mpz_add(result,result,buf);
    mpz_pow_ui(buf,curve_parameters.X,2);
    mpz_mul_ui(buf,buf,18);
    mpz_add(result,result,buf);
    mpz_mul_ui(buf,curve_parameters.X,6);
    mpz_add(result,result,buf);
    mpz_add_ui(result,result,1);
    
    //isprime
    if(mpz_probab_prime_p(result,25)==0){
        mpz_clear(buf);
        mpz_clear(result);
        return 0;
    }else{
        mpz_set(curve_parameters.order,result);
        mpz_clear(buf);
        mpz_clear(result);
        return 1;
    }
}


void generate_bn_trace(){
    mpz_t buf;
    mpz_init(buf);
    
    mpz_pow_ui(buf,curve_parameters.X,2);
    mpz_mul_ui(buf,buf,6);
    mpz_add_ui(curve_parameters.trace_t,buf,1);
    
    mpz_clear(buf);
}

void set_bn_curve_parameter(){
    mpz_set_ui(curve_parameters.curve_b,4);
}

void weil(){
    mpz_t t2,t6,t12,p2,p6,buf;
    mpz_init(t2);
    mpz_init(t6);
    mpz_init(t12);
    mpz_init(p2);
    mpz_init(p6);
    mpz_init(buf);
    
    //EFp_total
    mpz_add_ui(buf,curve_parameters.prime,1);
    mpz_sub(curve_parameters.EFp_total,buf,curve_parameters.trace_t);
    
    //t2←α^2+β^2
    mpz_pow_ui(t2,curve_parameters.trace_t,2);
    mpz_mul_ui(buf,curve_parameters.prime,2);
    mpz_sub(t2,t2,buf);
    mpz_pow_ui(p2,curve_parameters.prime,2);
    
    //α^6+β^6
    mpz_pow_ui(t6,t2,3);
    mpz_mul(buf,t2,p2);
    mpz_mul_ui(buf,buf,3);
    mpz_sub(t6,t6,buf);
    mpz_pow_ui(p6,p2,3);
    
    //α^12+β^12
    mpz_pow_ui(t12,t6,2);
    mpz_mul_ui(buf,p6,2);
    mpz_sub(t12,t12,buf);
    
    //EFp12_total
    mpz_pow_ui(buf,p6,2);
    mpz_sub(buf,buf,t12);
    mpz_add_ui(curve_parameters.EFpd_total,buf,1);
    
    mpz_clear(t2);
    mpz_clear(t6);
    mpz_clear(t12);
    mpz_clear(p2);
    mpz_clear(p6);
    mpz_clear(buf);
}


static void print_bn_blscurve() {
    gmp_printf("\nelliptic curve\n");
    gmp_printf("E:y^2=x^3-%Zd\n",curve_parameters.curve_b);
}

static void print_kss16curve() {
    gmp_printf("\nKSS16 curve\n");
    gmp_printf("E:y^2=x^3 + %Zd\n",curve_parameters.curve_a);
}

void print_curve_parameters(){
    printf("====================================================================================\n");
    printf("bn12\n\n");
    gmp_printf("parameters\n");
    gmp_printf("X     (%dbit length) : %Zd \n",(int)mpz_sizeinbase(curve_parameters.X,2),curve_parameters.X);
    gmp_printf("prime (%dbit length) : %Zd \n",(int)mpz_sizeinbase(curve_parameters.prime,2),curve_parameters.prime);
    gmp_printf("order (%dbit length) : %Zd \n",(int)mpz_sizeinbase(curve_parameters.order,2),curve_parameters.order);
    gmp_printf("trace (%dbit length) : %Zd \n",(int)mpz_sizeinbase(curve_parameters.trace_t,2),curve_parameters.trace_t);
    
    if (mpz_cmp_ui(curve_parameters. curve_b, 0) > 0) {
         print_bn_blscurve();
    }
    if (mpz_cmp_ui(curve_parameters.curve_a, 0) > 0) {
        print_kss16curve();
    }
//    gmp_printf("\nmodulo polynomial\n");
//    gmp_printf("Fp2  : f(x) = x^2+1\n");
//    gmp_printf("Fp6  : f(x) = x^3-(alpha+1)\n");
//    gmp_printf("Fp12 : f(x) = x^2-beta\n");
    
//    gmp_printf("\nnumber of the total rational points\n");
//    gmp_printf("EFp total   : %Zd\n",curve_parameters.EFp_total);
//    gmp_printf("EFp12 total : %Zd\n",curve_parameters.EFpd_total);
    
//    gmp_printf("\ncubic root of 1\n");
//    gmp_printf("epsilon1 : %Zd\n",epsilon1);
//    gmp_printf("epsilon2 : %Zd\n",epsilon2);
    
}

void bls12_generate_X(){
    int i;
    mpz_t buf,set_2;
    mpz_init(buf);
    mpz_init(set_2);
    mpz_set_ui(set_2,2);
    
    //bls12_X_binary
    bls12_X_binary[77]=-1;
    bls12_X_binary[50]=1;
    bls12_X_binary[33]=1;
    
    //bls12_X
    mpz_init(bls12_X);
    mpz_set_ui(bls12_X,0);
    for(i=bls12_X_length; i>=0; i--){
        if(bls12_X_binary[i]==1){
            mpz_pow_ui(buf,set_2,i);
            mpz_add(bls12_X,bls12_X,buf);
        }else if(bls12_X_binary[i]==-1){
            mpz_pow_ui(buf,set_2,i);
            mpz_sub(bls12_X,bls12_X,buf);
        }
    }
    
    mpz_clear(buf);
    mpz_clear(set_2);
}

int  bls12_generate_prime(){
    mpz_t result,buf1,buf2,modtest;
    mpz_init(result);
    mpz_init(buf1);
    mpz_init(buf2);
    mpz_init(modtest);
    
    mpz_sub_ui(result,bls12_X,1);
    mpz_pow_ui(result,result,2);
    
    mpz_pow_ui(buf1,bls12_X,4);
    mpz_pow_ui(buf2,bls12_X,2);
    mpz_sub(buf1,buf1,buf2);
    mpz_add_ui(buf1,buf1,1);
    
    mpz_mul(result,result,buf1);
    
    //check div3
    mpz_mod_ui(modtest,result,3);
    if(mpz_cmp_ui(modtest,0)!=0){
        mpz_init(result);
        mpz_init(buf1);
        mpz_init(buf2);
        mpz_init(modtest);
        return 0;
    }
    
    mpz_tdiv_q_ui(result,result,3);
    mpz_add(result,result,bls12_X);
    
    //isprime
    if(mpz_probab_prime_p(result,25)==0){
        mpz_init(result);
        mpz_init(buf1);
        mpz_init(buf2);
        mpz_init(modtest);
        return 0;
    }
    
    mpz_set(curve_parameters.prime,result);
    
    mpz_init(result);
    mpz_init(buf1);
    mpz_init(buf2);
    mpz_init(modtest);
    return 1;
}

int  bls12_generate_order(){
    mpz_t buf1,buf2;
    mpz_init(buf1);
    mpz_init(buf2);
    
    mpz_pow_ui(buf1,bls12_X,4);
    mpz_pow_ui(buf2,bls12_X,2);
    mpz_sub(curve_parameters.order,buf1,buf2);
    mpz_add_ui(curve_parameters.order,curve_parameters.order,1);
    
    mpz_clear(buf1);
    mpz_clear(buf2);
    
    return 0;
}

void bls12_generate_trace(){
    mpz_add_ui(curve_parameters.trace_t,bls12_X,1);
}

void bls12_set_curve_parameter(){
    mpz_set_ui(curve_parameters.curve_b,4);
}

void bls12_weil(){
    mpz_t t2,t6,t12,p2,p6,buf;
    mpz_init(t2);
    mpz_init(t6);
    mpz_init(t12);
    mpz_init(p2);
    mpz_init(p6);
    mpz_init(buf);
    
    //EFp_total
    mpz_add_ui(buf,curve_parameters.prime,1);
    mpz_sub(curve_parameters.EFp_total,buf,curve_parameters.trace_t);
    
    //t2←α^2+β^2
    mpz_pow_ui(t2,curve_parameters.trace_t,2);
    mpz_mul_ui(buf,curve_parameters.prime,2);
    mpz_sub(t2,t2,buf);
    mpz_pow_ui(p2,curve_parameters.prime,2);
    
    //α^6+β^6
    mpz_pow_ui(t6,t2,3);
    mpz_mul(buf,t2,p2);
    mpz_mul_ui(buf,buf,3);
    mpz_sub(t6,t6,buf);
    mpz_pow_ui(p6,p2,3);
    
    //α^12+β^12
    mpz_pow_ui(t12,t6,2);
    mpz_mul_ui(buf,p6,2);
    mpz_sub(t12,t12,buf);
    
    //EFp12_total
    mpz_pow_ui(buf,p6,2);
    mpz_sub(buf,buf,t12);
    mpz_add_ui(curve_parameters.EFpd_total,buf,1);
    
    mpz_clear(t2);
    mpz_clear(t6);
    mpz_clear(t12);
    mpz_clear(p2);
    mpz_clear(p6);
    mpz_clear(buf);
}

void bls12_print_parameters(){
    printf("====================================================================================\n");
    printf("bls12\n\n");
    gmp_printf("parameters\n");
    gmp_printf("X     (%dbit length) : %Zd \n",(int)mpz_sizeinbase(bls12_X,2),bls12_X);
    gmp_printf("prime (%dbit length) : %Zd \n",(int)mpz_sizeinbase(curve_parameters.prime,2),curve_parameters.prime);
    gmp_printf("order (%dbit length) : %Zd \n",(int)mpz_sizeinbase(curve_parameters.order,2),curve_parameters.order);
    gmp_printf("trace (%dbit length) : %Zd \n",(int)mpz_sizeinbase(curve_parameters.trace_t,2),curve_parameters.trace_t);
    
    gmp_printf("\nelliptic curve\n");
    gmp_printf("E:y^2=x^3+4\n",curve_parameters.curve_b);
    
    gmp_printf("\nmodulo polynomial\n");
    gmp_printf("Fp2  : f(x) = x^2+1\n");
    gmp_printf("Fp6  : f(x) = x^3-(alpha+1)\n");
    gmp_printf("Fp12 : f(x) = x^2-beta\n");
    
}
//---------------------------------------------------------------------


void generate_kss16_motherparam(void){
    //c1 = 2
    // 2^ -2^32-2^18+2^8+1
    X_bit_binary_kss16[35]=1;
    X_bit_binary_kss16[32]=-1;
    X_bit_binary_kss16[18]=-1;
    X_bit_binary_kss16[8]=1;
    X_bit_binary_kss16[0]=1;
    //2^49+2^26+2^15-2^7-1
//        X_bit_binary_kss16[49]=1;
//        X_bit_binary_kss16[26]=1;
//        X_bit_binary_kss16[15]=1;
//        X_bit_binary_kss16[7]=-1;
//        X_bit_binary_kss16[0]=-1;
    
    mpz_t tmp,set_2;
    mpz_init(tmp);
    mpz_init(set_2);
    mpz_set_ui(set_2,2);
    
    int i;
    for(i=KSS16_X_length;i>=0;i--){
        printf("%d",X_bit_binary_kss16[i]);
        if(X_bit_binary_kss16[i]==1){
            mpz_pow_ui(tmp,set_2,i);
            mpz_add(curve_parameters.X,curve_parameters.X,tmp);
        }else if(X_bit_binary_kss16[i]==-1){
            mpz_pow_ui(tmp,set_2,i);
            mpz_sub(curve_parameters.X,curve_parameters.X,tmp);
        }
    }
    printf("\n");
    mpz_out_str(stdout,10,curve_parameters.X);
    printf("\n");
    return;
}

void generate_kss16_parameters(void){
    
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
    
    mpz_mul(xpow2,curve_parameters.X,curve_parameters.X);
    mpz_mul(xpow4,xpow2,xpow2);
    mpz_mul(xpow5,xpow4,curve_parameters.X);
    mpz_mul(xpow6,xpow5,curve_parameters.X);
    mpz_mul(xpow8,xpow6,xpow2);
    mpz_mul(xpow9,xpow8,curve_parameters.X);
    mpz_mul(xpow10,xpow9,curve_parameters.X);
    
    //t=1/35(2x^5+41x+35)
    mpz_mul_ui(tmp1,curve_parameters.X,41);
    mpz_add_ui(tmp1,tmp1,35);
    mpz_mul_ui(tmp2,xpow5,2);
    mpz_add(t_tmp,tmp1,tmp2);
    
    mpz_div_ui(curve_parameters.trace_t,t_tmp,35);
    
    //r=x^8+48x^4+625
    mpz_mul_ui(tmp1,xpow4,48);
    mpz_add_ui(r_tmp,xpow8,625);
    mpz_add(tmp2,tmp1,r_tmp);
    mpz_tdiv_q_ui(curve_parameters.order,tmp2,61250);
    //     mpz_tdiv_q_ui(order_r,tmp2,49);
    //     mpz_set(order_r,tmp2);
    //    gmp_printf ("order = %Zd\n",order_r);
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
    mpz_mul_ui(tmp1,curve_parameters.X,2398);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_add_ui(p_tmp,p_tmp,3125);
    
    mpz_div_ui(curve_parameters.prime,p_tmp,980);
    
    mpz_add_ui(curve_parameters.EFp_total,curve_parameters.prime,1);
    mpz_sub(curve_parameters.EFp_total,curve_parameters.EFp_total,curve_parameters.trace_t);
    
    if(mpz_probab_prime_p(curve_parameters.prime,25)==0){
        gmp_printf("p:%Zd\n",curve_parameters.prime);
        printf("not  prime number!\n");
        exit(0);
    }

    mpz_set_ui(curve_parameters.curve_a, 1);
    
//    struct EFp P,ANS;
//    int legendre;
//    struct Fp rhs,tmp_ax,x;
//    mpz_init(tmp_a);
//    Fp_init(&rhs);
//    Fp_init(&tmp_ax);
//    EFp_init(&P);
//    EFp_init(&ANS);
//    Fp_init(&x);
//    mpz_init(tmp_a);
//    mpz_set_ui(tmp_a,0);
//
//    for(;;){
//        mpz_add_ui(tmp_a,tmp_a,1);
//        Fp_set_ui(&x,1);
//        legendre=0;
//        while(legendre !=1){
//            mpz_powm_ui(rhs.x0,x.x0,3,PRIME_P);
//            //            gmp_printf("tmp %Zd\n",tmp_a);
//            mpz_mul(tmp_ax.x0,x.x0,tmp_a);
//            Fp_add(&rhs, &rhs, &tmp_ax);
//            if((legendre = mpz_legendre(rhs.x0,PRIME_P))==1){
//                //gmp_printf("a in while = %Zd\n",rhs.x0);
//                Fp_printf(&rhs);
//                Fp_sqrt(&P.y,&rhs);
//                Fp_set(&P.x,&x);
//                EFp_SCM_BIN(&ANS,&P,order_EFp);
//                //                printf("SCM  ==\n");
//                //                EFp_printf(&ANS);
//                if(ANS.infity == TRUE){
//                    mpz_set(a_x,tmp_a);
//                    // mpz_clear(tmp_a);
//                    Fp_clear(&rhs);
//                    Fp_clear(&x);
//                    EFp_clear(&P);
//                    EFp_clear(&ANS);
//                    return;
//                }
//            }
//            Fp_add_ui(&x,&x,1);
//        }
//    }
//    return;
}

