//
//  Fp16.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/26/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_Fp16.h>

// #pragma mark Fp8 methods
void Fp16_init(struct Fp16 *A){
    Fp8_init(&A->x0);
    Fp8_init(&A->x1);
}
void Fp16_set(struct Fp16 *ANS,struct Fp16 *A){
    Fp8_set(&ANS->x0,&A->x0);
    Fp8_set(&ANS->x1,&A->x1);
}
void Fp16_set_ui(struct Fp16 *A,signed long int B){
    Fp8_set_ui(&A->x0,B);
    Fp8_set_ui(&A->x1,B);
}
void Fp16_random(struct Fp16 *A){
    Fp8_random(&A->x0);
    Fp8_random(&A->x1);
}
void Fp16_clear(struct Fp16 *A){
    Fp8_clear(&A->x0);
    Fp8_clear(&A->x1);
}
void Fp16_printf(struct Fp16 *A){
    gmp_printf("(%Zd,\n%Zd,\n%Zd,\n%Zd,\n",A->x0.x0.x0.x0.x0,A->x0.x0.x0.x1.x0,A->x0.x0.x1.x0.x0,A->x0.x0.x1.x1.x0);
    gmp_printf("%Zd,\n%Zd,\n%Zd,\n%Zd,\n",A->x0.x1.x0.x0.x0,A->x0.x1.x0.x1.x0,A->x0.x1.x1.x0.x0,A->x0.x1.x1.x1.x0);
    gmp_printf("%Zd,\n%Zd,\n%Zd,\n%Zd,\n",A->x1.x0.x0.x0.x0,A->x1.x0.x0.x1.x0,A->x1.x0.x1.x0.x0,A->x1.x0.x1.x1.x0);
    gmp_printf("%Zd,\n%Zd,\n%Zd,\n%Zd)\n",A->x1.x1.x0.x0.x0,A->x1.x1.x0.x1.x0,A->x1.x1.x1.x0.x0,A->x1.x1.x1.x1.x0);
}
void Fp16_add(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp8_add(&tmp.x0,&A->x0,&B->x0);
    Fp8_add(&tmp.x1,&A->x1,&B->x1);
    
    Fp16_set(ANS,&tmp);
    
    Fp16_clear(&tmp);
}
void Fp16_add_ui(struct Fp16 *ANS,struct Fp16 *A,unsigned long int B){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp8_add_ui(&tmp.x0,&A->x0,B);
    Fp8_add_ui(&tmp.x1,&A->x1,B);
    
    Fp16_set(ANS,&tmp);
    
    Fp16_clear(&tmp);
}
void Fp16_sub(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp8_sub(&tmp.x0,&A->x0,&B->x0);
    Fp8_sub(&tmp.x1,&A->x1,&B->x1);
    
    Fp16_set(ANS,&tmp);
    
    Fp16_clear(&tmp);
}
void Fp16_mul(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B){
    //x^2-v=0
    struct Fp8 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
    Fp8_init(&tmp1);
    Fp8_init(&tmp2);
    Fp8_init(&tmp3);
    Fp8_init(&tmp4);
    Fp8_init(&tmp5);
    Fp8_init(&tmp6);
    
    struct Fp16 t_ans;
    Fp16_init(&t_ans);
    
    Fp8_mul(&tmp1,&A->x0,&B->x0);//a*c
    Fp8_mul(&tmp2,&A->x1,&B->x1);//b*d
    Fp8_mul_v(&tmp3,&tmp2);//b*d*v
    Fp8_add(&t_ans.x0,&tmp1,&tmp3);//a*c+b*d*v
    Fp8_add(&tmp4,&A->x0,&A->x1);//a+b
    Fp8_add(&tmp5,&B->x0,&B->x1);//c+d
    Fp8_mul(&tmp6,&tmp4,&tmp5);//(a+b)(c+d)
    Fp8_sub(&t_ans.x1,&tmp6,&tmp1);
    Fp8_sub(&t_ans.x1,&t_ans.x1,&tmp2);
    
    Fp16_set(ANS,&t_ans);
    
    Fp8_clear(&tmp1);
    Fp8_clear(&tmp2);
    Fp8_clear(&tmp3);
    Fp8_clear(&tmp4);
    Fp8_clear(&tmp5);
    Fp8_clear(&tmp6);
    Fp16_clear(&t_ans);
}
void Fp16_mul_v(struct Fp16 *ANS,struct Fp16 *A){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp8_mul_v(&tmp.x0,&A->x1);
    Fp8_set(&tmp.x1,&A->x0);
    
    Fp16_set(ANS,&tmp);
    Fp16_clear(&tmp);
}
void Fp16_mul_ui(struct Fp16 *ANS,struct Fp16 *A,unsigned long int B){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp8_mul_ui(&tmp.x0,&A->x0,B);
    Fp8_mul_ui(&tmp.x1,&A->x1,B);
    
    Fp16_set(ANS,&tmp);
    
    Fp16_clear(&tmp);
}
void Fp16_mul_Fp(struct Fp16 *ANS,struct Fp16 *A,struct Fp *B){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp8_mul_Fp(&tmp.x0,&A->x0,B);
    Fp8_mul_Fp(&tmp.x1,&A->x1,B);
    
    Fp16_set(ANS,&tmp);
    
    Fp16_clear(&tmp);
}
void Fp16_mul_mpz(struct Fp16 *ANS,struct Fp16 *A,mpz_t B){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp8_mul_mpz(&tmp.x0,&A->x0,B);
    Fp8_mul_mpz(&tmp.x1,&A->x1,B);
    
    Fp16_set(ANS,&tmp);
    Fp16_clear(&tmp);
}
void Fp16_invert(struct Fp16 *ANS,struct Fp16 *A){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    // tmp=A^(q^6)=(x0,-x1)
    Fp8_set(&tmp.x0,&A->x0);
    Fp8_neg(&tmp.x1,&A->x1);
    
    struct Fp8 c,a,b;
    Fp8_init(&c);
    Fp8_init(&a);
    Fp8_init(&b);
    
    Fp8_mul(&a,&A->x0,&A->x0); // a=x0^2
    Fp8_mul(&b,&A->x1,&A->x1); // b=x1^2
    Fp8_mul_v(&b,&b); // b=x1^2*v
    Fp8_sub(&c,&a,&b); // c=x0^2-x1^2*v mod q
    
    Fp8_invert(&c,&c);
    
    // ANS=A^{-1}=(c)^{-1}*A^(p^6) A which c is Fp16-element and tmp is a vector A Fp16
    Fp8_mul(&tmp.x0,&tmp.x0,&c);
    Fp8_mul(&tmp.x1,&tmp.x1,&c);
    
    Fp16_set(ANS,&tmp);
    
    Fp8_clear(&c);
    Fp8_clear(&a);
    Fp8_clear(&b);
    Fp16_clear(&tmp);
}
void Fp16_div(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B){
    struct Fp16 tmp,t_ans;
    Fp16_init(&tmp);
    Fp16_init(&t_ans);
    
    Fp16_invert(&tmp,B);
    Fp16_mul(&t_ans,A,&tmp);
    
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&tmp);
    Fp16_clear(&t_ans);
}
void Fp16_pow(struct Fp16 *ANS,struct Fp16 *A,mpz_t B){
    int i,length;
    length= (int)mpz_sizeinbase(B,2);
    char B_binary[length];
    mpz_get_str(B_binary,2,B);
    struct Fp16 tmp;
    Fp16_init(&tmp);
    Fp16_set(&tmp,A);
    for(i=1;B_binary[i]!='\0';i++){
        Fp16_mul(&tmp,&tmp,&tmp);
        if(B_binary[i]=='1'){
            Fp16_mul(&tmp,&tmp,A);
        }
    }
    Fp16_set(ANS,&tmp);
    Fp16_clear(&tmp);
}
void Fp16_sqrt(struct Fp16 *ANS,struct Fp16 *A){
    struct Fp16 n,y,x,b,t,tmp_Fp4;
    Fp16_init(&n);
    Fp16_init(&y);
    Fp16_init(&x);
    Fp16_init(&b);
    Fp16_init(&t);
    Fp16_init(&tmp_Fp4);
    Fp16_set(&n,A);
    
    mpz_t tmp_mpz,q,e,r,set_1,set_2;
    mpz_init(tmp_mpz);
    mpz_init(q);
    mpz_init(e);
    mpz_init(r);
    mpz_init(set_1);
    mpz_init(set_2);
    mpz_set_ui(set_1,1);
    mpz_set_ui(set_2,2);
    
    while(Fp16_legendre(&n)!=-1){
        Fp16_random(&n);
    }
    mpz_pow_ui(q,params.prime,16);
    mpz_sub_ui(q,q,1);
    mpz_set_ui(e,0);
    while(mpz_odd_p(q)==0){
        mpz_add_ui(e,e,1);
        mpz_div_ui(q,q,2);
    }
    Fp16_pow(&y,&n,q);
    
    mpz_set(r,e);
    
    mpz_sub_ui(tmp_mpz,q,1);
    mpz_div_ui(tmp_mpz,tmp_mpz,2);
    
    Fp16_pow(&x,A,tmp_mpz);
    Fp16_pow(&tmp_Fp4,&x,set_2);
    Fp16_mul(&b,&tmp_Fp4,A);
    Fp16_mul(&x,&x,A);
    
    int m;
    
    while(Fp16_cmp_mpz(&b,set_1)==1){
        m=-1;
        Fp16_set(&tmp_Fp4,&b);
        while(Fp16_cmp_mpz(&tmp_Fp4,set_1)==1){
            m++;
            mpz_pow_ui(tmp_mpz,set_2,m);
            Fp16_pow(&tmp_Fp4,&b,tmp_mpz);
            // Fp16_printf(&tmp_Fp4);
        }
        mpz_sub_ui(tmp_mpz,r,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,params.prime);
        // gmp_printf("%Zd,%Zd,%d\n",tmp_mpz,r,m);
        Fp16_pow(&t,&y,tmp_mpz);
        Fp16_pow(&y,&t,set_2);
        // gmp_printf("%Zd,%Zd,\n",y.x0.x0.x0,y.x0.x1.x0);
        mpz_set_ui(r,m);
        Fp16_mul(&x,&x,&t);
        Fp16_mul(&b,&b,&y);
        
    }
    
    Fp16_set(ANS,&x);
    
    Fp16_clear(&n);
    Fp16_clear(&y);
    Fp16_clear(&x);
    Fp16_clear(&b);
    Fp16_clear(&t);
    Fp16_clear(&tmp_Fp4);
    mpz_clear(tmp_mpz);
    mpz_clear(q);
    mpz_clear(e);
    mpz_clear(r);
    mpz_clear(set_1);
}
int Fp16_legendre(struct Fp16 *a){
    mpz_t i,cmp;
    struct Fp16 tmp;
    Fp16_init(&tmp);
    mpz_init(i);
    mpz_init(cmp);
    mpz_set_ui(cmp,1);
    mpz_pow_ui(i,params.prime,16);
    mpz_sub_ui(i,i,1);
    mpz_tdiv_q_ui(i,i,2);
    Fp16_pow(&tmp,a,i);
    
    if((Fp16_cmp_mpz(&tmp,cmp))==0){
        Fp16_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return 1;
    }else{
        Fp16_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return -1;
    }
}
int Fp16_cmp(struct Fp16 *A,struct Fp16 *B){
    if(Fp8_cmp(&A->x0,&B->x0)==0 && Fp8_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}
int Fp16_cmp_mpz(struct Fp16 *A,mpz_t B){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    if(Fp8_cmp_mpz(&A->x0,B)==0 && Fp8_cmp(&A->x1,&tmp.x1)==0){
        Fp16_clear(&tmp);
        return 0;
    }
    Fp16_clear(&tmp);
    return 1;
}
void Fp16_neg(struct Fp16 *ans,struct Fp16 *a){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    Fp8_neg(&tmp.x0,&a->x0);
    Fp8_neg(&tmp.x1,&a->x1);
    Fp16_set(ans,&tmp);
    Fp16_clear(&tmp);
}

void Fp16_frobenius_map(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 tmp_ans;
    Fp16_init(&tmp_ans);
    Fp16_set(&tmp_ans, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &tmp_ans.x0.x0.x0.x0);
    
    
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &tmp_ans.x0.x0.x0.x1);
    if (mpz_cmp_ui(tmp_ans.x0.x0.x0.x1.x0, 0) != 0) {
        Fp_neg(&ans_tmp8.x0.x0.x0.x1, &ans_tmp8.x0.x0.x0.x1);
    }
    
    Fp_mul(&ans_tmp8.x0.x0.x1.x0, &tmp_ans.x0.x0.x1.x0, &pm1d4);
    
    Fp_mul(&ans_tmp8.x0.x0.x1.x1, &tmp_ans.x0.x0.x1.x1, &m_pm1d4);
    //  Fp_neg(&ans_tmp8.x0.x0.x1.x1, &ans_tmp8.x0.x0.x1.x1);
    
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &tmp_ans.x0.x1.x0.x1, &m_cpm5d8);
    //    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &tmp_ans.x0.x1.x0.x1, &pm5d8);
    //    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &ans_tmp8.x0.x1.x0.x0, &qnr_c1);
    //    Fp_neg(&ans_tmp8.x0.x1.x0.x0, &ans_tmp8.x0.x1.x0.x0);
    
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &tmp_ans.x0.x1.x0.x0, &pm5d8);
    
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &tmp_ans.x0.x1.x1.x1, &m_cpm1d4pm5d8);
    //    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &tmp_ans.x0.x1.x1.x1, &pm1d4);
    //    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &ans_tmp8.x0.x1.x1.x0, &pm5d8);
    //    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &ans_tmp8.x0.x1.x1.x0, &qnr_c1);
    //    Fp_neg(&ans_tmp8.x0.x1.x1.x0, &ans_tmp8.x0.x1.x1.x0);
    
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &tmp_ans.x0.x1.x1.x0, &pm1d4pm5d8);
    //    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &ans_tmp8.x0.x1.x1.x1, &pm5d8);
    
    
    
    //from 9-16
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &tmp_ans.x1.x0.x1.x0, &pm1d4);
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &ans_tmp8.x1.x0.x0.x0, &pm13d16);
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &ans_tmp8.x1.x0.x0.x0, &qnr_c1);
    //     Fp_neg(&ans_tmp8.x1.x0.x0.x0, &ans_tmp8.x1.x0.x0.x0);
    
    //    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &tmp_ans.x1.x0.x1.x1, &m_cpm1d4);
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &tmp_ans.x1.x0.x1.x1, &pm1d4);
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &ans_tmp8.x1.x0.x0.x1, &qnr_c1);
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &ans_tmp8.x1.x0.x0.x1, &pm13d16);
    Fp_neg(&ans_tmp8.x1.x0.x0.x1, &ans_tmp8.x1.x0.x0.x1);
    
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &tmp_ans.x1.x0.x0.x1, &m_cpm13d16);
    //    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &tmp_ans.x1.x0.x0.x1, &pm13d16);
    //    Fp_mul(&ans_tmp8.x1.x0.x1.x0,&ans_tmp8.x1.x0.x1.x0, &qnr_c1);
    //    Fp_neg(&ans_tmp8.x1.x0.x1.x0,&ans_tmp8.x1.x0.x1.x0);
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &tmp_ans.x1.x0.x0.x0, &pm13d16);
    
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &tmp_ans.x1.x1.x1.x1, &m_ccpm1d4pm5d8p13d16);
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &tmp_ans.x1.x1.x1.x0, &cpm1d4pm5d6pm13d16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &tmp_ans.x1.x1.x0.x0, &cpm5d8pm13d16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &tmp_ans.x1.x1.x0.x1, &m_cpm5d8pm13d16);
    
    Fp8_set(&tmp_ans.x0, &ans_tmp8.x0);
    Fp8_set(&tmp_ans.x1, &ans_tmp8.x1);
    
    Fp16_set(ANS,&tmp_ans);
    
    Fp16_clear(&tmp_ans);
}
