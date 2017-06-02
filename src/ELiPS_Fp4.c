//
//  Fp4.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/26/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_Fp4.h>

struct Fp4 z_inv2;

#pragma mark Fp4 methods

void Fp4_init(struct Fp4 *A){
    Fp2_init(&A->x0);
    Fp2_init(&A->x1);
}
void Fp4_set(struct Fp4 *ANS,struct Fp4 *A){
    Fp2_set(&ANS->x0,&A->x0);
    Fp2_set(&ANS->x1,&A->x1);
}
void Fp4_set_ui(struct Fp4 *A,signed long int B){
    Fp2_set_ui(&A->x0,B);
    Fp2_set_ui(&A->x1,B);
}
void Fp4_random(struct Fp4 *A){
    Fp2_random(&A->x0);
    Fp2_random(&A->x1);
}
void Fp4_clear(struct Fp4 *A){
    Fp2_clear(&A->x0);
    Fp2_clear(&A->x1);
}
void Fp4_printf(struct Fp4 *A){
    gmp_printf("(%Zd,\n%Zd,\n%Zd,\n%Zd)\n",A->x0.x0.x0,A->x0.x1.x0,A->x1.x0.x0,A->x1.x1.x0);
}
void Fp4_add(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    
    Fp2_add(&tmp.x0,&A->x0,&B->x0);
    Fp2_add(&tmp.x1,&A->x1,&B->x1);
    
    Fp4_set(ANS,&tmp);
    
    Fp4_clear(&tmp);
}
void Fp4_add_ui(struct Fp4 *ANS,struct Fp4 *A,unsigned long int B){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    
    Fp2_add_ui(&tmp.x0,&A->x0,B);
    Fp2_add_ui(&tmp.x1,&A->x1,B);
    
    Fp4_set(ANS,&tmp);
    
    Fp4_clear(&tmp);
}
void Fp4_sub(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    
    Fp2_sub(&tmp.x0,&A->x0,&B->x0);
    Fp2_sub(&tmp.x1,&A->x1,&B->x1);
    
    Fp4_set(ANS,&tmp);
    
    Fp4_clear(&tmp);
}
void Fp4_mul(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B){
    //x^2-v=0
    struct Fp2 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
    Fp2_init(&tmp1);
    Fp2_init(&tmp2);
    Fp2_init(&tmp3);
    Fp2_init(&tmp4);
    Fp2_init(&tmp5);
    Fp2_init(&tmp6);
    
    struct Fp4 t_ans;
    Fp4_init(&t_ans);
    
    Fp2_mul(&tmp1,&A->x0,&B->x0);//a*c
    Fp2_mul(&tmp2,&A->x1,&B->x1);//b*d
    Fp4_mul_basis(&tmp3,&tmp2);//b*d*v
    Fp2_add(&t_ans.x0,&tmp1,&tmp3);//a*c+b*d*v
    Fp2_add(&tmp4,&A->x0,&A->x1);//a+b
    Fp2_add(&tmp5,&B->x0,&B->x1);//c+d
    Fp2_mul(&tmp6,&tmp4,&tmp5);//(a+b)(c+d)
    Fp2_sub(&t_ans.x1,&tmp6,&tmp1);
    Fp2_sub(&t_ans.x1,&t_ans.x1,&tmp2);
    
    Fp4_set(ANS,&t_ans);
    
    Fp2_clear(&tmp1);
    Fp2_clear(&tmp2);
    Fp2_clear(&tmp3);
    Fp2_clear(&tmp4);
    Fp2_clear(&tmp5);
    Fp2_clear(&tmp6);
    Fp4_clear(&t_ans);
}
void Fp4_mul_v(struct Fp4 *ANS,struct Fp4 *A){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    
    Fp4_mul_basis(&tmp.x0,&A->x1);
    Fp2_set(&tmp.x1,&A->x0);
    
    Fp4_set(ANS,&tmp);
    Fp4_clear(&tmp);
}
void Fp4_mul_ui(struct Fp4 *ANS,struct Fp4 *A,unsigned long int B){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    
    Fp2_mul_ui(&tmp.x0,&A->x0,B);
    Fp2_mul_ui(&tmp.x1,&A->x1,B);
    
    Fp4_set(ANS,&tmp);
    
    Fp4_clear(&tmp);
}
void Fp4_mul_mpz(struct Fp4 *ANS,struct Fp4 *A,mpz_t B){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    
    Fp2_mul_mpz(&tmp.x0,&A->x0,B);
    Fp2_mul_mpz(&tmp.x1,&A->x1,B);
    
    Fp4_set(ANS,&tmp);
    
    Fp4_clear(&tmp);
}
void Fp4_mul_Fp(struct Fp4 *ANS,struct Fp4 *A,struct Fp *B){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    
    Fp2_mul_Fp(&tmp.x0,&A->x0,B);
    Fp2_mul_Fp(&tmp.x1,&A->x1,B);
    
    Fp4_set(ANS,&tmp);
    
    Fp4_clear(&tmp);
}
void Fp4_invert(struct Fp4 *ANS,struct Fp4 *A){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    
    // tmp=A^(q^6)=(x0,-x1)
    Fp2_set(&tmp.x0,&A->x0);
    Fp2_neg(&tmp.x1,&A->x1);
    
    struct Fp2 c,a,b;
    Fp2_init(&c);
    Fp2_init(&a);
    Fp2_init(&b);
    
    Fp2_mul(&a,&A->x0,&A->x0); // a=x0^2
    Fp2_mul(&b,&A->x1,&A->x1); // b=x1^2
    Fp4_mul_basis(&b,&b); // b=x1^2*v
    Fp2_sub(&c,&a,&b); // c=x0^2-x1^2*v mod q
    
    Fp2_invert(&c,&c);
    
    // ANS=A^{-1}=(c)^{-1}*A^(p^6) A which c is Fp4-element and tmp is a vector A Fp4
    Fp2_mul(&tmp.x0,&tmp.x0,&c);
    Fp2_mul(&tmp.x1,&tmp.x1,&c);
    
    Fp4_set(ANS,&tmp);
    
    Fp2_clear(&c);
    Fp2_clear(&a);
    Fp2_clear(&b);
    Fp4_clear(&tmp);
}
void Fp4_div(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B){
    struct Fp4 tmp,t_ans;
    Fp4_init(&tmp);
    Fp4_init(&t_ans);
    
    Fp4_invert(&tmp,B);
    Fp4_mul(&t_ans,A,&tmp);
    
    Fp4_set(ANS,&t_ans);
    
    Fp4_clear(&tmp);
    Fp4_clear(&t_ans);
}
void Fp4_pow(struct Fp4 *ANS,struct Fp4 *A,mpz_t B){
    int i,length;
    length= (int)mpz_sizeinbase(B,2);
    char B_binary[length];
    mpz_get_str(B_binary,2,B);
    struct Fp4 tmp;
    Fp4_init(&tmp);
    Fp4_set(&tmp,A);
    for(i=1;B_binary[i]!='\0';i++){
        Fp4_mul(&tmp,&tmp,&tmp);
        if(B_binary[i]=='1'){
            Fp4_mul(&tmp,&tmp,A);
        }
    }
    Fp4_set(ANS,&tmp);
    Fp4_clear(&tmp);
}
void Fp4_sqrt(struct Fp4 *ANS,struct Fp4 *A){
    struct Fp4 n,y,x,b,t,tmp_Fp4;
    Fp4_init(&n);
    Fp4_init(&y);
    Fp4_init(&x);
    Fp4_init(&b);
    Fp4_init(&t);
    Fp4_init(&tmp_Fp4);
    Fp4_set(&n,A);
    
    mpz_t tmp_mpz,q,e,r,set_1,set_2;
    mpz_init(tmp_mpz);
    mpz_init(q);
    mpz_init(e);
    mpz_init(r);
    mpz_init(set_1);
    mpz_init(set_2);
    mpz_set_ui(set_1,1);
    mpz_set_ui(set_2,2);
    
    while(Fp4_legendre(&n)!=-1){
        Fp4_random(&n);
    }
    mpz_pow_ui(q,params.prime,12);
    mpz_sub_ui(q,q,1);
    mpz_set_ui(e,0);
    while(mpz_odd_p(q)==0){
        mpz_add_ui(e,e,1);
        mpz_div_ui(q,q,2);
    }
    Fp4_pow(&y,&n,q);
    
    mpz_set(r,e);
    
    mpz_sub_ui(tmp_mpz,q,1);
    mpz_div_ui(tmp_mpz,tmp_mpz,2);
    
    Fp4_pow(&x,A,tmp_mpz);
    Fp4_pow(&tmp_Fp4,&x,set_2);
    Fp4_mul(&b,&tmp_Fp4,A);
    Fp4_mul(&x,&x,A);
    
    int m;
    
    while(Fp4_cmp_mpz(&b,set_1)==1){
        m=-1;
        Fp4_set(&tmp_Fp4,&b);
        while(Fp4_cmp_mpz(&tmp_Fp4,set_1)==1){
            m++;
            mpz_pow_ui(tmp_mpz,set_2,m);
            Fp4_pow(&tmp_Fp4,&b,tmp_mpz);
        }
        mpz_sub_ui(tmp_mpz,r,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,params.prime);
        // gmp_printf("%Zd,%Zd,%d\n",tmp_mpz,r,m);
        Fp4_pow(&t,&y,tmp_mpz);
        Fp4_pow(&y,&t,set_2);
        // gmp_printf("%Zd,%Zd,\n",y.x0.x0.x0,y.x0.x1.x0);
        mpz_set_ui(r,m);
        Fp4_mul(&x,&x,&t);
        Fp4_mul(&b,&b,&y);
    }
    
    Fp4_set(ANS,&x);
    
    Fp4_clear(&n);
    Fp4_clear(&y);
    Fp4_clear(&x);
    Fp4_clear(&b);
    Fp4_clear(&t);
    Fp4_clear(&tmp_Fp4);
    mpz_clear(tmp_mpz);
    mpz_clear(q);
    mpz_clear(e);
    mpz_clear(r);
    mpz_clear(set_1);
}
int Fp4_legendre(struct Fp4 *a){
    mpz_t i,cmp;
    struct Fp4 tmp;
    Fp4_init(&tmp);
    mpz_init(i);
    mpz_init(cmp);
    mpz_set_ui(cmp,1);
    mpz_pow_ui(i,params.prime,4);
    mpz_sub_ui(i,i,1);
    mpz_tdiv_q_ui(i,i,2);
    Fp4_pow(&tmp,a,i);
    
    if((Fp4_cmp_mpz(&tmp,cmp))==0){
        Fp4_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return 1;
    }else{
        Fp4_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return -1;
    }
}
int Fp4_cmp(struct Fp4 *A,struct Fp4 *B){
    if(Fp2_cmp(&A->x0,&B->x0)==0 && Fp2_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}
int Fp4_cmp_mpz(struct Fp4 *A,mpz_t B){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    if(Fp2_cmp_mpz(&A->x0,B)==0 && Fp2_cmp(&A->x1,&tmp.x1)==0){
        Fp4_clear(&tmp);
        return 0;
    }
    Fp4_clear(&tmp);
    return 1;
}
void Fp4_neg(struct Fp4 *ans,struct Fp4 *a){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    Fp2_neg(&tmp.x0,&a->x0);
    Fp2_neg(&tmp.x1,&a->x1);
    Fp4_set(ans,&tmp);
    Fp4_clear(&tmp);
}

void Fp4_mul_beta_inv(struct Fp4 *ANS)
{
    struct Fp4 tmp;
    Fp4_init(&tmp);
    Fp4_set_ui(&tmp, 0);
    
    mpz_set(tmp.x1.x1.x0,kss_curve_const.a);
    
    mpz_t c_inv,c;
    mpz_init(c_inv);
    mpz_init(c);
    mpz_set_ui(c,c1);
    mpz_invert(c_inv,c,params.prime);
    mpz_mul(tmp.x1.x1.x0,tmp.x1.x1.x0,c_inv);
    
    Fp4_set(ANS, &tmp);
    
    mpz_clear(c);
    mpz_clear(c_inv);
    Fp4_clear(&tmp);
}

void Fp4_frobenius_map(struct Fp4 *ANS, struct Fp4 *A){
    struct Fp4 t_ans;
    Fp4_init(&t_ans);
    
    Fp2_frobenius_map(&t_ans.x0,&A->x0);
    Fp2_frobenius_map(&t_ans.x1,&A->x1);
    Fp2_mul_Fp(&t_ans.x1,&t_ans.x1,&pm1d4);
    
    Fp4_set(ANS,&t_ans);
    
    Fp4_clear(&t_ans);
}
