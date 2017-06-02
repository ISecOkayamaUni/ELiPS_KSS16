//
//  Fp8.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/26/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_Fp8.h>


#pragma mark Fp8 methods definitions
void Fp8_init(struct Fp8 *A){
    Fp4_init(&A->x0);
    Fp4_init(&A->x1);
}
void Fp8_set(struct Fp8 *ANS,struct Fp8 *A){
    Fp4_set(&ANS->x0,&A->x0);
    Fp4_set(&ANS->x1,&A->x1);
}
void Fp8_set_ui(struct Fp8 *A,signed long int B){
    Fp4_set_ui(&A->x0,B);
    Fp4_set_ui(&A->x1,B);
}
void Fp8_random(struct Fp8 *A){
    Fp4_random(&A->x0);
    Fp4_random(&A->x1);
}
void Fp8_clear(struct Fp8 *A){
    Fp4_clear(&A->x0);
    Fp4_clear(&A->x1);
}
void Fp8_printf(struct Fp8 *A){
    gmp_printf("(%Zd,\n%Zd,\n%Zd,\n%Zd\n",A->x0.x0.x0.x0,A->x0.x0.x1.x0,A->x0.x1.x0.x0,A->x0.x1.x1.x0);
    gmp_printf("%Zd,\n%Zd,\n%Zd,\n%Zd)\n",A->x1.x0.x0.x0,A->x1.x0.x1.x0,A->x1.x1.x0.x0,A->x1.x1.x1.x0);
}
void Fp8_add(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    Fp4_add(&tmp.x0,&A->x0,&B->x0);
    Fp4_add(&tmp.x1,&A->x1,&B->x1);
    
    Fp8_set(ANS,&tmp);
    
    Fp8_clear(&tmp);
}
void Fp8_add_ui(struct Fp8 *ANS,struct Fp8 *A,unsigned long int B){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    Fp4_add_ui(&tmp.x0,&A->x0,B);
    Fp4_add_ui(&tmp.x1,&A->x1,B);
    
    Fp8_set(ANS,&tmp);
    
    Fp8_clear(&tmp);
}
void Fp8_sub(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    Fp4_sub(&tmp.x0,&A->x0,&B->x0);
    Fp4_sub(&tmp.x1,&A->x1,&B->x1);
    
    Fp8_set(ANS,&tmp);
    
    Fp8_clear(&tmp);
}
void Fp8_mul(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B){
    //x^2-v=0
    struct Fp4 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    Fp4_init(&tmp4);
    Fp4_init(&tmp5);
    Fp4_init(&tmp6);
    
    struct Fp8 t_ans;
    Fp8_init(&t_ans);
    
    Fp4_mul(&tmp1,&A->x0,&B->x0);//a*c
    Fp4_mul(&tmp2,&A->x1,&B->x1);//b*d
    Fp4_mul_v(&tmp3,&tmp2);//b*d*v
    Fp4_add(&t_ans.x0,&tmp1,&tmp3);//a*c+b*d*v
    Fp4_add(&tmp4,&A->x0,&A->x1);//a+b
    Fp4_add(&tmp5,&B->x0,&B->x1);//c+d
    Fp4_mul(&tmp6,&tmp4,&tmp5);//(a+b)(c+d)
    Fp4_sub(&t_ans.x1,&tmp6,&tmp1);
    Fp4_sub(&t_ans.x1,&t_ans.x1,&tmp2);
    
    Fp8_set(ANS,&t_ans);
    
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp4_clear(&tmp4);
    Fp4_clear(&tmp5);
    Fp4_clear(&tmp6);
    Fp8_clear(&t_ans);
}
void Fp8_mul_v(struct Fp8 *ANS,struct Fp8 *A){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    Fp4_mul_v(&tmp.x0,&A->x1);
    Fp4_set(&tmp.x1,&A->x0);
    
    Fp8_set(ANS,&tmp);
    Fp8_clear(&tmp);
}
void Fp8_mul_ui(struct Fp8 *ANS,struct Fp8 *A,unsigned long int B){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    Fp4_mul_ui(&tmp.x0,&A->x0,B);
    Fp4_mul_ui(&tmp.x1,&A->x1,B);
    
    Fp8_set(ANS,&tmp);
    
    Fp8_clear(&tmp);
}
void Fp8_mul_Fp(struct Fp8 *ANS,struct Fp8 *A,struct Fp *B){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    Fp4_mul_Fp(&tmp.x0,&A->x0,B);
    Fp4_mul_Fp(&tmp.x1,&A->x1,B);
    
    Fp8_set(ANS,&tmp);
    
    Fp8_clear(&tmp);
}
void Fp8_mul_mpz(struct Fp8 *ANS,struct Fp8 *A,mpz_t B){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    Fp4_mul_mpz(&tmp.x0,&A->x0,B);
    Fp4_mul_mpz(&tmp.x1,&A->x1,B);
    
    Fp8_set(ANS,&tmp);
    
    Fp8_clear(&tmp);
}
void Fp8_invert(struct Fp8 *ANS,struct Fp8 *A){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    // tmp=A^(q^6)=(x0,-x1)
    Fp4_set(&tmp.x0,&A->x0);
    Fp4_neg(&tmp.x1,&A->x1);
    
    struct Fp4 c,a,b;
    Fp4_init(&c);
    Fp4_init(&a);
    Fp4_init(&b);
    
    Fp4_mul(&a,&A->x0,&A->x0); // a=x0^2
    Fp4_mul(&b,&A->x1,&A->x1); // b=x1^2
    Fp4_mul_v(&b,&b); // b=x1^2*v
    Fp4_sub(&c,&a,&b); // c=x0^2-x1^2*v mod q
    
    Fp4_invert(&c,&c);
    
    // ANS=A^{-1}=(c)^{-1}*A^(p^6) A which c is Fp8-element and tmp is a vector A Fp8
    Fp4_mul(&tmp.x0,&tmp.x0,&c);
    Fp4_mul(&tmp.x1,&tmp.x1,&c);
    
    Fp8_set(ANS,&tmp);
    
    Fp4_clear(&c);
    Fp4_clear(&a);
    Fp4_clear(&b);
    Fp8_clear(&tmp);
}
void Fp8_div(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B){
    struct Fp8 tmp,t_ans;
    Fp8_init(&tmp);
    Fp8_init(&t_ans);
    
    Fp8_invert(&tmp,B);
    Fp8_mul(&t_ans,A,&tmp);
    
    Fp8_set(ANS,&t_ans);
    
    Fp8_clear(&tmp);
    Fp8_clear(&t_ans);
}
void Fp8_pow(struct Fp8 *ANS,struct Fp8 *A,mpz_t B){
    int i,length;
    length= (int)mpz_sizeinbase(B,2);
    char B_binary[length];
    mpz_get_str(B_binary,2,B);
    struct Fp8 tmp;
    Fp8_init(&tmp);
    Fp8_set(&tmp,A);
    for(i=1;B_binary[i]!='\0';i++){
        Fp8_mul(&tmp,&tmp,&tmp);
        if(B_binary[i]=='1'){
            Fp8_mul(&tmp,&tmp,A);
        }
    }
    Fp8_set(ANS,&tmp);
    Fp8_clear(&tmp);
}
void Fp8_sqrt(struct Fp8 *ANS,struct Fp8 *A){
    struct Fp8 n,y,x,b,t,tmp_Fp4;
    Fp8_init(&n);
    Fp8_init(&y);
    Fp8_init(&x);
    Fp8_init(&b);
    Fp8_init(&t);
    Fp8_init(&tmp_Fp4);
    Fp8_set(&n,A);
    
    mpz_t tmp_mpz,q,e,r,set_1,set_2;
    mpz_init(tmp_mpz);
    mpz_init(q);
    mpz_init(e);
    mpz_init(r);
    mpz_init(set_1);
    mpz_init(set_2);
    mpz_set_ui(set_1,1);
    mpz_set_ui(set_2,2);
    
    while(Fp8_legendre(&n)!=-1){
        Fp8_random(&n);
    }
    mpz_pow_ui(q,params.prime,12);
    mpz_sub_ui(q,q,1);
    mpz_set_ui(e,0);
    while(mpz_odd_p(q)==0){
        mpz_add_ui(e,e,1);
        mpz_div_ui(q,q,2);
    }
    Fp8_pow(&y,&n,q);
    
    mpz_set(r,e);
    
    mpz_sub_ui(tmp_mpz,q,1);
    mpz_div_ui(tmp_mpz,tmp_mpz,2);
    
    Fp8_pow(&x,A,tmp_mpz);
    Fp8_pow(&tmp_Fp4,&x,set_2);
    Fp8_mul(&b,&tmp_Fp4,A);
    Fp8_mul(&x,&x,A);
    
    int m;
    
    while(Fp8_cmp_mpz(&b,set_1)==1){
        m=-1;
        Fp8_set(&tmp_Fp4,&b);
        while(Fp8_cmp_mpz(&tmp_Fp4,set_1)==1){
            m++;
            mpz_pow_ui(tmp_mpz,set_2,m);
            Fp8_pow(&tmp_Fp4,&b,tmp_mpz);
        }
        mpz_sub_ui(tmp_mpz,r,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,params.prime);
        // gmp_printf("%Zd,%Zd,%d\n",tmp_mpz,r,m);
        Fp8_pow(&t,&y,tmp_mpz);
        Fp8_pow(&y,&t,set_2);
        // gmp_printf("%Zd,%Zd,\n",y.x0.x0.x0,y.x0.x1.x0);
        mpz_set_ui(r,m);
        Fp8_mul(&x,&x,&t);
        Fp8_mul(&b,&b,&y);
    }
    
    Fp8_set(ANS,&x);
    
    Fp8_clear(&n);
    Fp8_clear(&y);
    Fp8_clear(&x);
    Fp8_clear(&b);
    Fp8_clear(&t);
    Fp8_clear(&tmp_Fp4);
    mpz_clear(tmp_mpz);
    mpz_clear(q);
    mpz_clear(e);
    mpz_clear(r);
    mpz_clear(set_1);
}
int Fp8_legendre(struct Fp8 *a){
    mpz_t i,cmp;
    struct Fp8 tmp;
    Fp8_init(&tmp);
    mpz_init(i);
    mpz_init(cmp);
    mpz_set_ui(cmp,1);
    mpz_pow_ui(i,params.prime,8);
    mpz_sub_ui(i,i,1);
    mpz_tdiv_q_ui(i,i,2);
    Fp8_pow(&tmp,a,i);
    
    if((Fp8_cmp_mpz(&tmp,cmp))==0){
        Fp8_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return 1;
    }else{
        Fp8_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return -1;
    }
}
int Fp8_cmp(struct Fp8 *A,struct Fp8 *B){
    if(Fp4_cmp(&A->x0,&B->x0)==0 && Fp4_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}
int Fp8_cmp_mpz(struct Fp8 *A,mpz_t B){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    if(Fp4_cmp_mpz(&A->x0,B)==0 && Fp4_cmp(&A->x1,&tmp.x1)==0){
        Fp8_clear(&tmp);
        return 0;
    }
    Fp8_clear(&tmp);
    return 1;
}
void Fp8_neg(struct Fp8 *ans,struct Fp8 *a){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    Fp4_neg(&tmp.x0,&a->x0);
    Fp4_neg(&tmp.x1,&a->x1);
    Fp8_set(ans,&tmp);
    Fp8_clear(&tmp);
}

void Fp8_frobenius_map(struct Fp8 *ANS, struct Fp8 *A){
    struct Fp8 tmp_ans;
    struct Fp4 ans_tmp4;
    Fp4_init(&ans_tmp4);
    Fp8_init(&tmp_ans);
    
    Fp4_frobenius_map(&tmp_ans.x0,&A->x0);
    Fp4_frobenius_map(&tmp_ans.x1,&A->x1);
    Fp4_mul_Fp(&tmp_ans.x1,&tmp_ans.x1,&pm5d8);
    
    Fp4_mul_basis(&tmp_ans.x1.x0, &tmp_ans.x1.x0);
    Fp4_mul_basis(&tmp_ans.x1.x1, &tmp_ans.x1.x1);
    
    Fp8_set(ANS,&tmp_ans);
    
    //    Fp_clear(&pm5d8);
    Fp4_clear(&ans_tmp4);
    //    Fp_clear(&set_c1);
    Fp8_clear(&tmp_ans);
}
