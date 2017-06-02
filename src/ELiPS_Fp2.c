//
//  Fp2.c
//  KSS16
//
//  Created by Khandaker Md. Al-Amin on 4/26/17.
//  Copyright © 2017 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_KSS16/ELiPS_Fp2.h>

#pragma mark Fp2 Methods
void Fp2_init(struct Fp2 *A){
    Fp_init(&A->x0);
    Fp_init(&A->x1);
}
void Fp2_set(struct Fp2 *ANS,struct Fp2 *A){
    //    Fp2_printf(A);
    Fp_set(&ANS->x0,&A->x0);
    Fp_set(&ANS->x1,&A->x1);
}
void Fp2_set_ui(struct Fp2 *A,signed long int B){
    Fp_set_ui(&A->x0,B);
    Fp_set_ui(&A->x1,B);
}
void Fp2_random(struct Fp2 *A){
    Fp_random(&A->x0);
    Fp_random(&A->x1);
}
void Fp2_clear(struct Fp2 *A){
    Fp_clear(&A->x0);
    Fp_clear(&A->x1);
}
void Fp2_printf(struct Fp2 *A){
    gmp_printf("(%Zd,\n%Zd)\n",A->x0.x0,A->x1.x0);
}
void Fp2_add(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
    struct Fp2 tmp;
    Fp2_init(&tmp);
    
    Fp_add(&tmp.x0,&A->x0,&B->x0);
    Fp_add(&tmp.x1,&A->x1,&B->x1);
    
    Fp2_set(ANS,&tmp);
    
    Fp2_clear(&tmp);
}
void Fp2_add_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int B){
    struct Fp2 tmp;
    Fp2_init(&tmp);
    
    Fp_add_ui(&tmp.x0,&A->x0,B);
    Fp_add_ui(&tmp.x1,&A->x1,B);
    
    Fp2_set(ANS,&tmp);
    
    Fp2_clear(&tmp);
}
void Fp2_sub(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
    struct Fp2 tmp;
    Fp2_init(&tmp);
    
    Fp_sub(&tmp.x0,&A->x0,&B->x0);
    Fp_sub(&tmp.x1,&A->x1,&B->x1);
    
    Fp2_set(ANS,&tmp);
    
    Fp2_clear(&tmp);
}
void Fp2_mul(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
    struct Fp tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&tmp3);
    Fp_init(&tmp4);
    Fp_init(&tmp5);
    Fp_init(&tmp6);
    
    struct Fp2 t_ans;
    Fp2_init(&t_ans);
    
    Fp_mul(&tmp1,&A->x0,&B->x0);//a*c
    Fp_mul(&tmp2,&A->x1,&B->x1);//b*d
    Fp_mul_ui(&tmp3,&tmp2,c1);//b*d*v
    Fp_add(&t_ans.x0,&tmp1,&tmp3);//a*c+b*d*v
    Fp_add(&tmp4,&A->x0,&A->x1);//a+b
    Fp_add(&tmp5,&B->x0,&B->x1);//c+d
    Fp_mul(&tmp6,&tmp4,&tmp5);//(a+b)(c+d)
    Fp_sub(&t_ans.x1,&tmp6,&tmp1);
    Fp_sub(&t_ans.x1,&t_ans.x1,&tmp2);
    
    Fp2_set(ANS,&t_ans);
    
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&tmp3);
    Fp_clear(&tmp4);
    Fp_clear(&tmp5);
    Fp_clear(&tmp6);
    Fp2_clear(&t_ans);
    
}
void Fp4_mul_basis(struct Fp2 *ANS,struct Fp2 *A){
    //(a,b)(1,1)=(a-b,a+b)
    struct Fp2 tmp;
    Fp2_init(&tmp);
    
    Fp_mul_ui(&tmp.x0, &A->x1, c1);
    Fp_set(&tmp.x1,&A->x0);
    
    Fp2_set(ANS,&tmp);
    Fp2_clear(&tmp);
}
void Fp2_mul_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int B){
    struct Fp2 tmp;
    Fp2_init(&tmp);
    Fp_mul_ui(&tmp.x0,&A->x0,B);
    Fp_mul_ui(&tmp.x1,&A->x1,B);
    Fp2_set(ANS,&tmp);
    Fp2_clear(&tmp);
}
void Fp2_mul_Fp(struct Fp2 *ANS,struct Fp2 *A,struct Fp *B){
    struct Fp2 tmp;
    Fp2_init(&tmp);
    
    Fp_mul(&tmp.x0,&A->x0,B);
    Fp_mul(&tmp.x1,&A->x1,B);
    
    Fp2_set(ANS,&tmp);
    
    Fp2_clear(&tmp);
}
void Fp2_mul_mpz(struct Fp2 *ANS,struct Fp2 *A,mpz_t B){
    struct Fp2 tmp;
    Fp2_init(&tmp);
    
    Fp_mul_mpz(&tmp.x0,&A->x0,B);
    Fp_mul_mpz(&tmp.x1,&A->x1,B);
    
    Fp2_set(ANS,&tmp);
    
    Fp2_clear(&tmp);
}
void Fp2_invert(struct Fp2 *ANS,struct Fp2 *A){
    struct Fp2 tmp;
    Fp2_init(&tmp);
    
    // tmp=A^(q^6)=(x0,-x1)
    Fp_set(&tmp.x0,&A->x0);
    Fp_neg(&tmp.x1,&A->x1);
    
    struct Fp c,a,b;
    Fp_init(&c);
    Fp_init(&a);
    Fp_init(&b);
    
    Fp_mul(&a,&A->x0,&A->x0); // a=x0^2
    Fp_mul(&b,&A->x1,&A->x1); // b=x1^2
    Fp_mul_ui(&b,&b,c1); // b=x1^2*v
    Fp_sub(&c,&a,&b); // c=x0^2-x1^2*v mod q
    
    mpz_invert(c.x0,c.x0,params.prime);
    
    // ANS=A^{-1}=(c)^{-1}*A^(p^6) A which c is Fp2-element and tmp is a vector A Fp2
    Fp_mul(&tmp.x0,&tmp.x0,&c);
    Fp_mul(&tmp.x1,&tmp.x1,&c);
    
    Fp2_set(ANS,&tmp);
    
    Fp_clear(&c);
    Fp_clear(&a);
    Fp_clear(&b);
    Fp2_clear(&tmp);
}
void Fp2_div(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
    struct Fp2 tmp,t_ans;
    Fp2_init(&tmp);
    Fp2_init(&t_ans);
    
    Fp2_invert(&tmp,B);
    Fp2_mul(&t_ans,A,&tmp);
    
    Fp2_set(ANS,&t_ans);
    
    Fp2_clear(&tmp);
    Fp2_clear(&t_ans);
}
void Fp2_pow(struct Fp2 *ANS,struct Fp2 *A,mpz_t B){
    int i;
    int r;//bit数
    r= (int)mpz_sizeinbase(B,2);
    // printf("r= %d\n",r);
    
    struct Fp2 answer_tmp;
    Fp2_init(&answer_tmp);
    Fp2_set(&answer_tmp,A);
    
    struct Fp2 in_tmp;
    Fp2_init(&in_tmp);
    Fp2_set(&in_tmp,A);
    
    for(i=r-2;i>=0;i--){
        if(mpz_tstbit(B,i)==1){
            Fp2_mul(&answer_tmp,&answer_tmp,&answer_tmp);//a*2
            Fp2_mul(&answer_tmp,&answer_tmp,&in_tmp);//*a
        }else{
            Fp2_mul(&answer_tmp,&answer_tmp,&answer_tmp);//a*2
        }
    }
    
    Fp2_set(ANS,&answer_tmp);
    
    Fp2_clear(&answer_tmp);
    Fp2_clear(&in_tmp);
}
void Fp2_sqrt(struct Fp2 *ANS,struct Fp2 *A){
    struct Fp2 n,y,x,b,t,tmp_Fp2;
    Fp2_init(&n);
    Fp2_init(&y);
    Fp2_init(&x);
    Fp2_init(&b);
    Fp2_init(&t);
    Fp2_init(&tmp_Fp2);
    Fp2_set(&n,A);
    
    mpz_t tmp_mpz,q,e,r,set_1,set_2;
    mpz_init(tmp_mpz);
    mpz_init(q);
    mpz_init(e);
    mpz_init(r);
    mpz_init(set_1);
    mpz_init(set_2);
    mpz_set_ui(set_1,1);
    mpz_set_ui(set_2,2);
    
    
    while(Fp2_legendre(&n)!=-1){
        Fp2_random(&n);
    }
    
    mpz_pow_ui(q,params.prime,2);
    mpz_sub_ui(q,q,1);
    mpz_set_ui(e,0);
    
    while(mpz_odd_p(q)==0){
        mpz_add_ui(e,e,1);
        mpz_div_ui(q,q,2);
    }
    
    Fp2_pow(&y,&n,q);
    mpz_set(r,e);
    mpz_sub_ui(tmp_mpz,q,1);
    mpz_div_ui(tmp_mpz,tmp_mpz,2);
    Fp2_pow(&x,A,tmp_mpz);
    Fp2_pow(&tmp_Fp2,&x,set_2);
    Fp2_mul(&b,&tmp_Fp2,A);
    Fp2_mul(&x,&x,A);
    
    int m;
    
    while(Fp2_cmp_mpz(&b,set_1)==1){
        m=-1;
        Fp2_set(&tmp_Fp2,&b);
        while(Fp2_cmp_mpz(&tmp_Fp2,set_1)==1){
            m++;
            mpz_pow_ui(tmp_mpz,set_2,m);
            Fp2_pow(&tmp_Fp2,&b,tmp_mpz);
        }
        mpz_sub_ui(tmp_mpz,r,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,params.prime);
        Fp2_pow(&t,&y,tmp_mpz);
        Fp2_pow(&y,&t,set_2);
        mpz_set_ui(r,m);
        Fp2_mul(&x,&x,&t);
        Fp2_mul(&b,&b,&y);
    }
    
    Fp2_set(ANS,&x);
    
    Fp2_clear(&n);
    Fp2_clear(&y);
    Fp2_clear(&x);
    Fp2_clear(&b);
    Fp2_clear(&t);
    Fp2_clear(&tmp_Fp2);
    mpz_clear(tmp_mpz);
    mpz_clear(q);
    mpz_clear(e);
    mpz_clear(r);
    mpz_clear(set_1);
}
int Fp2_cmp(struct Fp2 *A,struct Fp2 *B){
    if(Fp_cmp(&A->x0,&B->x0)==0 && Fp_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}
int Fp2_cmp_mpz(struct Fp2 *A,mpz_t B){
    struct Fp2 tmp;
    Fp2_init(&tmp);
    if(Fp_cmp_mpz(&A->x0,B)==0 && Fp_cmp(&A->x1,&tmp.x1)==0){
        Fp2_clear(&tmp);
        return 0;
    }
    Fp2_clear(&tmp);
    return 1;
}
int Fp2_legendre(struct Fp2 *a){
    mpz_t i;
    struct Fp2 tmp;
    Fp2_init(&tmp);
    mpz_init(i);
    
    mpz_pow_ui(i,params.prime,2);
    mpz_sub_ui(i,i,1);
    mpz_div_ui(i,i,2);
    
    Fp2_pow(&tmp,a,i);
    
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,1);
    
    if((Fp2_cmp_mpz(&tmp,cmp))==0){
        Fp2_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return 1;
    }else{
        Fp2_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return -1;
    }
}
void Fp2_neg(struct Fp2 *ans,struct Fp2 *a){
    struct Fp2 tmp;
    Fp2_init(&tmp);
    
    Fp_neg(&tmp.x0,&a->x0);
    Fp_neg(&tmp.x1,&a->x1);
    
    Fp2_set(ans,&tmp);
    
    Fp2_clear(&tmp);
}

void Fp2_frobenius_map(struct Fp2 *ANS, struct Fp2 *A){
    struct Fp2 t_ans;
    Fp2_init(&t_ans);
    
    Fp_set(&t_ans.x0,&A->x0);
    if (mpz_cmp_ui(A->x1.x0,0)==0) {
        Fp_set(&t_ans.x1,&A->x1);
    }
    else{
        mpz_sub(t_ans.x1.x0,params.prime,A->x1.x0);
    }
    
    
    Fp2_set(ANS,&t_ans);
    
    Fp2_clear(&t_ans);
}
