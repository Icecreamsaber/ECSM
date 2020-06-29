#include <gmpxx.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <math.h>
using namespace std;

static inline void Complete_addition(mpz_class *p, mpz_class *a, mpz_class *b3, mpz_class *X1, mpz_class *Y1, mpz_class *Z1, mpz_class *X2, mpz_class *Y2, mpz_class *Z2, mpz_class XYZ[3]);
static inline void Affine_add(mpz_class *P, mpz_class *Q, mpz_class *p);
static inline void Affine_dbl(mpz_class *P, mpz_class *r, mpz_class *a, mpz_class *p);
static inline void Affine_dbl_change(mpz_class *P, mpz_class *a, mpz_class *p);
static inline void Extended_Affine(mpz_class *P, mpz_class *Q, mpz_class *p);
static inline void Affine_DQ(mpz_class *P, mpz_class *P2, mpz_class *P4, mpz_class *a, mpz_class *p);
static inline void Affine_PQDQ(mpz_class P[2], mpz_class Q[2], mpz_class *a, mpz_class *p);
static inline void Affine_DPQ(mpz_class P[2], mpz_class Q[2], mpz_class *p);
static inline void Affine_DQ_LR(mpz_class P[2], mpz_class *a, mpz_class *p);
static inline void Affine_PQRT(mpz_class P[2], mpz_class Q[2], mpz_class R[2], mpz_class T[2], mpz_class *p);
static inline void Affine_RLML(mpz_class P[2], mpz_class Q[2], mpz_class R[2], mpz_class T[2], mpz_class *a, mpz_class *p);
static inline void Affine_PQR(mpz_class P[2], mpz_class Q[2], mpz_class R[2], mpz_class *p);
static inline void Affine_ODPQ(mpz_class P[2], mpz_class Q[2], mpz_class *a, mpz_class *p);
static inline void Affine_QPQ(mpz_class P[2], mpz_class Q[2], mpz_class *a, mpz_class *p);
static inline void Affine_add_(mpz_class *P, mpz_class *Q, mpz_class *R, mpz_class *p);

static inline void Complete_addition(mpz_class *p, mpz_class *a, mpz_class *b3, mpz_class *X1, mpz_class *Y1, mpz_class *Z1, mpz_class *X2, mpz_class *Y2, mpz_class *Z2, mpz_class XYZ[3]){
	mpz_class t0, t1, t2, t3, t4, t5;
	t0 = *X1 * *X2 % *p; t1 = *Y1 * *Y2 % *p; t2 = *Z1 * *Z2 % *p;
	t3 = *X1 + *Y1; t4 = *X2 + *Y2; t3 = t3 * t4 % *p;
	t4 = t0 + t1; t3 = t3 - t4; t4 = *X1 + *Z1;
	t5 = *X2 + *Z2; t4 = t4 * t5 % *p; t5 = t0 + t2;
	t4 = t4 - t5; t5 = *Y1 + *Z1;XYZ[0] = *Y2 + *Z2;
	t5 = t5 * XYZ[0] % *p; XYZ[0] = t1 + t2; t5 = t5 - XYZ[0];
	XYZ[2] = *a * t4 % *p; XYZ[0] = *b3 * t2 % *p; XYZ[2] = XYZ[0] + XYZ[2];
	XYZ[0] = t1 - XYZ[2]; XYZ[2] = t1 + XYZ[2]; XYZ[1] = XYZ[0] * XYZ[2] % *p;
	t1 = t0 + t0; t1 = t1 + t0; t2 = *a * t2 % *p;
	t4 = *b3 * t4 % *p; t1 = t1 + t2; t2 = t0 - t2;
	t2 = *a * t2 % *p; t4 = t4 + t2; t0 = t1 * t4 % *p;
	XYZ[1] = (XYZ[1] + t0) % *p; t0 = t5 * t4 % *p; XYZ[0] = t3 * XYZ[0] % *p;
	XYZ[0] = (XYZ[0] - t0) % *p; t0 = t3 * t1 % *p; XYZ[2] = t5 * XYZ[2] % *p;
	XYZ[2] = (XYZ[2] + t0) % *p; 
	return;
}

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
void Joye_RL(mpz_class P[3], string scalar, mpz_class A[3], mpz_class *a, mpz_class *b3, mpz_class *p){
	//Initialization
	mpz_class *R[2];
	R[0] = new mpz_class[3];
	R[1] = new mpz_class[3];
	R[0][0] = mpz_class("0"); R[0][1] = mpz_class("1"); R[0][2] = mpz_class("0");
	R[1][0] = mpz_class("0"); R[1][1] = mpz_class("1"); R[1][2] = mpz_class("0");
	for (int i = 0; i < 3; ++i) A[i] = P[i];
	//cout<<R[0][0]<<' '<<R[0][1]<<' '<<R[0][2]<<endl;
	//cout<<R[1][0]<<' '<<R[1][1]<<' '<<R[1][2]<<endl; 
	//cout<<A[0]<<' '<<A[1]<<' '<<A[2]<<endl;  	
	//Main Loop
	for (int i = scalar.size()-1; i > 1; --i){
		Complete_addition(p, a, b3, &R[scalar[i]-'0'][0], &R[scalar[i]-'0'][1], &R[scalar[i]-'0'][2], &A[0], &A[1], &A[2], R[scalar[i]-'0']);
		Complete_addition(p, a, b3, &A[0], &A[1], &A[2], &A[0], &A[1], &A[2], A);
	}
	//Aggregation and Final correction
	A[0] *= (scalar[1]-'0'-1);
	A[2] *= (scalar[1]-'0'-1);
	Complete_addition(p, a, b3, &A[0], &A[1], &A[2], &R[0][0], &R[0][1], &R[0][2], A);
	Complete_addition(p, a, b3, &R[1][0], &R[1][1], &R[1][2], &R[1][0], &R[1][1], &R[1][2], R[1]);
	Complete_addition(p, a, b3, &A[0], &A[1], &A[2], &R[1][0], &R[1][1], &R[1][2], A);
	Complete_addition(p, a, b3, &P[0], &P[1], &P[2], &A[0], &A[1], &A[2], A);
	A[1] = pow(-1, int(scalar[0])) * A[1];
	return; 
}
#pragma GCC pop_options

void generalization(mpz_class P, mpz_class R[3]){
	mpz_invert(R[2].get_mpz_t(), R[2].get_mpz_t(), P.get_mpz_t());
	for (int i = 0; i < 2; ++i){
		R[i] = R[i] * R[2] % P;
		while(R[i] < 0) R[i] += P;
	}
	return;
}

static inline void Affine_add(mpz_class *P, mpz_class *Q, mpz_class *p){
	mpz_class t0;
	t0 = Q[0]-P[0];
	mpz_invert(t0.get_mpz_t(), t0.get_mpz_t(), (*p).get_mpz_t());
	Q[1] -= P[1];
	t0 = t0 * Q[1] % (*p);
	Q[1] = (t0 * t0 - P[0] - Q[0]) % (*p);
	Q[0] = ((P[0] - Q[1]) * t0 - P[1]) % (*p);
	t0 = Q[1];
	Q[1] = Q[0];
	Q[0] = t0;
	return;
}

static inline void Affine_add_(mpz_class *P, mpz_class *Q, mpz_class *R, mpz_class *p){
	mpz_class t0;
	t0 = Q[0]-P[0];
	mpz_invert(t0.get_mpz_t(), t0.get_mpz_t(), (*p).get_mpz_t());
	R[1] = Q[1]-P[1];
	t0 = t0 * R[1] % (*p);
	R[0] = (t0 * t0 - P[0] - Q[0]) % (*p);
	R[1] = ((P[0] - R[0]) * t0 - P[1]) % (*p);
	return;
}

static inline void Affine_dbl(mpz_class *P, mpz_class *r, mpz_class *a, mpz_class *p){
	mpz_class t0;
	t0 = (3 * (P[0] * P[0]) + (*a)) % (*p);//
	r[0] = 2 * P[1];
	mpz_invert(r[0].get_mpz_t(), r[0].get_mpz_t(), (*p).get_mpz_t());
	t0 = t0 * r[0] % (*p);
	r[0] = (t0 * t0 - 2 * P[0]) % (*p);
	r[1] = ((P[0] - r[0]) * t0  - P[1]) % (*p);
	return;
}

static inline void Affine_dbl_change(mpz_class *P, mpz_class *a, mpz_class *p){
	mpz_class t0, t1, t2;
	t0 = (3 * (P[0] * P[0]) + (*a)) % (*p);
	t1 = 2 * P[1];
	mpz_invert(t1.get_mpz_t(), t1.get_mpz_t(), (*p).get_mpz_t());
	t0 = t0 * t1 % (*p);
	t1 = (t0 * t0 - 2 * P[0]) % (*p);
	t2 = ((P[0] - t1) * t0  - P[1]) % (*p);
	P[0] = t1;
	P[1] = t2;
	return;
}

static inline void Extended_Affine(mpz_class *P, mpz_class *Q, mpz_class *p){
	mpz_class t0, t1, t2;
	t0 = Q[0]-P[0];
	Q[0] = t0;
	mpz_invert(t0.get_mpz_t(), t0.get_mpz_t(), (*p).get_mpz_t());
	Q[1] -= P[1];
	t1 = (Q[0] + 2 * P[0]) * Q[0] % (*p);
	Q[0] = P[1] * Q[0] % (*p);
	t2 = (Q[1] * Q[1] % (*p) * t0 - t1) % (*p) * t0 % (*p); 
	t1 = ((P[0] - t2) * Q[1] - Q[0]) % (*p) * t0 % (*p);
	Q[0] = t2;
	Q[1] = t1;
	return;
}

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
void New_2ary_RL(mpz_class P[2], string scalar, mpz_class A[2], mpz_class *a, mpz_class *p){
	//Initialization
	mpz_class *R[2];
	R[0] = new mpz_class[2];
	R[1] = new mpz_class[2];
	R[0][0] = P[0]; R[0][1] = -P[1];
	R[1][0] = P[0]; R[1][1] = P[1];
	Affine_dbl(P, A, a, p);
	Affine_add(A, R[scalar[scalar.size()-1]-'0'], p);
	//Main Loop
	for (int i = scalar.size()-2; i >= 1; --i){
		Affine_add(A, R[scalar[i]-'0'], p);
		Affine_dbl_change(A, a, p);
	}
	//Final Correction
	P[1] = -P[1];
	Affine_add(P, R[scalar[scalar.size()-1]-'0'], p);
	A[1] = -A[1];
	Affine_add(R[0], A, p);
	Affine_dbl_change(R[1], a, p);
	Extended_Affine(R[1], A, p);
	A[1] = pow(-1, int(scalar[0])) * A[1];
	P[1] = -P[1];
	return;
}
#pragma GCC pop_options

void generalization_aff(mpz_class P[2], mpz_class p){
	while(P[0] < 0) P[0] += p;
	while(P[1] < 0) P[1] += p;
	return;
}

/*
static inline void Affine_DQ(mpz_class P[2], mpz_class P2[2], mpz_class P4[2], mpz_class *a, mpz_class *p){
	mpz_class t0, t1, t2, t3, t4;
	t0 = P[0] * P[0] % (*p); t1 = 3 * t0 + (*a); t2 = 2 * P[1] * P[1] % (*p);
	t3 = t2 * t2 % (*p); t4 = (P[0] + t2) * (P[0] + t2) % (*p) - t0 - t3;
	t0 = (3 * t4 - t1 * t1 % (*p)) * t1 % (*p) - 2 * t3; t2 = 2 * P[1] * t0 % (*p);
	mpz_invert(t4.get_mpz_t(), t2.get_mpz_t(), (*p).get_mpz_t());
	t2 = t0 * t1 % (*p) * t4 % (*p); P2[0] = (t2 * t2 - 2 * P[0]) % (*p); 
	P2[1] = ((P[0] - P2[0]) * t2 - P[1]) % (*p); t2 = 3 * P2[0] * P2[0] % (*p) + (*a);
	t0 = 2 * t2 * t3 % (*p) * t4 % (*p); P4[0] = (t0 * t0 - 2 * P2[0]) % (*p);
	P4[1] = ((P2[0] - P4[0]) * t0 - P2[1]) % (*p);
}
*/
static inline void Affine_DQ(mpz_class P[2], mpz_class P2[2], mpz_class P4[2], mpz_class *a, mpz_class *p){
	mpz_class t0, t1, t2, t3;
	t0 = P[0] * P[0] % (*p); t1 = P[1] * P[1] % (*p); t1 = 2 * t1; t2 = t1 * t1 % (*p);
	t1 = t1 + P[0]; t1 = t1 * t1  % (*p); t1 = t1 - t0; t1 = t1 - t2; t1 = 3 * t1;
	t0 = 3 * t0; t0 = t0 + (*a); t3 = t0 * t0 % (*p); t1 = t1 - t3; t1 = t1 * t0 % (*p);
	t2 = 2 * t2; t1 = t1 - t2; t3 = 2 * t1; t3 = t3 * P[1] % (*p);
	mpz_invert(t3.get_mpz_t(), t3.get_mpz_t(), (*p).get_mpz_t());
	t0 = t0 * t1 % (*p); t0 = t0 * t3 % (*p); P2[0] = t0 * t0 % (*p); P2[0] = P2[0] - P[0]; P2[0] = (P2[0] - P[0]) % (*p);
	P2[1] = P[0] - P2[0]; P2[1] = P2[1] * t0 % (*p); P2[1] = (P2[1] - P[1]) % (*p); t3 = t2 * t3 % (*p);
	t0 = P2[0] * P2[0] % (*p); t0 = 3 * t0; t0 = t0 + (*a); t0 = t0 * t3 % (*p); P4[0] = t0 * t0 % (*p); P4[0] = P4[0] - P2[0];
	P4[0] = (P4[0] - P2[0]) % (*p); P4[1] = P2[0] - P4[0]; P4[1] = P4[1] * t0 % (*p); P4[1] = (P4[1] - P2[1]) % (*p);
	return;
}

static inline void Affine_DQ_LR(mpz_class P[2], mpz_class *a, mpz_class *p){
	mpz_class t0, t1, t2, t3;
	t0 = P[0] * P[0] % (*p); t1 = P[1] * P[1] % (*p); t1 = 2 * t1; t2 = t1 * t1 % (*p);
	t1 = t1 + P[0]; t1 = t1 * t1  % (*p); t1 = t1 - t0; t1 = t1 - t2; t1 = 3 * t1;
	t0 = 3 * t0; t0 = t0 + (*a); t3 = t0 * t0 % (*p); t1 = t1 - t3; t1 = t1 * t0 % (*p);
	t2 = 2 * t2; t1 = t1 - t2; t3 = 2 * t1; t3 = t3 * P[1] % (*p);
	mpz_invert(t3.get_mpz_t(), t3.get_mpz_t(), (*p).get_mpz_t());
	t0 = t0 * t1 % (*p); t0 = t0 * t3 % (*p); t1 = t0 * t0 % (*p); t1 = t1 - P[0]; t1 = t1 - P[0];
	P[0] = P[0] - t1; P[0] = P[0] * t0 % (*p); P[0] = P[0] - P[1]; t3 = t2 * t3 % (*p);
	t0 = t1 * t1 % (*p); t0 = 3 * t0; t0 = t0 + (*a); t0 = t0 * t3 % (*p); P[1] = t0 * t0 % (*p); P[1] = P[1] - t1;
	P[1] = (P[1] - t1) % (*p); t1 = t1 - P[1]; t1 = t1 * t0 % (*p); t1 = (t1 - P[0]) % (*p);
	P[0] = P[1]; P[1] = t1;
	return;
}

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
void New_2bit_2ary_RL(mpz_class P[2], string scalar, mpz_class A[2], mpz_class *a, mpz_class *p){
	//Initialization
	mpz_class *R[2];
	R[0] = new mpz_class[2];
	R[1] = new mpz_class[2];
	R[0][0] = P[0]; R[0][1] = -P[1];
	R[1][0] = P[0]; R[1][1] = P[1];
	mpz_class A1[2];
	Affine_DQ(P, A, A1, a, p);
	Affine_add(A, R[scalar[scalar.size()-1]-'0'], p);
	//cout<<A[0]<<' '<<A[1]<<endl;
	//cout<<A1[0]<<' '<<A1[1]<<endl;
	//cout<<R[0][0]<<' '<<R[0][1]<<endl;
	//cout<<R[1][0]<<' '<<R[1][1]<<endl;
	//Main Loop
	for (int i = scalar.size()-2; i >= 1; i -= 2){
		Affine_add(A, R[scalar[i]-'0'], p);
		Affine_add(A1, R[scalar[i-1]-'0'], p);
		Affine_DQ(A1, A, A1, a, p);
	}
	//Final Correction
	P[1] = -P[1];
	Affine_add(P, R[scalar[scalar.size()-1]-'0'], p);
	A[1] = -A[1];
	Affine_add(R[0], A, p);
	Affine_dbl_change(R[1], a, p);
	Extended_Affine(R[1], A, p);
	A[1] = pow(-1, int(scalar[0])) * A[1];
	P[1] = -P[1];
	return;
}
#pragma GCC pop_options
//CONTINUE

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
void Joye_LR(mpz_class P[3], string scalar, mpz_class A[3], mpz_class *a, mpz_class *b3, mpz_class *p){
	//Initialization
	mpz_class *R[2];
	R[1] = new mpz_class[3];
	Complete_addition(p, a, b3, &P[0], &P[1], &P[2], &P[0], &P[1], &P[2], R[1]);
	R[0] = P;
	A[0] = (scalar[1]-'0'-1) * P[0];
	A[1] = P[1];
	A[2] = (scalar[1]-'0'-1) * P[2];
	//Main Loop
	for (int i = 2; i < scalar.size(); ++i){
		Complete_addition(p, a, b3, &A[0], &A[1], &A[2], &A[0], &A[1], &A[2], A);
		Complete_addition(p, a, b3, &A[0], &A[1], &A[2], &R[scalar[i]-'0'][0], &R[scalar[i]-'0'][1], &R[scalar[i]-'0'][2], A);
	}
	//Final Correction
	Complete_addition(p, a, b3, &A[0], &A[1], &A[2], &R[0][0], &R[0][1], &R[0][2], A);
	A[1] = pow(-1, int(scalar[0])) * A[1];
	return; 
}
#pragma GCC pop_options

static inline void Affine_ODPQ(mpz_class P[2], mpz_class Q[2], mpz_class *a, mpz_class *p){
	mpz_class t0, t1, t2, t3;
	t0 = 2 * P[1]; t0 = t0 * t0 % (*p); Q[0] = Q[0] + P[0]; Q[0] = Q[0] + P[0];
	t1 = t0 * Q[0] % (*p); t1 = -t1; Q[0] = Q[0] - P[0]; Q[0] = Q[0] - P[0];
	t2 = P[0] * P[0] % (*p); t2 = 3 * t2; t2 = t2 + (*a); t3 = t2 * t2 % (*p);
	t1 = t1 + t3; t3 = t1 * P[1] % (*p); t3 = 2 * t3;
	mpz_invert(t3.get_mpz_t(), t3.get_mpz_t(), (*p).get_mpz_t());
	t2 = t2 * t3 % (*p); t2 = t2 * t1 % (*p); t3 = t3 * P[1] % (*p); t3 = t3 * t0 % (*p);
	t3 = 2 * t3; t1 = t2 * t2 % (*p); t1 = t1 - P[0]; t1 = t1 - P[0]; t0 = P[0] - t1; t0 = t0 * t2 % (*p);
	t0 = t0 - P[1]; t0 = t0 - Q[1]; t0 = t0 * t3 % (*p); P[0] = t0 * t0 % (*p); P[0] = P[0] - Q[0];
	P[0] = (P[0] - t1) % (*p); P[1] = Q[0] - P[0]; P[1] = P[1] * t0 % (*p); P[1] = (P[1] - Q[1]) % (*p); 
	return;
}

static inline void Affine_DPQ(mpz_class P[2], mpz_class Q[2], mpz_class *p){
	mpz_class t0, t1, t2;
	Q[1] = Q[1] - P[1]; t0 = Q[1] * Q[1] % (*p); t1 = 2 * P[0]; t1 = t1 + Q[0];
	Q[0] = Q[0] - P[0]; t2 = Q[0] * Q[0] % (*p); t1 = t1 * t2 % (*p);
	t0 = t0 - t1; t1 = t0 * Q[0] % (*p); mpz_invert(t1.get_mpz_t(), t1.get_mpz_t(), (*p).get_mpz_t());
	t0 = t0 * t1 % (*p); t0 = t0 * Q[1] % (*p);
	t1 = (-2) * t1; t1 = t1 * Q[0] % (*p); t1 = t1 * t2 % (*p); t1 = t1 * P[1] % (*p); t1 = t1 - t0;
	t2 = t1 + t0; t0 = t1 - t0; t0 = t0 * t2 % (*p); t0 = t0 + Q[0]; t0 = t0 + P[0];
	Q[0] = Q[0] + P[0]; P[0] = P[0] - t0; P[0] = P[0] * t1 % (*p); P[0] = P[0] - P[1];
	Q[1] = Q[1] + P[1]; P[1] = P[0] % (*p); P[0] = t0 % (*p);
	return;
}

static inline void Affine_PQR(mpz_class P[2], mpz_class Q[2], mpz_class R[2], mpz_class *p){
	mpz_class t0, t1, t2, t3;
	Q[1] = Q[1] - P[1]; t0 = Q[1] * Q[1] % (*p); t1 = Q[0] - P[0]; t2 = t1 * t1 % (*p);
	t3 = P[0] + Q[0]; t3 = t3 + R[0]; t3 = t3 * t2 % (*p); t3 = t0 - t3; t0 = t1 * t3 % (*p);
	mpz_invert(t0.get_mpz_t(), t0.get_mpz_t(), (*p).get_mpz_t());
	Q[1] = Q[1] * t0 % (*p); Q[1] = Q[1] * t3 % (*p); t3 = Q[1] * Q[1] % (*p); t3 = t3 - P[0]; t3 = t3 - Q[0];
	Q[0] = P[0] - t3; Q[0] = Q[0] * Q[1] % (*p); Q[0] = Q[0] - P[1]; Q[0] = Q[0] - R[1];
	t0 = t0 * t1 % (*p); t0 = t0 * t2 % (*p); t0 = t0 * Q[0] % (*p);
	Q[0] = t0 * t0 % (*p); Q[0] = Q[0] - R[0]; Q[0] = (Q[0] - t3) % (*p); Q[1] = R[0] - Q[0]; Q[1] = Q[1] * t0 % (*p); 
	Q[1] = (Q[1] - R[1]) % (*p);
	return;
}

static inline void Affine_PQDQ(mpz_class P[2], mpz_class Q[2], mpz_class *a, mpz_class *p){
	mpz_class t0, t1, t2;
	t0 = Q[0] - P[0]; t0 = 2 * t0; t0 = t0 * Q[1] % (*p);
	mpz_invert(t0.get_mpz_t(), t0.get_mpz_t(), (*p).get_mpz_t());
	t1 = Q[1] - P[1]; t1 = 2 * t1; t1 = t1 * t0 % (*p); t1 = t1 * Q[1] % (*p);
	t2 = t1 * t1 % (*p); t2 = t2 - P[0]; t2 = t2 - Q[0];
	P[0] = Q[0] - P[0]; t0 = t0 * P[0] % (*p); P[0] = Q[0] - P[0]; P[0] = P[0] - t2; P[0] = P[0] * t1 % (*p); P[0] = P[0] - P[1];
	P[1] = P[0] % (*p); P[0] = t2 % (*p); t1 = Q[0] * Q[0] % (*p); t1 = 3 * t1; t1 = t1 + (*a); t1 = t1 * t0 % (*p);
	t0 = t1 * t1 % (*p); t0 = t0 - Q[0]; t0 = t0 - Q[0]; Q[0] = Q[0] - t0; Q[0] = Q[0] * t1 % (*p); Q[0] = Q[0] - Q[1];
	Q[1] = Q[0] % (*p); Q[0] = t0 % (*p);
	return; 
}

static inline void Affine_PQRT(mpz_class P[2], mpz_class Q[2], mpz_class R[2], mpz_class T[2], mpz_class *p){
	mpz_class t0, t1;
	P[0] = P[0] - Q[0]; R[0] = R[0] - T[0]; t0 = P[0] * R[0] % (*p);
	mpz_invert(t0.get_mpz_t(), t0.get_mpz_t(), (*p).get_mpz_t());
	P[1] = P[1] - Q[1]; t1 = t0 * R[0] % (*p); t1 = t1 * P[1] % (*p);
	t0 = t0 * P[0] % (*p); P[0] = P[0] + Q[0]; P[1] = t1 * t1 % (*p); P[1] = P[1] - P[0]; P[1] = P[1] - Q[0]; P[0] = P[1] % (*p);
	P[1] = Q[0] - P[0]; P[1] = P[1] * t1 % (*p); P[1] = (P[1] - Q[1]) % (*p); R[1] = R[1] - T[1]; t0 = t0 * R[1] % (*p);
	R[0] = R[0] + T[0]; t1 = t0 * t0 % (*p); R[0] = t1 - R[0]; R[0] = (R[0] - T[0]) % (*p); R[1] = T[0] - R[0]; R[1] = R[1] * t0 % (*p);
	R[1] = (R[1] - T[1]) % (*p);
	return;
}

static inline void Affine_RLML(mpz_class P[2], mpz_class Q[2], mpz_class R[2], mpz_class T[2], mpz_class *a, mpz_class *p){
	mpz_class t0, t1, t2, t3;
	t0 = T[0] * T[0] % (*p); t1 = T[1] * T[1] % (*p); t1 = 2 * t1; t2 = t1 * t1 % (*p); t1 = t1 + T[0];
	t1 = t1 * t1 % (*p); t1 = t1 - t0; t1 = t1 - t2; t1 = 3 * t1; t0 = 3 * t0; t0 = t0 + (*a);
	t3 = t0 * t0 % (*p); t1 = t1 - t3; t1 = t1 * t0 % (*p); t2 = 2 * t2; t1 = t1 - t2; Q[0] = Q[0] - P[0];
	R[0] = R[0] - T[0]; t3 = T[1] * Q[0] % (*p); t3 = t3 * R[0] % (*p); t3 = t3 * t1 % (*p); t3 = 2 * t3;
	mpz_invert(t3.get_mpz_t(), t3.get_mpz_t(), (*p).get_mpz_t());
	Q[1] = Q[1] - P[1]; Q[1] = Q[1] * T[1] % (*p); Q[1] = Q[1] * t3 % (*p); Q[1] = Q[1] * R[0] % (*p);
	Q[1] = Q[1] * t1 % (*p); Q[1] = 2 * Q[1]; t3 = t3 * Q[0] % (*p); R[1] = R[1] - T[1]; R[1] = R[1] * T[1] % (*p);
	R[1] = R[1] * t3 % (*p); R[1] = R[1] * t1 % (*p); R[1] = 2 * R[1]; t3 = t3 * R[0] % (*p); t0 = t0 * t3 % (*p);
	t0 = t0 * t1 % (*p); t2 = t2 * t3 % (*p); t1 = Q[1] * Q[1] % (*p); t1 = t1 - Q[0]; t1 = t1 - P[0]; t1 = t1 - P[0];
	P[0] = P[0] - t1; P[0] = P[0] * Q[1] % (*p); P[0] = P[0] - P[1]; P[1] = P[0] % (*p); P[0] = t1 % (*p); 
	t1 = R[1] * R[1] % (*p); t1 = t1 - R[0]; t1 = t1 - T[0]; t1 = t1 - T[0]; R[0] = T[0] - t1; R[0] = R[0] * R[1] % (*p);
	R[0] = R[0] - T[1]; R[1] = R[0] % (*p); R[0] = t1 % (*p); Q[0] = t0 * t0 % (*p); Q[0] = Q[0] - T[0]; Q[0] = (Q[0] - T[0]) % (*p);;
	Q[1] = T[0] - Q[0]; Q[1] = Q[1] * t0 % (*p); Q[1] = (Q[1] - T[1]) % (*p); t1 = Q[0] * Q[0] % (*p); t1 = 3 *  t1;
	t1 = t1 + (*a); t1 = t1 * t2 % (*p); T[0] = t1 * t1 % (*p); T[0] = T[0] - Q[0]; T[0] = (T[0] - Q[0]) % (*p);
	T[1] = Q[0] - T[0]; T[1] = T[1] * t1 % (*p); T[1] = (T[1] - Q[1]) % (*p);
	return;
}

static inline void Affine_QPQ(mpz_class P[2], mpz_class Q[2], mpz_class *a, mpz_class *p){
	mpz_class t0, t1, t2, t3, t4, t5;
	t0 = P[0] * P[0] % (*p); t1 = P[1] * P[1] % (*p); t1 = 2 * t1; t2 = t1 * t1 % (*p); t3 = t1 + P[0];
	t3 = t3 * t3 % (*p); t3 = t3 - t0; t3 = t3 - t2; t0 = 3 * t0; t0 = t0 + (*a); 
	t4 = t0 * t0 % (*p); 
	t4 = t4 - t3; 
	t4 = t4 - t3; t3 = t3 - t4; t3 = t0 * t3 % (*p); t3 = t3 - t2; t3 = t3 - t2; t1 = t1 * Q[0] % (*p); 
	t1 = 2 * t1; t1 = t1 + t4; t1 = t1 + t4; t5 = t3 * t3 % (*p); t1 = t1 * t5 % (*p); t1 = 4 * t1; t5 = 4 * t2; 
	t5 = (*a) * t5; 
	t4 = t4 * t4 % (*p); 
	t4 = 3 * t4; 
	t4 = t4 + t5; 
	t4 = t4 * t4 % (*p);
	t4 = t4 - t1; 
	t5 = t3 * t4 % (*p); t5 = t5 * P[1] % (*p); t5 = 2 * t5; mpz_invert(t5.get_mpz_t(), t5.get_mpz_t(), (*p).get_mpz_t());
	t1 = t2 * t4 % (*p); t1 = t1 * t5 % (*p); t1 = 2 * t1;
	t5 = t3 * t5 % (*p); t0 = t0 * t4 % (*p); t0 = t0 * t5 % (*p);
	t2 = t2 * t2 % (*p); t5 = t2 * t5 % (*p); t5 = t5 * P[1] % (*p); t5 = 32 * t5; 
	t2 = t0 * t0 % (*p); t2 = t2 - P[0]; t2 = t2 - P[0];
	P[0] = P[0] - t2; P[0] = P[0] * t0 % (*p); P[0] = P[0] - P[1];
	P[1] = t2 * t2 % (*p); P[1] = 3 * P[1]; 
	P[1] = P[1] + (*a); t1 = t1 * P[1] % (*p); P[1] = t1 * t1 % (*p); P[1] = P[1] - t2; P[1] = P[1] - t2; 
	t2 = t2 - P[1]; 
	t2 = t2 * t1 % (*p); t2 = t2 - P[0]; 
	P[0] = P[0] * P[0] % (*p); P[0] = 4 * P[0]; t5 = t5 * P[0] % (*p);
	t2 = t2 - Q[1]; 
	t2 = t2 * t5 % (*p); P[0] = t2 * t2 % (*p); P[0] = P[0] - Q[0]; P[0] = (P[0] - P[1]) % (*p); P[1] = Q[0] - P[0]; 
	P[1] = P[1] * t2 % (*p); P[1] = (P[1] - Q[1]) % (*p);
	return;
}

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
void New_2ary_LR(mpz_class P[2], string scalar, mpz_class A[2], mpz_class *a, mpz_class *p){
	//Initialization
	Affine_dbl(P, A, a, p);
	mpz_class *R[2];
	R[0] = new mpz_class[2];
	R[0][0] = A[0]; R[0][1] = -A[1];
	P[1] = -P[1];
	R[1] = P;
	//cout<<A[0]<<' '<<A[1]<<endl;
	//cout<<R[0][0]<<' '<<R[0][1]<<endl;
	//cout <<R[1][0]<<' '<<R[1][1]<<endl;
	//Main Loop
	for (int i = 1; i < scalar.size()-1; ++i){
		Affine_dbl_change(A, a, p);
		Affine_add(R[scalar[i]-'0'], A, p);
	}
	//Final Correction
	Affine_dbl_change(A, a, p);
	Affine_add(R[0], A, p);
	Extended_Affine(R[scalar[scalar.size()-1]-'0'], A, p);
	A[1] = pow(-1, int(scalar[0])) * A[1];
	P[1] = -P[1];
	return;
}
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
void New_2bit_2ary_LR(mpz_class P[2], string scalar, mpz_class A[2], mpz_class *a, mpz_class *p){
	//Initialization
	mpz_class *R[4];
	P[1] = -P[1];
	R[1] = P;
	R[0] = new mpz_class[2];
	R[2] = new mpz_class[2];
	Affine_DQ(R[1], R[0], R[2], a, p);
	R[3] = R[0];
	A[0] = R[0][0]; A[1] = -R[0][1]; 
	//Main Loop
	for (int i = 1; i < scalar.size()-1; i += 2){
		Affine_DQ_LR(A, a, p);
		Affine_add(R[scalar[i]-'0' + 2], A, p);
		Affine_add(R[scalar[i+1]-'0'], A, p);
	}
	//Final Correction
	Affine_dbl_change(A, a, p);
	Affine_add(R[0], A, p);
	Extended_Affine(R[scalar[scalar.size()-1]-'0'], A, p);
	A[1] = pow(-1, int(scalar[0])) * A[1];
	P[1] = -P[1];
	return;
}
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
void EXNew_2bit_2ary_LR(mpz_class P[2], string scalar, mpz_class A[2], mpz_class *a, mpz_class *p){
	//Initialization
	mpz_class *R[4];
	R[0] = new mpz_class[2];
	R[1] = new mpz_class[2];
	R[2] = new mpz_class[2];
	Affine_DQ(P, A, R[2], a, p);
	Affine_add_(R[2], P, R[1], p); R[1][1] = -R[1][1];
	Affine_add_(R[2], A, R[0], p); R[0][1] = -R[0][1];
	Affine_add(A, P, p); P[1] = -P[1]; R[3] = P;
	R[2][1] = -R[2][1];
	//A[0] += (*p);
	//A[1] += (*p);
	//cout << A[0] << ' ' << A[1] <<endl;
	//cout << R[0][0] << ' ' << R[0][1] <<endl;
	//cout << R[1][0] << ' ' << R[1][1] <<endl;
	//cout << R[2][0] << ' ' << R[2][1] <<endl;
	//cout << R[3][0] << ' ' << R[3][1] <<endl;
	//Main Loop
	for (int i = 1; i < scalar.size()-1; i += 2){
		//cout << R[2 * (scalar[i]-'0') + (scalar[i+1]-'0')][0] << ' ' << R[2 * (scalar[i]-'0') + (scalar[i+1]-'0')][1] <<endl;
		Affine_QPQ(A, R[2 * (scalar[i]-'0') + (scalar[i+1]-'0')], a, p);
	}
	//Final Correction
	R[2][1] = -R[2][1];
	Affine_add(R[2], R[0], p);
	Affine_add(R[2], R[1], p);
	Affine_dbl_change(A, a, p);
	Affine_add(R[0], A, p);
	Extended_Affine(R[scalar[scalar.size()-1]-'0'], A, p);
	A[1] = pow(-1, int(scalar[0])) * A[1];
	R[2][1] = -R[2][1];
	Affine_add(R[2], P, p);
	return;
}
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
void New_2ary_RLPQDQ(mpz_class P[2], string scalar, mpz_class A[2], mpz_class *a, mpz_class *p){
	//Initialization
	mpz_class *R[2];
	R[0] = new mpz_class[2];
	R[1] = new mpz_class[2];
	R[0][0] = P[0]; R[0][1] = -P[1];
	R[1][0] = P[0]; R[1][1] = P[1];
	Affine_dbl(P, A, a, p);
	Affine_add(A, R[scalar[scalar.size()-1]-'0'], p);
	//cout<<A[0]<<' '<<A[1]<<endl;
	//cout<<R[0][0]<<' '<<R[0][1]<<endl;
	//cout <<R[1][0]<<' '<<R[1][1]<<endl;
	//Main Loop
	for (int i = scalar.size()-2; i >= 1; --i){
		Affine_PQDQ(R[scalar[i]-'0'], A, a, p);
	}
	//Final Correction
	P[1] = -P[1];
	Affine_add(P, R[scalar[scalar.size()-1]-'0'], p);
	A[1] = -A[1];
	Affine_add(R[0], A, p);
	Affine_dbl_change(R[1], a, p);
	Extended_Affine(R[1], A, p);
	A[1] = pow(-1, int(scalar[0])) * A[1];
	P[1] = -P[1];
	return;
}
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
void New_2ary_LRODPQ(mpz_class P[2], string scalar, mpz_class A[2], mpz_class *a, mpz_class *p){
	//Initialization
	Affine_dbl(P, A, a, p);
	mpz_class *R[2];
	R[0] = new mpz_class[2];
	R[0][0] = A[0]; R[0][1] = -A[1];
	P[1] = -P[1];
	R[1] = P;
	//cout<<R[0][0]<<' '<<R[0][1]<<endl;
	//cout <<R[1][0]<<' '<<R[1][1]<<endl;
	//Main Loop
	for (int i = 1; i < scalar.size()-1; ++i){
		Affine_ODPQ(A, R[scalar[i]-'0'], a, p);
	}
	//Final Correction
	Affine_DPQ(A, R[0], p);
	Extended_Affine(R[scalar[scalar.size()-1]-'0'], A, p);
	A[1] = pow(-1, int(scalar[0])) * A[1];
	P[1] = -P[1];
	return;
}
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
void New_2ary_LRDPQ(mpz_class P[2], string scalar, mpz_class A[2], mpz_class *a, mpz_class *p){
	//Initialization
	Affine_dbl(P, A, a, p);
	mpz_class *R[2];
	R[0] = new mpz_class[2];
	R[0][0] = A[0]; R[0][1] = -A[1];
	P[1] = -P[1];
	R[1] = P;
	//cout<<R[0][0]<<' '<<R[0][1]<<endl;
	//cout <<R[1][0]<<' '<<R[1][1]<<endl;
	//Main Loop
	for (int i = 1; i < scalar.size()-1; ++i){
		Affine_DPQ(A, R[scalar[i]-'0'], p);
	}
	//Final Correction
	Affine_DPQ(A, R[0], p);
	Extended_Affine(R[scalar[scalar.size()-1]-'0'], A, p);
	A[1] = pow(-1, int(scalar[0])) * A[1];
	P[1] = -P[1];
	return;
}
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
void New_2bit_2ary_RLPQRT(mpz_class P[2], string scalar, mpz_class A0[2], mpz_class *a, mpz_class *p){
	//Initialization
	mpz_class *R[4];
	mpz_class A1[2];
	R[0] = new mpz_class[2];
	R[1] = new mpz_class[2];
	R[2] = new mpz_class[2];
	R[3] = new mpz_class[2];
	R[0][0] = P[0]; R[0][1] = -P[1];
	R[1][0] = P[0]; R[1][1] = P[1];
	Affine_DQ(P, A0, A1, a, p);
	Affine_add(A0, R[scalar[scalar.size()-1]-'0'], p);
	Affine_add(A0, R[scalar[scalar.size()-2]-'0'], p);
	R[scalar[scalar.size()-3]-'0'+2][0] = A1[0];
	R[scalar[scalar.size()-3]-'0'+2][1] = A1[1];
	R[!(scalar[scalar.size()-3]-'0')+2][0] = A0[0];
	R[!(scalar[scalar.size()-3]-'0')+2][1] = -A0[1];
	Affine_DQ(A1, A0, A1, a, p);
	//cout<<A[0]<<' '<<A[1]<<endl;
	//cout<<A1[0]<<' '<<A1[1]<<endl;
	//cout<<R[0][0]<<' '<<R[0][1]<<endl;
	//cout<<R[1][0]<<' '<<R[1][1]<<endl;
	//Main Loop
	for (int i = scalar.size()-4; i >= 1; i -= 2){
		Affine_PQRT(R[scalar[i]-'0'], A0, R[scalar[i-1]-'0'+2], A1, p);
		Affine_DQ(A1, A0, A1, a, p);
	}
	//Final Correction
	Affine_add(R[2], R[0], p);
	Affine_add(R[3], R[1], p);
	Affine_dbl(P, A1, a, p);
	Affine_add(A1, R[!(scalar[scalar.size()-3]-'0')], p);
	P[1] = -P[1];
	Affine_add(P, R[scalar[scalar.size()-1]-'0'], p);
	A0[1] = -A0[1];
	Affine_dbl_change(R[1], a, p);
	Affine_add(R[1], A0, p);
	Extended_Affine(R[0], A0, p);
	A0[1] = pow(-1, int(scalar[0])) * A0[1];
	P[1] = -P[1];
	return;
}
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
void New_2bit_2ary_RLML(mpz_class P[2], string scalar, mpz_class A0[2], mpz_class *a, mpz_class *p){
	//Initialization
	mpz_class *R[4];
	mpz_class A1[2];
	R[0] = new mpz_class[2];
	R[1] = new mpz_class[2];
	R[2] = new mpz_class[2];
	R[3] = new mpz_class[2];
	R[0][0] = P[0]; R[0][1] = -P[1];
	R[1][0] = P[0]; R[1][1] = P[1];
	Affine_DQ(P, A0, A1, a, p);
	Affine_add(A0, R[scalar[scalar.size()-1]-'0'], p);
	Affine_add(A0, R[scalar[scalar.size()-2]-'0'], p);
	R[scalar[scalar.size()-3]-'0'+2][0] = A1[0];
	R[scalar[scalar.size()-3]-'0'+2][1] = A1[1];
	R[!(scalar[scalar.size()-3]-'0')+2][0] = A0[0];
	R[!(scalar[scalar.size()-3]-'0')+2][1] = -A0[1];
	Affine_DQ(A1, A0, A1, a, p);
	//cout<<A[0]<<' '<<A[1]<<endl;
	//cout<<A1[0]<<' '<<A1[1]<<endl;
	//cout<<R[0][0]<<' '<<R[0][1]<<endl;
	//cout<<R[1][0]<<' '<<R[1][1]<<endl;
	//Main Loop
	for (int i = scalar.size()-4; i >= 1; i -= 2){
		Affine_RLML(R[scalar[i]-'0'], A0, R[scalar[i-1]-'0'+2], A1, a, p);
	}
	//Final Correction
	Affine_add(R[2], R[0], p);
	Affine_add(R[3], R[1], p);
	Affine_dbl(P, A1, a, p);
	Affine_add(A1, R[!(scalar[scalar.size()-3]-'0')], p);
	P[1] = -P[1];
	Affine_add(P, R[scalar[scalar.size()-1]-'0'], p);
	A0[1] = -A0[1];
	Affine_dbl_change(R[1], a, p);
	Affine_add(R[1], A0, p);
	Extended_Affine(R[0], A0, p);
	A0[1] = pow(-1, int(scalar[0])) * A0[1];
	P[1] = -P[1];
	return;
}
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
void New_2bit_2ary_LRPQR(mpz_class P[2], string scalar, mpz_class A[2], mpz_class *a, mpz_class *p){
	//Initialization
	mpz_class *R[4];
	P[1] = -P[1];
	R[1] = P;
	R[0] = new mpz_class[2];
	R[2] = new mpz_class[2];
	Affine_DQ(R[1], R[0], R[2], a, p);
	R[3] = R[0];
	A[0] = R[0][0]; A[1] = -R[0][1]; 
	//Main Loop
	for (int i = 1; i < scalar.size()-1; i += 2){
		Affine_DQ_LR(A, a, p);
		Affine_PQR(R[scalar[i]-'0' + 2], A, R[scalar[i+1]-'0'], p);
	}
	//Final Correction
	Affine_dbl_change(A, a, p);
	Affine_add(R[0], A, p);
	Extended_Affine(R[scalar[scalar.size()-1]-'0'], A, p);
	A[1] = pow(-1, int(scalar[0])) * A[1];
	P[1] = -P[1];
	return;
}
#pragma GCC pop_options

mpz_class* MDIC(mpz_class* Inv, int n, mpz_class p){
	int M_count = 0;
	for(int i = 0; i < n-1; i+=2){
	 	Inv[n] = Inv[i]*Inv[i+1]%p;
	 	M_count += 1;
	 	n += 1;
	 }
	cout<< "Memory: " << n << endl;
	mpz_invert(Inv[n-1].get_mpz_t(), Inv[n-1].get_mpz_t(), p.get_mpz_t());
	int c = 2;
	for(int i = n-1; i-c >= 0; --i){
	  	mpz_class temp = Inv[i-c];
	  	Inv[i-c] = Inv[i]*Inv[i-c+1]%p;
	  	Inv[i-c+1] = Inv[i]*temp%p;
	  	M_count += 2;
	  	c+=1;
	}
	cout << "Multiple times: " << M_count << endl;
	return Inv;
}