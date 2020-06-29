#include <fstream>
#include <time.h>
#include <windows.h>
#include "mpz_add.h"

void read_data(string name, mpz_class* K){
	ifstream f(name);
	string scalar;
	int index = 0;
	while(f >> scalar){
		K[index] = scalar;
		index += 1;
	}
	f.close();
	return;
}

void Test(float *R_clock, int bit, mpz_class a, mpz_class b, mpz_class b3, mpz_class p, mpz_class x, mpz_class y, mpz_class order, mpz_class* testarray, int testnum){
	mpz_class mis_count = 0;
	double clock_[13] = {0.0};
	LARGE_INTEGER start;
	LARGE_INTEGER end;
	LARGE_INTEGER frequency;
	QueryPerformanceFrequency(&frequency); 
	for (int i = 0; i < testnum; ++i){
		string scalar = testarray[i].get_str(2);
		//scalar = "0" + scalar;
		if(scalar[0] == '-'){
			scalar = scalar.substr(1);
			while(scalar.length() < bit-1){
				scalar = "0" + scalar;
			}
			scalar = "1" + scalar;
		}
		else{
			while(scalar.length() < bit-1){
				scalar = "0" + scalar;
			}
			scalar = "0" + scalar;
		}

		mpz_class P[3] = {x, y, mpz_class("1")};
		mpz_class P_[2] = {x, y};
		mpz_class PR[2][3];
		mpz_class PR_[11][2];	
		//mpz_class P[3] = {x, y, mpz_class("1")};
		//mpz_class k_P[3];
		QueryPerformanceCounter(&start);
		Joye_RL(P, scalar, PR[0], &a, &b3, &p);
		QueryPerformanceCounter(&end);
		clock_[0] += (double)(end.QuadPart-start.QuadPart)/frequency.QuadPart;	
		generalization(p, PR[0]);
		//clock1 += (clock() - start) * 1.0 / CLOCKS_PER_SEC * 1000;	
		//cout << "Joye's RL + CA: " << PR[0] << ' ' << PR[1] << endl;
		//cout<<endl;
		//cout <<"Alg4 + Complete_addition "<< k_P[0] << ' ' << k_P[1]<< endl;
		
		QueryPerformanceCounter(&start);
		New_2ary_RL(P_, scalar, PR_[0], &a, &p);
		QueryPerformanceCounter(&end);
		clock_[1] += (double)(end.QuadPart-start.QuadPart)/frequency.QuadPart;
		generalization_aff(PR_[0], p);
		//cout << "RL + (extended) affine: " << PR_[0] << ' ' << PR_[1] << endl;
		//cout<<endl;
		
		//clock2 += (clock() - start) * 1.0 / CLOCKS_PER_SEC * 1000;
		//generalization_aff(k_P_, p);
	    //cout <<"New_2ary_RL "<< k_P_[0] << ' ' << k_P_[1]<< endl;
		
		QueryPerformanceCounter(&start);
		New_2ary_RLPQDQ(P_, scalar, PR_[1], &a, &p);
		QueryPerformanceCounter(&end);
		clock_[2] += (double)(end.QuadPart-start.QuadPart)/frequency.QuadPart;
		generalization_aff(PR_[1], p);
		//cout << "RL + (extended) affine + PQ2Q: " << PR_[0] << ' ' << PR_[1] << endl;
		//cout<<endl;
		string temp;
		if(scalar.size()%2 != 0){
			temp = scalar.substr(1);
			temp = scalar.substr(0, 1) + '0' + temp;
		}
		else{
			temp = scalar;
		}
		QueryPerformanceCounter(&start);
	    New_2bit_2ary_RL(P_, temp, PR_[2], &a, &p);
	    QueryPerformanceCounter(&end);
	    clock_[3] += (double)(end.QuadPart-start.QuadPart)/frequency.QuadPart;
		generalization_aff(PR_[2], p);
		//cout << "RL2bit + (extended) affine: " << PR_[0] << ' ' << PR_[1] << endl;
		//cout<<endl;
	    
	    //clock3 += (clock() - start) * 1.0 / CLOCKS_PER_SEC * 1000;
	    //generalization_aff(k_P_2, p);
	    //cout <<"New_2bit_2ary_RL "<< k_P_2[0] << ' ' << k_P_2[1]<< endl;
	   	QueryPerformanceCounter(&start);
	   	New_2bit_2ary_RLPQRT(P_, temp, PR_[3], &a, &p);
	   	QueryPerformanceCounter(&end);
		clock_[4] += (double)(end.QuadPart-start.QuadPart)/frequency.QuadPart;
		generalization_aff(PR_[3], p);
		//cout << "RL2bit + (extended) affine + PQRT: " << PR_[0] << ' ' << PR_[1] << endl;
		//cout<<endl;

		QueryPerformanceCounter(&start);
		New_2bit_2ary_RLML(P_, temp, PR_[4], &a, &p);
		QueryPerformanceCounter(&end);
		clock_[5] += (double)(end.QuadPart-start.QuadPart)/frequency.QuadPart;
		generalization_aff(PR_[4], p);
		//cout << "RL2bit + (extended) affine + RL2bitML: " << PR_[0] << ' ' << PR_[1] << endl;
		//cout<<endl;

		QueryPerformanceCounter(&start);
		Joye_LR(P, scalar, PR[1], &a, &b3, &p);
		QueryPerformanceCounter(&end);	
		clock_[6] += (double)(end.QuadPart-start.QuadPart)/frequency.QuadPart;
		generalization(p, PR[1]);
		//cout << "Joye's LR + CA: " << PR[0] << ' ' << PR[1] << endl;
		//cout<<endl;

		QueryPerformanceCounter(&start);
		New_2ary_LR(P_, scalar, PR_[5], &a, &p);
		QueryPerformanceCounter(&end);
		clock_[7] += (double)(end.QuadPart-start.QuadPart)/frequency.QuadPart;
		generalization_aff(PR_[5], p);
		//cout << "LR + (extended) affine: " << PR_[0] << ' ' << PR_[1] << endl;
		//cout<<endl;

		int zero_c = 0;
		string temp_ = scalar;
		for(int i = 1; i < scalar.size(); ++i){
			if(scalar[i]-'0' == 0){
				zero_c += 1;
			}
			else{
				break;
			}
		}
		while(zero_c%2 != 0){
			zero_c += 1;
			temp_ = temp_.substr(0, 1) + '0' + temp_.substr(1);
		}
		QueryPerformanceCounter(&start);
		New_2ary_LRDPQ(P_, temp_, PR_[6], &a, &p);
		QueryPerformanceCounter(&end);
		clock_[8] += (double)(end.QuadPart-start.QuadPart)/frequency.QuadPart;
		generalization_aff(PR_[6], p);
		//cout << "LR + (extended) affine + 2PQ: " << PR_[0] << ' ' << PR_[1] << endl;
		//cout<<endl;

		QueryPerformanceCounter(&start);
		New_2ary_LRODPQ(P_, scalar, PR_[7], &a, &p);
		QueryPerformanceCounter(&end);
		clock_[9] += (double)(end.QuadPart-start.QuadPart)/frequency.QuadPart;
		generalization_aff(PR_[7], p);
		//cout << "LR + (extended) affine + O2PQ: " << PR_[0] << ' ' << PR_[1] << endl;
		//cout<<endl;

		QueryPerformanceCounter(&start);
		New_2bit_2ary_LR(P_, temp, PR_[8], &a, &p);
		QueryPerformanceCounter(&end);
		clock_[10] += (double)(end.QuadPart-start.QuadPart)/frequency.QuadPart;
		generalization_aff(PR_[8], p);
		//cout << "LR2bit + (extended) affine: " << PR_[0] << ' ' << PR_[1] << endl;
		//cout<<endl;

		QueryPerformanceCounter(&start);
		New_2bit_2ary_LRPQR(P_, temp, PR_[9], &a, &p);
		QueryPerformanceCounter(&end);
		clock_[11] += (double)(end.QuadPart-start.QuadPart)/frequency.QuadPart;
		generalization_aff(PR_[9], p);
		//cout << "LR2bit + (extended) affine + PQR: " << PR_[0] << ' ' << PR_[1] << endl;
		//cout<<endl;

		QueryPerformanceCounter(&start);
		EXNew_2bit_2ary_LR(P_, temp, PR_[10], &a, &p);
		QueryPerformanceCounter(&end);
		clock_[12] += (double)(end.QuadPart-start.QuadPart)/frequency.QuadPart;
		generalization_aff(PR_[10], p);
		//cout << "EXLR2bit + (extended) affine + QPQ: " << PR_[0] << ' ' << PR_[1] << endl;
		//cout<<endl;
		if(PR[0][0] != PR [1][0] || PR[0][1] != PR [1][1]){
			cout<<"ERROR: RLCA or LRCA!";	
		}
		for (int i = 0; i < 11; ++i){
			if(PR_[i][0] != PR[0][0] || PR_[i][1] != PR[0][1]){
				cout<<"ERROR: "<<"method "<<i<<' '<<testarray[i]<<endl;
	    		mis_count += 1;
			}
		}
	   	
	}
	//R_clock = new float[13];
	for (int i = 0; i < 13; ++i){
		R_clock[i] = clock_[i] * 1000.0; // CLOCKS_PER_SEC *1000;
	}
	cout<<mis_count<<endl;
	cout <<"Joye's RL + CA"<<"(time ms): "<<R_clock[0]<<endl;
	cout <<"RL + (extended) affine"<<"(time ms): "<<R_clock[1]<<endl;
	cout <<"RL + (extended) affine + PQ2Q"<<"(time ms): "<<R_clock[2]<<endl;
	cout <<"RL2bit + (extended) affine"<<"(time ms): "<<R_clock[3]<<endl;
	cout <<"RL2bit + (extended) affine + PQRT"<<"(time ms): "<<R_clock[4]<<endl;
	cout <<"RL2bit + (extended) affine + RL2bitML"<<"(time ms): "<<R_clock[5]<<endl;
	cout <<"Joye's LR + CA"<<"(time ms): "<<R_clock[6]<<endl;
	cout <<"LR + (extended) affine"<<"(time ms): "<<R_clock[7]<<endl;
	cout <<"LR + (extended) affine + 2PQ"<<"(time ms): "<<R_clock[8]<<endl;
	cout <<"LR + (extended) affine + O2PQ"<<"(time ms): "<<R_clock[9]<<endl;
	cout <<"LR2bit + (extended) affine"<<"(time ms): "<<R_clock[10]<<endl;
	cout <<"LR2bit + (extended) affine + PQR"<<"(time ms): "<<R_clock[11]<<endl;
	cout <<"EXLR2bit + (extended) affine + QPQ"<<"(time ms): "<<R_clock[12]<<endl;
	return;
}

int main(){
	mpz_class P_224_a, P_224_b, P_224_3b, P_224_p, P_224_x, P_224_y, P_224_order;
	mpz_class P_256_a, P_256_b, P_256_3b, P_256_p, P_256_x, P_256_y, P_256_order;
	mpz_class P_384_a, P_384_b, P_384_3b, P_384_p, P_384_x, P_384_y, P_384_order;

	P_224_a = "-3";
	P_224_b = "18958286285566608000408668544493926415504680968679321075787234672564";
	P_224_3b = 3 * P_224_b;
	P_224_p = "26959946667150639794667015087019630673557916260026308143510066298881";
	P_224_x = "19277929113566293071110308034699488026831934219452440156649784352033";
	P_224_y = "19926808758034470970197974370888749184205991990603949537637343198772";
	P_224_order = "26959946667150639794667015087019625940457807714424391721682722368061";

	P_256_a = "-3";
	P_256_b = "41058363725152142129326129780047268409114441015993725554835256314039467401291";
	P_256_3b = 3 * P_256_b;
	P_256_p = "115792089210356248762697446949407573530086143415290314195533631308867097853951";
	P_256_x = "48439561293906451759052585252797914202762949526041747995844080717082404635286";
	P_256_y = "36134250956749795798585127919587881956611106672985015071877198253568414405109";
	P_256_order = "115792089210356248762697446949407573529996955224135760342422259061068512044369";

	P_384_a = "-3";
	P_384_b = "27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575";
	P_384_3b = 3 * P_384_b;
	P_384_p = "39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319";
	P_384_x = "26247035095799689268623156744566981891852923491109213387815615900925518854738050089022388053975719786650872476732087";
	P_384_y = "8325710961489029985546751289520108179287853048861315594709205902480503199884419224438643760392947333078086511627871";
	P_384_order = "39402006196394479212279040100143613805079739270465446667946905279627659399113263569398956308152294913554433653942643";
	
	float R224[5][13] = {0};
	float R256[5][13] = {0};
	float R384[5][13] = {0};

	for(int i = 1; i < 5; ++i){
		int data_number = 10000 * i;
		mpz_class P224_scalar[data_number];
		mpz_class P256_scalar[data_number];
		mpz_class P384_scalar[data_number];

		read_data("P_224_"+ to_string(data_number) +"+-.txt", P224_scalar);
		read_data("P_256_"+ to_string(data_number) +"+-.txt", P256_scalar);
		read_data("P_384_"+ to_string(data_number) +"+-.txt", P384_scalar);
	/*
	for(int i = 0; i < 10000; ++i){
		cout<<P224_scalar[i]<<endl;
	}
	*/
		Test(R224[i-1], 224,P_224_a, P_224_b, P_224_3b, P_224_p, P_224_x, P_224_y, P_224_order, P224_scalar, data_number);
		Test(R256[i-1], 256,P_256_a, P_256_b, P_256_3b, P_256_p, P_256_x, P_256_y, P_256_order, P256_scalar, data_number);
		Test(R384[i-1], 384,P_384_a, P_384_b, P_384_3b, P_384_p, P_384_x, P_384_y, P_384_order, P384_scalar, data_number);
	/*
	mpz_class count = 0;
	for (int i = 0; i < 10000; ++i){
		string scalar = P224_scalar[i].get_str(2);
		scalar = "0" + scalar;

		mpz_class P[3] = {P_224_x, P_224_y, mpz_class("1")};
		mpz_class k_P[3];
		Joye_RL(P, scalar, k_P, &P_224_a, &P_224_3b, &P_224_p);	
		generalization(P_224_p, k_P);
		//cout <<"Alg4 + Complete_addition "<< k_P[0] << ' ' << k_P[1]<< endl;
		
		mpz_class P_[2] = {P_224_x, P_224_y};
		mpz_class k_P_[2];
		New_2ary_RL(P_, scalar, k_P_, &P_224_a, &P_224_p);
		generalization_aff(k_P_, P_224_p);
	    //cout <<"New_2ary_RL "<< k_P_[0] << ' ' << k_P_[1]<< endl;
		
		if(scalar.size()%2 == 0) scalar = '0'+ scalar;
		mpz_class P_2[2] = {P_224_x, P_224_y};
		mpz_class k_P_2[2];
	    New_2bit_2ary_RL(P_2, scalar, k_P_2, &P_224_a, &P_224_p);
	    generalization_aff(k_P_2, P_224_p);
	    //cout <<"New_2bit_2ary_RL "<< k_P_2[0] << ' ' << k_P_2[1]<< endl;
	   	
	    if (k_P[0] != k_P_[0] || k_P[0] != k_P_2[0] || k_P_[0] != k_P_2[0] || k_P[1] != k_P_[1] || k_P[1] != k_P_2[1] || k_P_[1] != k_P_2[1]){
	    	//cout << k_P[0] << ' ' << k_P[1]<< endl;
	    	cout<<"ERROR: "<<P224_scalar[i]<<endl;
	    	count += 1;
	    }
	}
	cout<<count;
    */
	}
	cout<<"begin!"<<endl;
	for(int i = 0; i < 13; ++i){
		for(int j = 0; j < 4; ++j){
			R224[4][i] += R224[j][i];
			R256[4][i] += R256[j][i];
			R384[4][i] += R384[j][i];
		}
		R224[4][i] /= 100000;
		R256[4][i] /= 100000;
		R384[4][i] /= 100000;
	}

	for(int i = 0; i < 5; ++i){
		for(int j = 0; j < 13; ++j){
			cout<<R224[i][j]<<' ';
		}
		cout<<endl;
	}
	for(int i = 0; i < 5; ++i){
		for(int j = 0; j < 13; ++j){
			cout<<R256[i][j]<<' ';
		}
		cout<<endl;
	}
	for(int i = 0; i < 5; ++i){
		for(int j = 0; j < 13; ++j){
			cout<<R384[i][j]<<' ';
		}
		cout<<endl;
	}
    char pause;
    pause = getchar();
    cout<<pause;
    return 0;
}