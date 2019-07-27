

#pragma once

#include <list>
#include <vector>
#include <NTL/ZZ.h>
#include <math/mpz_class.hh>
#include <math.h>
#include <algorithm>    // std::min
#include <float.h>
#include <limits.h>
using namespace std;
struct encnum
{
	int  exponent;
	mpz_class mantissa;
	//double plaintext;
};
class Paillier {
 public:
    Paillier(const std::vector<mpz_class> &pk, gmp_randstate_t state);
    std::vector<mpz_class> pubkey() const { return { n, g }; }
    encnum encrypt_f(double plaintext);
    std::vector  <std::vector<encnum>> encryptMatrix(const std::vector < std::vector<double> > &p);
    encnum encode( double plaintext) const;
    
    encnum sub_f(const encnum &a, const encnum &b) const;
    encnum add_f(const encnum &a, const encnum &b) const;
    std::vector  <std::vector<encnum>> addMatrix(const std::vector < std::vector<encnum> > &c0, const std::vector < std::vector<encnum> > &c1) ;
    std::vector  <std::vector<encnum>> multConstMatrixbyMatrix(const std::vector < std::vector<double> > &c, const std::vector < std::vector<encnum> > &mat) ;
    std::vector  <std::vector<encnum>> multMatrixbyConstMatrix(const std::vector < std::vector<encnum> > &mat , const std::vector < std::vector<double> > &c) ;
    
    encnum decrease_exponent_too(const encnum &y,const int new_exp) const;
    encnum constMult_f(const double &m, const encnum &c) const; 
	vector < vector<encnum> > multMatrixbyConst( double alpha,const vector < vector<encnum> > &mat) ;

	
    mpz_class encrypt(const mpz_class &plaintext);
    mpz_class add(const mpz_class &c0, const mpz_class &c1) const;
    mpz_class sub(const mpz_class &c0, const mpz_class &c1) const;
    mpz_class constMult(const mpz_class &m, const mpz_class &c) const;
    mpz_class constMult(long m, const mpz_class &c) const;
    mpz_class constMult(const mpz_class &c, long m) const { return constMult(m,c); };
    mpz_class scalarize(const mpz_class &c);
    void refresh(mpz_class &c);
    mpz_class random_encryption();

    mpz_class dot_product(const std::vector<mpz_class> &c, const std::vector<mpz_class> &v);
    mpz_class dot_product(const std::vector<mpz_class> &c, const std::vector<long> &v);
    void rand_gen(size_t niter = 100, size_t nmax = 1000);
    double logbase(double base, double x) const;

    void printMatrix(std::vector  <std::vector<double>> &Mat);
    
    const double precision = 0.00000001 ; //0.00000001 0.001
    bool acc_per = false;
    int max_exponent = 0;
	//int exponent;
	//int int_rep;
    int BASE = 10;
    double LOG2_BASE = log2(BASE);
    double FLOAT_MANTISSA_BITS = LDBL_MANT_DIG;
    encnum zero_enc;
   // double max_int = n.get_d()/3 -1;
    //double nsquare = n.get_d() * n.get_d() ;
 protected:
    /* Public key */
    const mpz_class n, g;


    
    /* Randomness state */
    gmp_randstate_t _randstate;

    /* Cached values */
    const uint nbits;
    const mpz_class n2;
    bool good_generator;
    
    /* Pre-computed randomness */
    std::list<mpz_class> rqueue;
};

class Paillier_priv : public Paillier {
 public:
    Paillier_priv(const std::vector<mpz_class> &sk, gmp_randstate_t state);
    std::vector<mpz_class> privkey() const { return { p, q, g, a }; }
    void find_crt_factors();
    
    // if a !=0, and if you are encrypting using the private key, use this function
    // 75% speedup
    mpz_class encrypt(const mpz_class &plaintext);
    

    // no speedup compared to the fast_encrypt
    mpz_class fast_encrypt_precompute(const mpz_class &plaintext);

    mpz_class decrypt(const mpz_class &ciphertext) const;
    static std::vector<mpz_class> keygen(gmp_randstate_t state, uint nbits = 1024, uint abits = 256);
    double decrypt_f(encnum x) const;
    std::vector< std::vector<double> > decryptMatrix(const std::vector < std::vector<encnum> > &c);

 protected:
    /* Private key, including g from public part; n=pq */
    const mpz_class p, q;
    const mpz_class a;      /* non-zero for fast mode */

    /* Cached values */
    const bool fast;
    const mpz_class p2, q2;
    mpz_class e_p2, e_q2;
    const mpz_class two_p, two_q;
    const mpz_class pinv, qinv;
    const mpz_class hp, hq;
};

class Paillier_priv_fast : public Paillier_priv {
public:
    Paillier_priv_fast(const std::vector<mpz_class> &sk, gmp_randstate_t state);
    void precompute_powers();
    mpz_class compute_g_star_power(const mpz_class &x);
    static std::vector<mpz_class> keygen(gmp_randstate_t state, uint nbits = 1024);
    
    mpz_class encrypt(const mpz_class &plaintext);
private:
    const mpz_class g_star_;
    const mpz_class phi_n;
    const mpz_class phi_n2;
    const uint phi_n2_bits;
    std::vector<mpz_class> g_star_powers_p_;
    std::vector<mpz_class> g_star_powers_q_;
};
