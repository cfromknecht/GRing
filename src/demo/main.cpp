#include <ring/GRingArray.hpp>

#include <cstdlib>
#include <sys/time.h>
#include <ctime>
#include <iostream>
#include <random>

const size_t DIM = 256;         // main security parameter
const size_t K = 14;            // lg( Q )
const size_t L = 6;            // Length of ABE id
const size_t S = 4;             // Gaussian parameter

std::default_random_engine generator;

long
sampleNormal( long Q, long mu ) {
  std::normal_distribution<float> distribution{(float)mu, 8.0};
  float sample = distribution( generator );
  return long(sample < 0 ? sample + Q : sample) % Q;
}

void 
printTimes( std::string action, const struct timeval t1, const struct 
    timeval t2 ) {
  long millis = t2.tv_usec - t1.tv_usec + 1000000.0*(t2.tv_sec - t1.tv_sec);
  std::cout << action << " time: " << millis << std::endl;
  std::cout << action << "s/sec: " << 1000000.0/millis << std::endl;
}

template< size_t N, size_t Q, size_t K, size_t LEN >
void 
printPolyAtIndex( std::string name, gring::GRingArray<N,Q,K,LEN>& array, size_t 
    index ) {
  std::cout << name << "[" << index << "]: " << std::endl;
  fmpz_mod_poly_print( array.data()[index] );
  std::cout << std::endl;
}

int main( void ) {
  srand( time( NULL ) );

  const size_t Q = 1 << K;
  const char* ID = "conner";

  struct timeval t1, t2;

  std::cout << "IBE" << std::endl;

  // Initialize Polynomial Ring Array (arrays of R^k_q_f(x) where f(x) = x^k + 1 
  // and k is a multiple of 2)
  gring::GRingArray<DIM,Q,K,2> A{}, Ta{}, skID{};

  // Make PK and SK
  gettimeofday( &t1, 0 );

  A.makeParityCheckForSecret( 0, Ta);

  gettimeofday( &t2, 0 );
  printTimes( "setup", t1, t2 );
  
  // Make Public Key for ID
  gettimeofday( &t1, 0 );

  fmpz_t q_z;
  fmpz_init_set_ui( q_z, Q );

  gring::ring_t y;
  fmpz_mod_poly_init2( y, q_z, DIM );
  skID.invertId( y, ID, 1 );
  Ta.mulmodGRingSection( 0, 1, skID.data(), 1, skID.data(), 0 );
//  skID *= Ta;

  gettimeofday( &t2, 0 );
  printTimes( "inversion", t1, t2 );

  // Encrypt
  gettimeofday( &t1, 0 );

  gring::ring_t s;
  fmpz_mod_poly_init2( s, q_z, DIM );

  A.makeUniformPoly( s );
  A.dualRegevEncode( s );
  fmpz_mod_poly_mulmod( y, y, s, A.FF() );
  A.ternaryPerturb( y );
  
  fmpz_mod_poly_clear( s );

  gettimeofday( &t2, 0 );
  printTimes( "encryption", t1, t2 );

  // Decrypt
  gettimeofday( &t1, 0 );

  A *= skID;
  for ( size_t i = 1; i < 2*K; ++i )
    fmpz_mod_poly_add( A.data()[0], A.data()[0], A.data()[i] );
  fmpz_mod_poly_sub( y, y, A.data()[0] );

  gettimeofday( &t2, 0 );
  printTimes( "decryption", t1, t2 );

  std::cout << "msg': " << std::endl;
  fmpz_mod_poly_print( y );
  std::cout << std::endl;


  std::cout << "ABE" << std::endl;

  gring::GRingArray<DIM,Q,K,L+2> abeA{}, abeTa{};
  gring::GRingArray<DIM,Q,K,L+1> H{}, e{};
  gring::GRingArray<DIM,Q,K,1> Bf{}, Sf{}, cg{}, D{};

  // Make PK and SK
  gettimeofday( &t1, 0 );

  abeA.makeParityCheckForSecret( 0, abeTa);
  abeA.uniformInit( 1, L + 2 );
  for ( size_t i = K; i < K*(L+2); ++i )
    fmpz_mod_poly_set( abeTa.data()[i], abeA.data()[i] );

  gettimeofday( &t2, 0 );
  printTimes( "setup", t1, t2 );

  std::cout << "making ciphertext ..." << std::endl;

  for ( size_t i = 0; i < K; ++i )
    fmpz_mod_poly_set( H.data()[i], abeA.data()[i] );
  size_t index = 0;
  for ( size_t i = K*(L+1); i < K*(L+2); ++i )
    fmpz_mod_poly_set( D.data()[index++], abeA.data()[i] );


  gring::ring_t coeffPoly;
  fmpz_mod_poly_init2( coeffPoly, q_z, 1 );
#pragma omp parallel for
  for ( size_t i = 1; i < L+1; ++i ) 
    for ( size_t j = 0; j < K; ++j ) {
      fmpz_mod_poly_set_coeff_ui( coeffPoly, 0, ID[i-1]*(1 << j ) );
      fmpz_mod_poly_add( H.data()[K*i + j], abeA.data()[K*i + j], coeffPoly );
    }

  std::cout << "making error vector ..." << std::endl;

  gring::ring_t e0, e1;
  fmpz_mod_poly_init2( e0, q_z, DIM );
  fmpz_mod_poly_init2( e1, q_z, DIM );

  e.ternaryPerturb( e0 );
  e.ternaryPerturb( e1 );

  e.identityInit( 0, 1 );
  e.ternaryInit( 1, L+1 );

  std::cout << "encoding ..." << std::endl;

  H.makeUniformPoly( s );
  H *= s;
  D *= s;

  H += e * e0;
  D += e1;

  std::cout << "eval pk, eval sim, eval ct ..." << std::endl;
  for ( size_t i = 1; i < L+1; ++i ) {
    Bf.plusGRingSection(0, 1, abeA.data(), i, Bf.data(), 0 );
    Sf.plusGRingSection(0, 1, e.data(), i, Sf.data(), 0 );
    cg.plusGRingSection(0, 1, H.data(), i, cg.data(), 0 );
  }

  std::cout << "Gaussian sample:" << std::endl;
  for ( size_t i = 0; i < 100; ++i )
    std::cout << sampleNormal( Q, 0 ) << " ";
  std::cout << std::endl;

  fmpz_clear( q_z );

  return 0;
}

