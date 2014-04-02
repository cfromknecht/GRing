#include <ring/GRingArray.hpp>

#include <cstdlib>
#include <sys/time.h>
#include <ctime>
#include <iostream>

const size_t DIM = 256;         // main security parameter
const size_t K = 14;            // lg( Q )
const size_t L = 10;            // Length of ABE id
const size_t S = 4;             // Gaussian parameter

int main( void ) {
  srand( time( NULL ) );

  const size_t Q = 1 << K;

  fmpz_t q_z;
  fmpz_init_set_ui( q_z, Q );

  gring::GRingArray<DIM,Q,K,1> Abar{}, GAR{}, Ta{}, rPrime{}, e1{}, e2{};
  gring::ring_t y, s, e_y;
  fmpz_mod_poly_init2( y, q_z, DIM );
  fmpz_mod_poly_init2( s, q_z, DIM );
  fmpz_mod_poly_init2( e_y, q_z, DIM );

  for ( size_t i = 0; i < DIM; ++i ) {
    fmpz_mod_poly_set_coeff_ui( s, i, rand() );
    fmpz_mod_poly_set_coeff_ui( e_y, i, (rand()&1) - (rand()&1) );
  }

  // Setup
  Ta.ternaryInit( 0, 1 );
  Abar.makeParityCheckForSecret( GAR, 0, 1, Ta );

  // Recipient
  gring::invertId<DIM,Q,K,1>( y, rPrime, "conner" );

  // Encrypt
  e1.ternaryInit( 0, 1 );
  e2.ternaryInit( 0, 1 );

  Abar *= s;
  GAR  *= s;  
  Abar += e1;
  Abar += e2;

  fmpz_mod_poly_mulmod( y, y, s, Abar.FF() );
  fmpz_mod_poly_add( y, y, e_y );

  // Decrypt
  Abar *= Ta;
  Abar *= rPrime;
  GAR  *= rPrime;
  Abar += GAR;

  for ( size_t i = 1; i < K; ++i )
    fmpz_mod_poly_add( Abar.data()[0], Abar.data()[0], Abar.data()[i] );
  fmpz_mod_poly_sub( y, y, Abar.data()[0] );

  std::cout << "msg': " << std::endl;
  fmpz_mod_poly_print( y );
  std::cout << std::endl;

  return 0;
}

