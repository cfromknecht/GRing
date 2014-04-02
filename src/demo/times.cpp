#include <cstdlib>
#include <sys/time.h>
#include <sstream>

#include <flint/flint.h>
#include <flint/fmpz_mod_poly.h>

#define NUM_TESTS 14

int main( int /*argc*/, char** /*argv*/ ) {
  srand( time( NULL ) );

  size_t n = 256;
  size_t k = 14;

  flint_rand_t state;
  flint_randinit( state );

  fmpz_t q_z;
  fmpz_init_set_ui( q_z, 1 << k );

  // Allocated on heap
  fmpz_mod_poly_t* heap1 = new fmpz_mod_poly_t[1];
  fmpz_mod_poly_t* heap2 = new fmpz_mod_poly_t[1];
  fmpz_mod_poly_t* heapr = new fmpz_mod_poly_t[1];
  fmpz_mod_poly_init( *heap1, q_z );
  fmpz_mod_poly_init( *heap2, q_z );
  fmpz_mod_poly_init( *heapr, q_z );

  fmpz_mod_poly_t* heapF = new fmpz_mod_poly_t[1];
  fmpz_mod_poly_init( *heapF, q_z );
  fmpz_mod_poly_set_coeff_ui( *heapF, 0, 1 );
  fmpz_mod_poly_set_coeff_ui( *heapF, n, 1 );

  // Allocated on stack
  fmpz_mod_poly_t stack1, stack2, stackr;
  fmpz_mod_poly_init( stack1, q_z );
  fmpz_mod_poly_init( stack2, q_z );
  fmpz_mod_poly_init( stackr, q_z );

  fmpz_mod_poly_t stackF;
  fmpz_mod_poly_init( stackF, q_z );
  fmpz_mod_poly_set_coeff_ui( stackF, 0, 1 );
  fmpz_mod_poly_set_coeff_ui( stackF, n, 1 );

  struct timeval tv1, tv2;
  long multS = 0;
  long multH = 0;
  long randS = 0;
  long randH = 0;
  long addH = 0;
  long addS = 0;
  long copyHH = 0;
  long copySS = 0;
  long copySH = 0;
  long copyHS = 0;

  for ( size_t i = 0; i < NUM_TESTS; i++ ) {
    // Random heap initialization
    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_randtest( *heap1, state, n );
    gettimeofday( &tv2, 0 );
    randH += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );
    for ( size_t j = 0; j < n; j++ )
      fmpz_mod_poly_set_coeff_ui( *heap2, j, rand() );
    
    // Random stack initialization
    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_randtest( stack1, state, n );
    gettimeofday( &tv2, 0 );
    randS += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );
    for ( size_t j = 0; j < n; j++ )
      fmpz_mod_poly_set_coeff_ui( stack2, j, rand() );

    // Addition
    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_add( *heapr, *heap1, *heap2 );
    gettimeofday( &tv2, 0 );
    addH += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );

    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_add( stackr, stack1, stack2 );
    gettimeofday( &tv2, 0 );
    addS += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );

    // Multiplication 
    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_mulmod( *heapr, *heap1, *heap2, *heapF );
    gettimeofday( &tv2, 0 );
    multH += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );

    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_mulmod( stackr, stack1, stack2, stackF );
    gettimeofday( &tv2, 0 );
    multS += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );

    // Copying
    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_set( *heap1, *heap2 );
    gettimeofday( &tv2, 0 );
    copyHH += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );

    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_set( stack1, stack2 );
    gettimeofday( &tv2, 0 );
    copySS += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );

    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_set( *heapr, stackr );
    gettimeofday( &tv2, 0 );
    copySH += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );

    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_set( stack1, stackr );
    gettimeofday( &tv2, 0 );
    copySH += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );

    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_set( stack1, *heap2 );
    gettimeofday( &tv2, 0 );
    copyHS += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );
  }

  printf( "Random Polynomials\n" );
  printf( "\trand init stack time: %f us\n", randS / (float) NUM_TESTS );
  printf( "\trand init heap time: %f us\n", randH / (float) NUM_TESTS );
  printf( "\tmult stack time: %f us\n", multS / (float) NUM_TESTS );
  printf( "\tmult heap time: %f us\n", multH / (float) NUM_TESTS );
  printf( "\tadd stack time: %f us\n", addS / (float) NUM_TESTS );
  printf( "\tadd heap time: %f us\n", addH / (float) NUM_TESTS );
  printf( "\tcopy hh time: %f us\n", copyHH / (float) NUM_TESTS );
  printf( "\tcopy ss time: %f us\n", copySS / (float) NUM_TESTS );
  printf( "\tcopy hs time: %f us\n", copyHS / (float) NUM_TESTS );
  printf( "\tcopy sh time: %f us\n", copySH / (float) NUM_TESTS );

  randS = 0;
  randH = 0;
  multS = 0;
  multH = 0;
  addS = 0;
  addH = 0;

  // Ternary tests
  for ( size_t i = 0; i < NUM_TESTS; i++ ) {
    // Random heap initialization
    gettimeofday( &tv1, 0 );
    for ( size_t j = 0; j < n; j++ )
      fmpz_mod_poly_set_coeff_ui( *heap1, j, (rand() & 1) - (rand() & 1) );
    gettimeofday( &tv2, 0 );
    randH += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );
    for ( size_t j = 0; j < n; j++ )
      fmpz_mod_poly_set_coeff_ui( *heap2, j, (rand() & 1) - (rand() & 1) );
    
    // Random stack initialization
    gettimeofday( &tv1, 0 );
    for ( size_t j = 0; j < n; j++ )
      fmpz_mod_poly_set_coeff_ui( stack1, j, (rand() & 1) - (rand() & 1) );
    gettimeofday( &tv2, 0 );
    randS += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );
    for ( size_t j = 0; j < n; j++ )
      fmpz_mod_poly_set_coeff_ui( stack2, j, (rand() & 1) - (rand() & 1) );

    // Addition
    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_add( *heapr, *heap1, *heap2 );
    gettimeofday( &tv2, 0 );
    addH += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );

    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_add( stackr, stack1, stack2 );
    gettimeofday( &tv2, 0 );
    addS += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );

    // Multiplication 
    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_mulmod( *heapr, *heap1, *heap2, *heapF );
    gettimeofday( &tv2, 0 );
    multH += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );

    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_mulmod( stackr, stack1, stack2, stackF );
    gettimeofday( &tv2, 0 );
    multS += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );
  }

  printf( "Random Ternary Polynomials\n" );
  printf( "\trand init stack time: %f us\n", randS / (float) NUM_TESTS );
  printf( "\trand init heap time: %f us\n", randH / (float) NUM_TESTS );
  printf( "\tmult stack time: %f us\n", multS / (float) NUM_TESTS );
  printf( "\tmult heap time: %f us\n", multH / (float) NUM_TESTS );
  printf( "\tadd stack time: %f us\n", addS / (float) NUM_TESTS );
  printf( "\tadd heap time: %f us\n", addH / (float) NUM_TESTS );


  multS = 0;
  multH = 0;
  addS = 0;
  addH = 0;

  // Ternary and Random tests
  for ( size_t i = 0; i < NUM_TESTS; i++ ) {
    // Random heap initialization
    for ( size_t j = 0; j < n; j++ ) {
      fmpz_mod_poly_set_coeff_ui( *heap1, j, rand() );
      fmpz_mod_poly_set_coeff_ui( *heap2, j, (rand() & 1) - (rand() & 1) );
    }
    
    // Random stack initialization
    gettimeofday( &tv1, 0 );
    for ( size_t j = 0; j < n; j++ )
      fmpz_mod_poly_set_coeff_ui( stack1, j, (rand() & 1) - (rand() & 1) );
    gettimeofday( &tv2, 0 );
    randS += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );
    for ( size_t j = 0; j < n; j++ )
      fmpz_mod_poly_set_coeff_ui( stack2, j, (rand() & 1) - (rand() & 1) );

    // Addition
    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_add( *heapr, *heap1, *heap2 );
    gettimeofday( &tv2, 0 );
    addH += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );

    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_add( stackr, stack1, stack2 );
    gettimeofday( &tv2, 0 );
    addS += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );

    // Multiplication 
    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_mulmod( *heapr, *heap1, *heap2, *heapF );
    gettimeofday( &tv2, 0 );
    multH += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );

    gettimeofday( &tv1, 0 );
    fmpz_mod_poly_mulmod( stackr, stack1, stack2, stackF );
    gettimeofday( &tv2, 0 );
    multS += tv2.tv_usec - tv1.tv_usec + 1000000 * (tv2.tv_sec - tv1.tv_sec );
  }

  printf( "Random and Ternary Polynomials\n" );
  printf( "\tmult stack time: %f us\n", multS / (float) NUM_TESTS );
  printf( "\tmult heap time: %f us\n", multH / (float) NUM_TESTS );
  printf( "\tadd stack time: %f us\n", addS / (float) NUM_TESTS );
  printf( "\tadd heap time: %f us\n", addH / (float) NUM_TESTS );

  flint_randclear( state );
}

