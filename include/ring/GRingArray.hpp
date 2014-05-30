#ifndef _G_RING_ARRAY_
#define _G_RING_ARRAY_

#include <flint/fmpz_mod_poly.h>
#include <openssl/sha.h>
#include <omp.h>
#include <sys/time.h>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <cstring>

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

namespace gring {

  typedef fmpz_mod_poly_t ring_t;

  template< size_t N, size_t Q, size_t K >
  class GRingArray {
  private:
    ring_t* _polys;
    ring_t& _FF;
    size_t _len;
    ring_t& makeF() const; 
    GRingArray() = delete;
  public:
    GRingArray( size_t len );
    GRingArray( const ring_t* polys, size_t length ) : _polys{polys}, 
      _FF{makeF()}, _len{length} {}
    GRingArray( const GRingArray& other );
    GRingArray& operator=( const GRingArray& other );
    ~GRingArray();
    // Accessors
    inline size_t NN() const { return N; }
    inline size_t QQ() const { return Q; }
    inline size_t KK() const { return K; }    
    inline size_t len() const { return _len; }
    inline size_t size() const { return K*_len; }
    inline ring_t* data() const { return _polys; }
    inline ring_t& FF() const { return _FF; }
    // Initializers
    inline void uniformInit() { uniformInit( 0, _len ); }
    void uniformInit( const size_t start, const size_t end );
    inline void ternaryInit() { ternaryInit( 0, _len ); }
    void ternaryInit( const size_t start, const size_t end );
    void identityInit( const size_t index );
    void identityInit( const size_t start, const size_t end );
    // Arithmetic Operations
    GRingArray operator-( const ring_t& s ) const;
    GRingArray operator-( const GRingArray& other ) const;
    void operator-=( const ring_t& s );
    void operator-=( const GRingArray& other );
    GRingArray operator+( const ring_t& s ) const;
    GRingArray operator+( const GRingArray& other ) const;
    void operator+=( const ring_t& s );
    void operator+=( const GRingArray& other );
    GRingArray operator*( const ring_t& s ) const;
    GRingArray operator*( const GRingArray& other ) const;
    void operator*=( const ring_t& s );
    void operator*=( const GRingArray& other );
    void negGRingSection ( const size_t start, const size_t end );
    void minusGRingSection( const size_t start, const size_t end, const ring_t&
        otherPoly, GRingArray& target, const size_t targetStart );
    void minusGRingSection( const size_t start, const size_t end, const 
        GRingArray& other, const size_t otherStart, GRingArray& target, const 
        size_t targetStart );
    void plusGRingSection( const size_t start, const size_t end, const ring_t&
        otherPoly, GRingArray& target, const size_t targetStart );
    void plusGRingSection( const size_t start, const size_t end, const 
        GRingArray& other, const size_t otherStart, GRingArray& target, const 
        size_t targetStart );
    void mulmodGRingSection( const size_t start, const size_t end, const ring_t& 
        otherPoly, GRingArray& target, const size_t targetStart );
    void mulmodGRingSection( const size_t start, const size_t end, const 
        GRingArray& other, const size_t otherStart, GRingArray& target, const 
        size_t targetStart );
    // Cryptographic Operations
    void dualRegevEncode();
    void dualRegevEncode( ring_t& s );
    ring_t* dualRegevEncrypt( std::string& data );
    void dualRegevEncrypt( ring_t* s, std::string& data );
    void ternaryPerturb( ring_t& r );
    // Factory Methods
    void makeUniformPoly( ring_t& u );
  };

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>::GRingArray( size_t length ) : _polys{new ring_t[K*length]}, 
      _FF{makeF()}, _len{length} {
    assert( (K % 2) == 0 );
    assert( (size_t(1) << K) == Q );

    fmpz_t q_z;
    fmpz_init_set_ui( q_z, Q );
#pragma omp parallel for
    for ( size_t i = 0; i < K*_len; ++i )
      fmpz_mod_poly_init2( _polys[i], q_z, N );
    fmpz_clear( q_z );
  }

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>::GRingArray( const GRingArray& other ) : 
      _polys{new ring_t[K*other.len()]}, _FF{makeF()}, _len{other.len()} {
    assert( (size_t(1) << K) >= Q );
    assert( Q > (size_t(1) << (K-1)) );

    ring_t* otherPolys = other.data();
    fmpz_t q_z;
    fmpz_init_set_ui( q_z, Q );
#pragma omp parallel for
    for ( size_t i = 0; i < K*_len; ++i ) {
      fmpz_mod_poly_init2( _polys[i], q_z, N );
      fmpz_mod_poly_set( _polys[i], otherPolys[i] );
    }
    fmpz_clear( q_z );
  }

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>::~GRingArray() {
#pragma omp parallel for
    for ( size_t i = 0; i < K*_len; ++i )
      fmpz_mod_poly_clear( _polys[i] );
    delete [] _polys;
    fmpz_mod_poly_clear( _FF );
  }

  template< size_t N, size_t Q, size_t K >
  void
  GRingArray<N,Q,K>::uniformInit( const size_t start, const size_t end ) {
    assert( end <= _len );
    assert( start < end );

#pragma omp parallel for 
    for ( size_t i = K*start; i < K*end; ++i )
      for ( size_t j = 0; j < N; ++j )
        fmpz_mod_poly_set_coeff_ui( _polys[i], j, rand() );
    return;
  }

  template< size_t N, size_t Q, size_t K >
  void
  GRingArray<N,Q,K>::ternaryInit( const size_t start, const size_t end ) {
    assert( end <= _len );
    assert( start < end );

    size_t randbuff = 0;
#pragma omp parallel for
    for ( size_t i = K*start; i < K*end; ++i )
      for ( size_t j = 0; j < N; ++j ) {
        if ( j % 32 == 0 ) randbuff = rand();
        fmpz_mod_poly_set_coeff_ui( _polys[i], j, (randbuff&1) + ((randbuff>>1)&1)*(Q-1) );
        randbuff >>= 2;
      }
    return;
  }

  template< size_t N, size_t Q, size_t K >
  void
  GRingArray<N,Q,K>::identityInit( const size_t index ) {
#pragma omp parallel for
    for ( size_t i = K*index; i < K*(index+1); ++i )
      fmpz_mod_poly_set_coeff_ui( _polys[i], 0, 1 );
    return;
  }

  template< size_t N, size_t Q, size_t K >
  void
  GRingArray<N,Q,K>::identityInit( const size_t start, const size_t end ) {
#pragma omp parallel for
    for ( size_t i = start; i < end; ++i )
      this->identityInit( i );
  }

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>& 
  GRingArray<N,Q,K>::operator=( const GRingArray& other ) {
    _polys = other.data();
    _len = other.len();
    return *this;
  }

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>
  GRingArray<N,Q,K>::operator-( const ring_t& s ) const {
    GRingArray rv{*this};
    rv -= s;
    return rv;
  }

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>
  GRingArray<N,Q,K>::operator-( const GRingArray& other ) const {
    GRingArray rv{*this};
    rv -= other;
    return rv;
  }

  template< size_t N, size_t Q, size_t K >
  void
  GRingArray<N,Q,K>::negGRingSection( const size_t start, const size_t end ) {
    assert( start < end );
    assert( end <= _len );

#pragma omp parallel for
    for ( size_t i = K*start; i < K*end; ++i )
      fmpz_mod_poly_neg( _polys[i], _polys[i] );
  }

  template< size_t N, size_t Q, size_t K >
  void 
  GRingArray<N,Q,K>::operator-=( const ring_t& s ) {
    minusGRingSection( 0, _len, s, *this, 0);
  }

  template< size_t N, size_t Q, size_t K >
  void 
  GRingArray<N,Q,K>::operator-=( const GRingArray& other ) {
    minusGRingSection( 0, _len, other, 0, *this, 0);
  }

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>
  GRingArray<N,Q,K>::operator+( const ring_t& s ) const {
    GRingArray rv{*this};
    rv += s;
    return rv;
  }

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>
  GRingArray<N,Q,K>::operator+( const GRingArray& other ) const {
    GRingArray rv{*this};
    rv += other;
    return rv;
  }

  template< size_t N, size_t Q, size_t K >
  void 
  GRingArray<N,Q,K>::operator+=( const ring_t& s ) {
    plusGRingSection( 0, _len, s, *this, 0);
  }

  template< size_t N, size_t Q, size_t K >
  void 
  GRingArray<N,Q,K>::operator+=( const GRingArray& other ) {
    plusGRingSection( 0, _len, other, 0, *this, 0);
  }

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>
  GRingArray<N,Q,K>::operator*( const ring_t& s ) const {
    GRingArray rv{*this};
    rv *= s;
    return rv;
  }

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>
  GRingArray<N,Q,K>::operator*( const GRingArray& other ) const {
    GRingArray rv{*this};
    rv *= other;
    return rv;
  }

  template< size_t N, size_t Q, size_t K >
  void 
  GRingArray<N,Q,K>::operator*=( const ring_t& s ) {
    mulmodGRingSection( 0, _len, s, *this, 0 );
  }

  template< size_t N, size_t Q, size_t K >
  void 
  GRingArray<N,Q,K>::operator*=( const GRingArray& other ) {
    mulmodGRingSection( 0, _len, other, 0, *this, 0 );
  }

  template< size_t N, size_t Q, size_t K >
  void 
  GRingArray<N,Q,K>::minusGRingSection( const size_t start, const size_t end, const ring_t& 
      otherPoly, GRingArray& target, const size_t targetStart ) {
    assert( end <= _len );
    assert( start < end );
    assert( targetStart + (end - start) <= target.len() );

    ring_t* targetPolys = target.data();
#pragma omp parallel for
    for ( size_t i = 0; i < K*(end - start); ++i )
      fmpz_mod_poly_sub( targetPolys[K*(targetStart - start) + i], _polys[i], 
          otherPoly );
  }

  template< size_t N, size_t Q, size_t K >
  void 
  GRingArray<N,Q,K>::minusGRingSection( const size_t start, const size_t end, const GRingArray& 
      other, const size_t otherStart, GRingArray& target, const size_t 
      targetStart ) {
    assert( end <= _len );
    assert( start < end );
    assert( otherStart + (end - start) <= other.len() );
    assert( targetStart + (end - start) <= target.len() );

    ring_t* otherPolys = other.data();
    ring_t* targetPolys = target.data();
#pragma omp parallel for
    for ( size_t i = 0; i < K*(end - start); ++i )
      fmpz_mod_poly_sub( targetPolys[K*targetStart + i], _polys[K*start + i], 
          otherPolys[K*otherStart + i] );
  }

  template< size_t N, size_t Q, size_t K >
  void 
  GRingArray<N,Q,K>::plusGRingSection( const size_t start, const size_t end, const ring_t& 
      otherPoly, GRingArray& target, const size_t targetStart ) {
    assert( end <= _len );
    assert( start < end );
    assert( targetStart + (end - start) <= target.len() );

    ring_t* targetPolys = target.data();
#pragma omp parallel for
    for ( size_t i = 0; i < K*(end - start); ++i )
      fmpz_mod_poly_add( targetPolys[K*targetStart + i], _polys[K*start + i], 
          otherPoly );
  }

  template< size_t N, size_t Q, size_t K >
  void 
  GRingArray<N,Q,K>::plusGRingSection( const size_t start, const size_t end, const GRingArray& 
      other, const size_t otherStart, GRingArray& target, const size_t 
      targetStart ) {
    assert( end <= _len );
    assert( start < end );
    assert( otherStart + (end - start) <= other.len() );
    assert( targetStart + (end - start) <= target.len() );

    ring_t* otherPolys = other.data();
    ring_t* targetPolys = target.data();
#pragma omp parallel for
    for ( size_t i = 0; i < K*(end - start); ++i )
      fmpz_mod_poly_add( targetPolys[K*targetStart + i], _polys[K*start + i], 
          otherPolys[K*otherStart + i] );
  }

  template< size_t N, size_t Q, size_t K >
  void 
  GRingArray<N,Q,K>::mulmodGRingSection( const size_t start, const size_t end, const ring_t& 
      otherPoly, GRingArray& target, const size_t targetStart ) {
    assert( end <= _len );
    assert( start < end );
    assert( targetStart + (end - start) <= target.len() );

    ring_t* targetPolys = target.data();
    ring_t& F = FF();
#pragma omp parallel for
    for ( size_t i = 0; i < K*(end - start); ++i )
      fmpz_mod_poly_mulmod( targetPolys[K*targetStart + i], _polys[K*start + i], 
          otherPoly, F );
  }


  template< size_t N, size_t Q, size_t K >
  void 
  GRingArray<N,Q,K>::mulmodGRingSection( const size_t start, const size_t end, 
      const GRingArray& other, const size_t otherStart, GRingArray& target, 
      const size_t targetStart ) {
    assert( end <= _len );
    assert( start < end );
    assert( otherStart + (end - start) <= other.len() );
    assert( targetStart + (end - start) <= target.len() );

    ring_t* otherPolys = other.data();
    ring_t* targetPolys = target.data();
    ring_t& F = FF();
#pragma omp parallel for
    for ( size_t i = 0; i < K*(end - start); ++i )
      fmpz_mod_poly_mulmod( targetPolys[K*targetStart + i], _polys[K*start + i], 
          otherPolys[K*otherStart + i], F );
  }

  template< size_t N, size_t Q, size_t K >
  void
  GRingArray<N,Q,K>::dualRegevEncode() {
    ring_t s;
    makeUniformPoly( s );
    dualRegevEncode( s );
    fmpz_mod_poly_clear( s );
  }

  template< size_t N, size_t Q, size_t K >
  void
  GRingArray<N,Q,K>::dualRegevEncode( ring_t& s ) {
    *this *= s;
#pragma omp parallel for
    for ( size_t i = 0; i < K*_len; ++i )
      ternaryPerturb( _polys[i] );
  }

  template< size_t N, size_t Q, size_t K >
  void
  GRingArray<N,Q,K>::dualRegevEncrypt( ring_t* s, std::string& mu ) {
    size_t blocksize = N/8;
    size_t randbuff = 0;
    ring_t& F = FF();
#pragma omp parallel for
    for ( size_t i = 0; i < K*_len; ++i ) {
      fmpz_mod_poly_mulmod( _polys[i], _polys[i], *s, F );
      fmpz_t coeff, cj;
      fmpz_init( coeff ); fmpz_init( cj );
      size_t byte, bit;
      for ( size_t j = 0; j < blocksize; ++j ) {
        if ( j % 4 == 0 ) 
          randbuff = rand();
        byte = mu[blocksize*i + j];
        for ( size_t k = 0; k < 8; ++k ) {
          bit = byte & 1;
          fmpz_set_ui(coeff, bit *(Q/2) );
          fmpz_mod_poly_get_coeff_fmpz( cj, _polys[i], j );
          fmpz_add( cj, cj, coeff );
          fmpz_set_ui( cj, (randbuff&1) + ((randbuff>>1)&1)*(Q-1) );
          fmpz_add( cj, cj, coeff );
          fmpz_mod_poly_set_coeff_fmpz( _polys[i], j, coeff );
          byte >>= 1;
          randbuff >>= 2;
        }
      }
    }
  }

  template< size_t N, size_t Q, size_t K >
  ring_t*
  GRingArray<N,Q,K>::dualRegevEncrypt( std::string& mu ) {
    ring_t* s = new ring_t[1];
    makeUniformPoly( *s );
    dualRegevEncrypt( s, mu );
    return s;
  }

  /* Instead of creating a random ternary error polynomial and adding, we can 
   * avoid creating ephemeral polynomials by perturbing each coefficient 
   * individually.
   */

  template< size_t N, size_t Q, size_t K >
  void 
  GRingArray<N,Q,K>::ternaryPerturb( ring_t& ri ) {
    fmpz_t coeff, eij;
    fmpz_init( eij ); fmpz_init( coeff );

    size_t randbuff = 0;
    for ( size_t j = 0; j < N; ++j ) {
      if ( j % 32 == 0 ) randbuff = rand();
      fmpz_set_ui( eij, (randbuff&1) + ((randbuff>>1)&1)*(Q-1) );
      fmpz_mod_poly_get_coeff_fmpz( coeff, ri, j );
      fmpz_add( coeff, coeff, eij );
      fmpz_mod_poly_set_coeff_fmpz( ri, j, coeff );
      randbuff >>= 2;
    }

    fmpz_clear( coeff ); fmpz_clear( eij );
  }

  template< size_t N, size_t Q, size_t K >
  void
  GRingArray<N,Q,K>::makeUniformPoly( ring_t& r ) {
    fmpz_t q_z;
    fmpz_init_set_ui( q_z, Q );
    fmpz_mod_poly_init2( r, q_z, N );
    for ( size_t i = 0; i < N; ++i )
      fmpz_mod_poly_set_coeff_ui( r, i, rand() );
    fmpz_clear( q_z );
  }

  template< size_t N, size_t Q, size_t K >
  ring_t&
  GRingArray<N,Q,K>::makeF() const {
    fmpz_t q_z;
    fmpz_init_set_ui( q_z, Q );
    ring_t* F = new ring_t[1];
    fmpz_mod_poly_init2( F[0], q_z, N );
    fmpz_mod_poly_set_coeff_ui( F[0], 0, 1 );
    fmpz_mod_poly_set_coeff_ui( F[0], N, 1 );
    return F[0];
  }
     
  float 
  erfinv( float x ) {
    float w, p;
    w = - log((1.0f-x)*(1.0f+x));
    if ( w < 5.000000f ) {
      w = w - 2.500000f;
      p = 2.81022636e-08f;
      p = 3.43273939e-07f + p*w;
      p = -3.5233877e-06f + p*w;
      p = -4.39150654e-06f + p*w;
      p = 0.00021858087f + p*w;
      p = -0.00125372503f + p*w;
      p = -0.00417768164f + p*w;
      p = 0.246640727f + p*w;
      p = 1.50140941f + p*w;
    } else {
      w = sqrtf(w) - 3.000000f;
      p = -0.000200214257f;
      p = 0.000100950558f + p*w;
      p = 0.00134934322f + p*w;
      p = -0.00367342844f + p*w;
      p = 0.00573950773f + p*w;
      p = -0.0076224613f + p*w;
      p = 0.00943887047f + p*w;
      p = 1.00167406f + p*w;
      p = 2.83297682f + p*w;
    }
    return p*x;
  }

} // namespace gring

#endif

