#ifndef _ABE_H_
#define _ABE_H_

#include <ring/GRingArray.hpp>
#include <ring/Trapdoor.hpp>

#include <sstream>
#include <vector>

namespace gring {

  class ABEId {
  private:
    std::vector< size_t > _x;
    size_t _L;
  public:
    ABEId( std::vector< size_t > id ) : _x{std::vector< size_t >{id.size()}}, 
        _L{id.size()} { _x = id; }
    std::vector< size_t > x() const { return _x; }
    size_t L() const { return _L; }
    void print() const;
  };

  void
  ABEId::print() const {
    std::cout << "ABEId: ";
    for ( size_t i = 0; i < this->x().size(); ++i )
      std::cout << this->x()[i] << " ";
    std::cout << " len: " << this->L() << std::endl;
  }

  template< size_t N, size_t Q, size_t K >
  class ABECtxt {
  private:
    GRingArray<N,Q,K>* _ctxt;
    ABEId* _id;
    ABECtxt() = delete;
  public:
    ABECtxt( ABEId& ID, const size_t l ) : _ctxt{new GRingArray<N,Q,K>{l}}, 
      _id{new ABEId{ID}} {}
    ABECtxt( const ABECtxt& other ) : 
      _ctxt{new GRingArray<N,Q,K>{*other.rings()}}, _id{other.id()} {}
    ABECtxt& operator=( const ABECtxt& other ) { 
      _ctxt = other.rings(); _id = other.id(); }
    ~ABECtxt() { delete _ctxt; }
    GRingArray<N,Q,K>* rings() const { return _ctxt; }
    inline ABEId* id() const { return this->_id; }
    inline size_t L() const { return _ctxt->len(); }
  };


  template< size_t N, size_t Q, size_t K >
  class ABEGate {
  private:
    std::vector< size_t > _weights;
    size_t _k;
    ABEGate() = delete;
  public:
    ABEGate( std::vector< size_t > w ) : _weights{w}, _k{w.size()} {}
    ABEGate( const ABEGate<N,Q,K>& other ) : _weights{other.weights()}, 
      _k{other.weights().size()} {}
    ABEGate<N,Q,K>* operator=( const ABEGate<N,Q,K>& other ) { 
      _weights = other.weights();
      _k = other.k(); 
    }
    virtual ~ABEGate() {}
    std::vector< size_t > weights() const { return _weights; }
    size_t k() const { return _k; }
    virtual GRingArray<N,Q,K>* evalPK( GRingArray<N,Q,K>* B ) = 0;
    virtual GRingArray<N,Q,K>* evalCt( GRingArray<N,Q,K>* B, ABECtxt<N,Q,K>* 
        ctxt ) = 0;
    virtual size_t evalId( ABEId& id ) = 0;
    virtual void print() const = 0;
  };

  template< size_t N, size_t Q, size_t K >
  class ABEAddGate : public ABEGate<N,Q,K> {
  public:
    ABEAddGate( std::vector< size_t > w ) : ABEGate<N,Q,K>{w} { assert( w.size() > 0 );}
    ~ABEAddGate() {}
    GRingArray<N,Q,K>* evalPK( GRingArray<N,Q,K>* B );
    GRingArray<N,Q,K>* evalCt( GRingArray<N,Q,K>* B, ABECtxt<N,Q,K>* ctxt );
    size_t evalId( ABEId& id );
    void print() const;
  };

  template< size_t N, size_t Q, size_t K >
  class ABEMulGate : public ABEGate<N,Q,K> {
  public:
    ABEMulGate( std::vector< size_t > w ) : ABEGate<N,Q,K>{w} 
      { assert( w.size() == 1 ); }
    ~ABEMulGate() {}
    GRingArray<N,Q,K>* evalPK( GRingArray<N,Q,K>* B );
    GRingArray<N,Q,K>* evalCt( GRingArray<N,Q,K>* B, ABECtxt<N,Q,K>* ctxt );
    GRingArray<N,Q,K>* calcR( GRingArray<N,Q,K>* B );
    size_t evalId( ABEId& id );
    void print() const;
  };

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>*
  ABEAddGate<N,Q,K>::evalPK( GRingArray<N,Q,K>* B ) {
    assert( this->weights().size() == B->len()/2 );

    size_t length = this->k();
    GRingArray<N,Q,K> R{4*K*length};
    solveRWeighted( &R, this->weights() );

    for ( size_t i = 0; i < length; ++i )
      for ( size_t j = 0; j < 2*K; ++j )
        B->mulmodGRingSection( 2*i, 2*(i+1), R, 4*K*i + 2*j, R, 4*K*i + 2*j );

    auto Bg = new GRingArray<N,Q,K>{2};
    for ( size_t i = 0; i < length; ++i )
      for ( size_t j = 0; j < 2*K; ++j )
        for ( size_t l = 1; l < 2*K; ++l )
          fmpz_mod_poly_add( R.data()[4*K*K*i + 2*K*j], R.data()[4*K*K*i + 2*K*j], 
              R.data()[4*K*K*i + 2*K*j + l] );
    for ( size_t i = 0; i < 2*K; ++i )
      for ( size_t j = 0; j < length; ++j )
        fmpz_mod_poly_add( Bg->data()[i], Bg->data()[i], R.data()[2*K*i + 4*K*K*j] );

    return Bg;
  }

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>*
  ABEAddGate<N,Q,K>::evalCt( GRingArray<N,Q,K>*, ABECtxt<N,Q,K>* c ) {
    auto cg = GRingArray<N,Q,K>{2*c->id()->x().size()};
    cg.plusGRingSection( 0, 2*c->id()->x().size(), *c->rings(), 2, cg, 0 );

    return evalPK( &cg );
  }

  template< size_t N, size_t Q, size_t K >
  size_t
  ABEAddGate<N,Q,K>::evalId( ABEId& id ) {
    assert( this->weights().size() == id.L() );

    std::cout << "ABEAddGate::evalId: ";

    std::cout << this->weights()[0] << "*" << id.x()[0];
    size_t fx = this->weights()[0] * id.x()[0];;
    for ( size_t i = 1; i < id.L(); ++i ) {
      std::cout << " + " << this->weights()[i] << "*" << id.x()[i];
      fx += this->weights()[i] * id.x()[i];
    }
    fx %= Q;
    std::cout << " = " << fx << std::endl;
    
    return fx;
  }

  template< size_t N, size_t Q, size_t K >
  void
  ABEAddGate<N,Q,K>::print() const {
    std::cout << "ABEAddGate: ";
    for ( size_t i = 0; i < this->weights().size(); ++i )
      std::cout << this->weights()[i] << " ";
    std::cout << std::endl;
  }
  
  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>*
  ABEMulGate<N,Q,K>::evalPK( GRingArray<N,Q,K>* B ) {
    size_t length = B->len()/2;
    auto R = this->calcR( B );

    GRingArray<N,Q,K> temp{4*K};
    GRingArray<N,Q,K> target{2};
    auto Bg = new GRingArray<N,Q,K>{2};
    for ( size_t j = 0; j < 2*K; ++j ) {
      B->mulmodGRingSection( 2*length-2, 2*length, *R, 4*K*(length-1) + 2*j, temp, 2*j );
      for ( size_t l = 0; l < 2*K; ++l )
        fmpz_mod_poly_add( Bg->data()[j], Bg->data()[j], temp.data()[2*K*j + l] );
    }

    return Bg;
  }

  template< size_t N, size_t Q, size_t K>
  GRingArray<N,Q,K>*
  ABEMulGate<N,Q,K>::calcR( GRingArray<N,Q,K>* B ) {
    size_t length = B->len()/2;

    auto R = new GRingArray<N,Q,K>{4*K*length};
    GRingArray<N,Q,K> temp{4*K};
    GRingArray<N,Q,K> target{2};

    solveRWeightedSection( R, this->weights()[0] );
    for ( size_t i = 1; i < length; ++i ) {
      for ( size_t j = 0; j < 2*K; ++j ) {
        B->mulmodGRingSection( 2*(i-1), 2*i, *R, 4*K*(i-1) + 2*j, temp, 2*j );
        fmpz_mod_poly_set( target.data()[j], temp.data()[2*K*j] );
        for ( size_t l = 1; l < 2*K; ++l )
          fmpz_mod_poly_add( target.data()[j], 
            target.data()[j], temp.data()[2*K*j + l] );
        fmpz_mod_poly_neg( target.data()[j], target.data()[j] );
      }
      solveRTarget( R, &target, i );
    }

    return R;
  }

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>*
  ABEMulGate<N,Q,K>::evalCt( GRingArray<N,Q,K>* B, ABECtxt<N,Q,K>* c ) {
    assert( B->len() == c->rings()->len()-4 );

    size_t length = B->len()/2;
    auto R = this->calcR( B );
    auto BRtemp = new GRingArray<N,Q,K>{R->len()};

    for ( size_t i = 0; i < length; ++i )
      for ( size_t j = 0; j < 2*K; ++j )
        c->rings()->mulmodGRingSection( 2*i+2, 2*i+4, *R, 4*K*i + 2*j, *BRtemp, 4*K*i + 2*j );

    for ( size_t i = 0; i < length; ++i )
      for ( size_t j = 0; j < 2*K; ++j )
        for ( size_t l = 1; l < 2*K; ++l )
          fmpz_mod_poly_add( BRtemp->data()[4*K*K*i + 2*K*j], BRtemp->data()[4*K*K*i + 2*K*j], 
              BRtemp->data()[4*K*K*i + 2*K*j + l] );

    auto Cg = new GRingArray<N,Q,K>{2};
    size_t coeff = 1;
    fmpz_t coeff_z;
    fmpz_init_set_ui( coeff_z, coeff );
    for ( size_t j = length-1; j < length; --j ) {
      for ( size_t i = 0; i < 2*K; ++i ) {
        fmpz_mod_poly_scalar_mul_fmpz( BRtemp->data()[2*K*i + 4*K*K*j], BRtemp->data()[2*K*i + 4*K*K*j], coeff_z );
        fmpz_mod_poly_add( Cg->data()[i], Cg->data()[i], BRtemp->data()[2*K*i + 4*K*K*j] );
      }
      
      coeff *= c->id()->x()[j];
      coeff %= Q;
      fmpz_set_ui( coeff_z, coeff );
    }

    GRingArray<N,Q,K> temp{4*K};
    auto BR = new GRingArray<N,Q,K>{2};
    for ( size_t j = 0; j < 2*K; ++j ) {
      B->mulmodGRingSection( 2*length-2, 2*length, *R, 4*K*(length-1) + 2*j, temp, 2*j );
      for ( size_t l = 0; l < 2*K; ++l )
        fmpz_mod_poly_add( BR->data()[j], BR->data()[j], temp.data()[2*K*j + l] );
    }

    return Cg;
  }

  template< size_t N, size_t Q, size_t K >
  size_t
  ABEMulGate<N,Q,K>::evalId( ABEId& id ) {
    std::cout << "ABEMulGate::evalId: " << this->weights()[0];
    size_t eval = this->weights()[0];
    for ( size_t i = 0; i < id.x().size(); ++i ) {
      std::cout << " * " << id.x()[i];
      eval *= id.x()[i];
    }
    eval %= Q;
    std::cout << " = " << eval << std::endl;
    
    return eval;
  }

  template< size_t N, size_t Q, size_t K >
  void
  ABEMulGate<N,Q,K>::print() const {
    std::cout << "ABEMulGate: ";
    for ( size_t i = 0; i < this->weights().size(); ++i )
      std::cout << this->weights()[i] << " ";
    std::cout << std::endl;
  }

  template< size_t N, size_t Q, size_t K >
  class ABEEvalCiphertext {
  private:
    GRingArray<N,Q,K>* _evalCtxt;
  public:
    ABEEvalCiphertext() : _evalCtxt{new GRingArray<N,Q,K>{4}} {}
    ABEEvalCiphertext( const ABEEvalCiphertext& other ) : 
      _evalCtxt{new GRingArray<N,Q,K>{*other.rings()}} {}
    ABEEvalCiphertext& operator=( const ABEEvalCiphertext& other ) { 
      _evalCtxt = other.rings(); }
    ~ABEEvalCiphertext() { delete _evalCtxt; }
    GRingArray<N,Q,K>* rings() const { return _evalCtxt; }
  };

  template< size_t N, size_t Q, size_t K >
  class ABEPublic {
  private:
    GRingArray<N,Q,K>* _pk;
    GRingArray<N,Q,K>* _ctxtKeys;
    GRingArray<N,Q,K>* _d;
    size_t _L;
    ABEPublic() = delete;
  public:
    ABEPublic( const size_t len ) : 
      _pk{new GRingArray<N,Q,K>{2}}, _ctxtKeys{new GRingArray<N,Q,K>{2*len}},
      _d{new GRingArray<N,Q,K>{2}}, _L{len} {}
    ABEPublic( const GRingArray<N,Q,K>* pkPolys, const GRingArray<N,Q,K>* 
        ctxtPolys, const GRingArray<N,Q,K>* dPolys, const size_t len ) : 
      _pk{new GRingArray<N,Q,K>{*pkPolys}}, _ctxtKeys{new 
        GRingArray<N,Q,K>{*ctxtPolys}}, _d{new GRingArray<N,Q,K>{*dPolys}}, 
      _L{len} {}
    ABEPublic( const ABEPublic& other ) : 
      _pk{new GRingArray<N,Q,K>{*other.pkRings()}}, 
      _ctxtKeys{other.ctxtRings()}, _d{other.D()}, _L{other.L()} {}
    ABEPublic& operator=( const ABEPublic& other ) { 
      _pk = other.pkRings();
      _ctxtKeys = other.ctxtKeys();
      _d = other.D();
      _L = other.L();
    }
    virtual ~ABEPublic() { delete _pk; delete _ctxtKeys; }
    GRingArray<N,Q,K>* pkRings() const { return _pk; }
    GRingArray<N,Q,K>* ctxtRings() const { return _ctxtKeys; }
    GRingArray<N,Q,K>* D() const { return _d; }
    size_t L() const { return _L; }
    ABECtxt<N,Q,K>* encrypt( ABEId& id, std::string msg );
  };

  template< size_t N, size_t Q, size_t K >
  ABECtxt<N,Q,K>*
  ABEPublic<N,Q,K>::encrypt( ABEId& id, std::string message ) {
    assert( this->L() == id.x().size() );
    fmpz_t q_z; fmpz_init_set_ui( q_z, Q );
    ring_t s;
    fmpz_mod_poly_init2( s, q_z, N );

    auto e0 = GRingArray<N,Q,K>{2};
    e0.ternaryInit();
    auto e1 = GRingArray<N,Q,K>{2};
    e1.ternaryInit();

    auto ctxt =  new ABECtxt<N,Q,K>{id, 2*(this->L()+2)};
    auto H = ctxt->rings();
    H->plusGRingSection( 0, 2, *this->pkRings(), 0, *H, 0 );
    H->plusGRingSection( 2, 2*(this->L()+1), *this->ctxtRings(), 0, *H, 2 );
    H->plusGRingSection( 2*(this->L()+1), 2*(this->L()+2), *this->D(), 0, *H, 
      2*(this->L()+1) );

    fmpz_t coeff_z, xj_z;
    fmpz_init( coeff_z ); 
    fmpz_init( xj_z );

    for ( size_t i = 0; i < this->L(); ++i )  {
      for ( size_t j = 0; j < 2*K; ++j ) {
        fmpz_set_ui( coeff_z, (id.x()[i] * (size_t(1) << j)) % Q );
        fmpz_mod_poly_get_coeff_fmpz( xj_z, H->data()[2*K*(i+1) + j], 0 );
        fmpz_add( xj_z, xj_z, coeff_z );
        fmpz_mod_poly_set_coeff_fmpz( H->data()[2*K*(i+1) + j], 0, xj_z );
      }
    }
    fmpz_clear( xj_z );
    fmpz_clear( coeff_z );

    H->makeUniformPoly( s );
    *H *= s;

    auto eSim = GRingArray<N,Q,K>{4*K*(this->L()+1)};
    eSim.identityInit( 0, 4*K );
    eSim.ternaryInit( 4*K, 4*K*(this->L()+1) );
    std::cout << "parallel for ..." << std::endl;
    for ( size_t i = 0; i < 2*K*(this->L()+1); ++i ) {
      e0.mulmodGRingSection( 0, 2, eSim, 2*i, eSim, 2*i );
      for ( size_t j = 1; j < 2*K; ++j )
        fmpz_mod_poly_add( eSim.data()[2*K*i], eSim.data()[2*K*i], 
            eSim.data()[2*K*i + j] );
      fmpz_mod_poly_add( H->data()[2*K + i], H->data()[2*K + i], eSim.data()[2*K*i] );
    }

    /*
    std::cout << "e0: " << std::endl;
    fmpz_mod_poly_print( e0.data()[0] );
    std::cout << std::endl;

    std::cout << "Rf: " << std::endl;
    fmpz_mod_poly_print( skF->data()[0] );
    std::cout << std::endl;

    ring_t first; fmpz_mod_poly_init2( first, q_z, N );
    fmpz_mod_poly_mulmod( first, skF->data()[0], e0.data()[0], skF->FF() );
    std::cout << "first: " << std::endl;
    fmpz_mod_poly_print( first );
    std::cout << std::endl;

    std::cout << "skF: " << std::endl;
    for ( size_t i = 0; i < K; ++i ) {
      std::cout << "poly " << i << std::endl;
      fmpz_mod_poly_print( skF->data()[i] );
      std::cout << std::endl;
    }

    auto Rfe0 = GRingArray<N,Q,K>{2};
    skF->mulmodGRingSection( 0, 2, e0, 0, Rfe0, 0 );

    std::cout << "re: " << std::endl;
    for ( size_t i = 0; i < K; ++i ) {
      std::cout << "poly " << i << std::endl;
      fmpz_mod_poly_print( Rfe0.data()[i] );
      std::cout << std::endl;
    }

    for ( size_t i = 1; i < 2*K; ++i )
      fmpz_mod_poly_add( Rfe0.data()[0], Rfe0.data()[0], Rfe0.data()[i] );
    std::cout << "Rfe0:" << std::endl;
    fmpz_mod_poly_print( Rfe0.data()[0] );
    std::cout << std::endl;


    std::cout << "ctxt A: " << std::endl;
    fmpz_mod_poly_print( H->data()[0] );
    std::cout << std::endl;
    */

    e0.plusGRingSection( 0, 2, *H, 0, *H, 0 );
    e1.plusGRingSection( 0, 2, *H, 2*(this->L()+1), *H, 2*(this->L()+1) );

    fmpz_mod_poly_clear( s );
    
    ring_t mu; 

    for ( size_t i = 0; i < 2*K; ++i ) {
      if ( (N/8)*i > message.size() ) break;
      fmpz_mod_poly_init2( mu, q_z, N );
//      std::cout << "msg_" << i << " " << message.substr( (N/8)*i, (N/8) ) << std::endl;
      encode<N,Q,K>( mu, message.substr( (N/8)*i, (N/8)*(i+1) ) );
      fmpz_mod_poly_add( H->data()[2*(this->L()+1)*K + i], H->data()[2*(this->L()+1)*K + i], mu );
    }

    fmpz_mod_poly_clear( mu );
    fmpz_clear( q_z );

    
    return ctxt;
  }

  template< size_t N, size_t Q, size_t K >
  class ABESecret : public ABEPublic<N,Q,K> {
  private:
    GRingArray<N,Q,K>* _skF;
    ABEGate<N,Q,K>* _g;
    ABESecret() = delete;
  public:  
    ABESecret( const ABEPublic<N,Q,K>* pub, const GRingArray<N,Q,K>* sk, 
      ABEGate<N,Q,K>* gate ) : ABEPublic<N,Q,K>{*pub}, 
      _skF{new GRingArray<N,Q,K>{*sk}}, _g{gate} {}
    ABESecret( const ABESecret& other ) : 
      _skF{new GRingArray<N,Q,K>{other.skRings()}}, _g{other.g()} {}
    ABESecret& operator=( const ABESecret& other ) { _skF = other.skRings(); 
      _g = other.g(); }
    ~ABESecret() { delete _skF; delete _g; }
    GRingArray<N,Q,K>* skRings() const { return _skF; }
    ABEGate<N,Q,K>* g() const { return _g; }
    std::string decrypt( ABECtxt<N,Q,K>* ctxt ) const;
  };

  template< size_t N, size_t Q, size_t K >
  std::string 
  ABESecret<N,Q,K>::decrypt( ABECtxt<N,Q,K>* ctxt ) const {
    if ( this->g()->evalId( *ctxt->id() ) ) return "invalid id";

    auto Cg = this->g()->evalCt( this->ctxtRings(), ctxt );
    auto eval = GRingArray<N,Q,K>{this->skRings()->len()};
    for ( size_t i = 0; i < this->skRings()->len()/4; ++i ) {
      ctxt->rings()->mulmodGRingSection( 0, 2, *this->skRings(), 2*i, eval, 4*i );
      Cg->mulmodGRingSection( 0, 2, *this->skRings(), this->skRings()->len()/2 + 2*i, eval, 4*i + 2 );
    }

    auto final = GRingArray<N,Q,K>{2};
    for ( size_t i = 0; i < 2*K; ++i )
      for ( size_t j = 0; j < 4*K; ++j )
        fmpz_mod_poly_add( final.data()[i], final.data()[i], eval.data()[4*K*i + j] );
    ctxt->rings()->minusGRingSection( ctxt->L()-2, ctxt->L(), final, 0, final, 0 );


    /*
    std::cout << "final: " << std::endl;
    fmpz_mod_poly_print( final.data()[0] );
    std::cout << std::endl;
    */
    
    std::stringstream ss;
    for ( size_t i = 0; i < 2*K; ++i ) {
//      std::cout << "dec_" << i << " " << decode<N,Q,K>( final.data()[i] ) << std::endl;
      ss << decode<N,Q,K>( final.data()[i] );
    }

    return ss.str();
  }

  template< size_t N, size_t Q, size_t K >
  class ABEMasterSecret : public ABEPublic<N,Q,K> {
  private:
    GRingArray<N,Q,K>* _msk;
    long* _samples;
    ABEMasterSecret() = delete;
  public:
    ABEMasterSecret( const size_t len, long* samples );
    ABEMasterSecret( const ABEMasterSecret& other ) : 
      ABEPublic<N,Q,K>{other.pkRings(), other.ctxtRings(), other.L()}, 
      _msk{new GRingArray<N,Q,K>{*other.mskRings()}}, _samples{other.samples()} {}
    ABEMasterSecret& operator=( const ABEMasterSecret& other );
    ~ABEMasterSecret() { delete _msk; }
    GRingArray<N,Q,K>* mskRings() const { return _msk; }
    ABEPublic<N,Q,K>* publicKey() const { 
      return new ABEPublic<N,Q,K>{this->pkRings(), this->ctxtRings(), this->D(),
        this->L()}; 
    }
    long* samples() const { return _samples; }
    ABESecret<N,Q,K>* keyGen( ABEGate<N,Q,K>* f ) const;
  };

  template< size_t N, size_t Q, size_t K >
  ABEMasterSecret<N,Q,K>::ABEMasterSecret( const size_t len, long* samps ) : 
      ABEPublic<N,Q,K>{len}, _msk{new GRingArray<N,Q,K>{2}}, _samples{samps} {
    trapGen( this->pkRings(), _msk, 0 );
    this->ctxtRings()->uniformInit();
    this->D()->uniformInit();
  }

  template< size_t N, size_t Q, size_t K >
  ABEMasterSecret<N,Q,K>&
  ABEMasterSecret<N,Q,K>::operator=( const ABEMasterSecret& other ) {
    this->_pk = other.pkRings();
    this->_ctxtKeys = other.ctxtRings();
    this->_L = other.L();
    this->_msk = other.mskRings();
    this->_samples = other.samples();
  }

  template< size_t N, size_t Q, size_t K >
  ABESecret<N,Q,K>*
  ABEMasterSecret<N,Q,K>::keyGen( ABEGate<N,Q,K>* g ) const {
    auto Bg = g->evalPK( this->ctxtRings() );
    auto Rg = extendRightSampleD( this->mskRings(), this->pkRings(), Bg, 
        this->D(), this->samples() );
    
    return new ABESecret<N,Q,K>{this, Rg, g};
  }

  template< size_t N, size_t Q, size_t K >
  class ABECircuit {
  private:
    std::vector< ABEGate<N,Q,K>* > _gates;
    std::vector< std::vector< size_t > > _wires;
    size_t _L;
    ABECircuit() = delete;
  public:
    ABECircuit( std::vector< ABEGate<N,Q,K>* > gs, 
        std::vector< std::vector< size_t > > ws, const size_t l ) : _gates{gs}, 
        _wires{ws}, _L{l} {}
    ABECircuit( const ABECircuit<N,Q,K>& other ) : _gates{other.gates()}, 
      _wires{other.wires()}, _L{other.L()} {}
    ABECircuit<N,Q,K>* operator=( const ABECircuit<N,Q,K>& other ) { 
      _gates = other.gates();
      _wires = other.wires(); 
      _L = other.L();
    }
    ~ABECircuit() {}
    std::vector< ABEGate<N,Q,K>* > gates() const { return _gates; }
    std::vector< std::vector< size_t > > wires() const { return _wires; }
    size_t L() const { return _L; }
    GRingArray<N,Q,K>* evalPK( GRingArray<N,Q,K>* B );
    GRingArray<N,Q,K>* evalCt( GRingArray<N,Q,K>* B, ABECtxt<N,Q,K>* ctxt );
    size_t evalId( ABEId& id );
  };

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>*
  ABECircuit<N,Q,K>::evalPK( GRingArray<N,Q,K>* B ) {

  }

  template< size_t N, size_t Q, size_t K >
  GRingArray<N,Q,K>*
  ABECircuit<N,Q,K>::evalCt( GRingArray<N,Q,K>*B, ABECtxt<N,Q,K>* ctxt ) {

  }

  template< size_t N, size_t Q, size_t K >
  size_t
  ABECircuit<N,Q,K>::evalId( ABEId& id ) {
    auto endGate = _gates.back();
    auto endWire = _wires.back();
    std::vector< size_t > endInput;
    for ( size_t i = 0; i < endWire.size(); ++i ) {
      if ( endWire[i] < _L ) endInput.push_back( id.x()[endWire[i]] );
    }
    auto inputId = ABEId{endInput};
    return endGate->evalId( inputId );
  }

} // namespace gring

#endif // _ABE_H_
