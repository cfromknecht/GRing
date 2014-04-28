#ifndef _IBE_H_
#define _IBE_H_

#include <ring/GRingArray.hpp>
#include <ring/Trapdoor.hpp>

namespace gring {

  /**
   * IBE Ciphertext
   */

  template< size_t N, size_t Q, size_t K >
  class IBECiphertext {
  private:
    GRingArray<N,Q,K>* _ctxt;
    ring_t* _ys;
  public:
    IBECiphertext() : _ctxt{new GRingArray<N,Q,K>{2}}, _ys{nullptr} {}
    IBECiphertext( const IBECiphertext<N,Q,K>& other ) : 
      _ctxt{new GRingArray<N,Q,K>{*other.ctxt()}}, _ys{nullptr} {}
    IBECiphertext( GRingArray<N,Q,K>* polys ) : 
      _ctxt{new GRingArray<N,Q,K>{*polys}}, _ys{nullptr} {}
    IBECiphertext& operator=( const IBECiphertext<N,Q,K>& other );
    ~IBECiphertext() { delete _ctxt; fmpz_mod_poly_clear( *_ys ); delete [] _ys; }
    GRingArray<N,Q,K>* rings() const { return _ctxt; }
    ring_t* YS() const { return _ys; }
    void setYS( ring_t* ys ) { _ys = ys; }
  };

  template< size_t N, size_t Q, size_t K >
  IBECiphertext<N,Q,K>&
  IBECiphertext<N,Q,K>::operator=( const IBECiphertext& other ) {
    this->_ctxt = other.ctxt();
    this->_ys = other.ys();
  }

  /**
   * IBE Secret
   */

  template< size_t N, size_t Q, size_t K >
  class IBESecret {
  private:
    GRingArray<N,Q,K>* _skID;
    std::string _ID;
    IBESecret() = delete;
  public:
    IBESecret( std::string id ) : _skID{new GRingArray<N,Q,K>{2}}, _ID{id} {}
    IBESecret( std::string id, GRingArray<N,Q,K>* polys ) : 
      _skID{polys}, _ID{id} {}
    IBESecret( const IBESecret& other ) : IBESecret{other.ID(), 
      new GRingArray<N,Q,K>{other.skID()}} {}
    IBESecret& operator=( const IBESecret& other );
    ~IBESecret() { delete _skID; }
    GRingArray<N,Q,K>* rings() const { return _skID; }
    std::string ID() const { return _ID; }
    std::string decrypt( const IBECiphertext<N,Q,K>* ctxt );
  };

  template< size_t N, size_t Q, size_t K >
  std::string
  IBESecret<N,Q,K>::decrypt( const IBECiphertext<N,Q,K>* ctxt ) {
    fmpz_t q_z;
    fmpz_init_set_ui( q_z, Q );
    ring_t plainRing;
    fmpz_mod_poly_init2( plainRing, q_z, N );

    GRingArray<N,Q,K> A_skID = *ctxt->rings() * *_skID;
    for ( size_t i = 0; i < 2*K; ++i )
      fmpz_mod_poly_add( plainRing, plainRing, A_skID.data()[i] );
    fmpz_mod_poly_sub( plainRing, *ctxt->YS(), plainRing );

    return decode<N,Q,K>( plainRing );
  }

  template< size_t N, size_t Q, size_t K >
  IBESecret<N,Q,K>&
  IBESecret<N,Q,K>::operator=( const IBESecret& other ) {
    this->_skID = other.skID();
    this->_ID = other.ID();
  }

  /**
   * IBE Public
   */

  template< size_t N, size_t Q, size_t K >
  class IBEPublic {
  protected:
    GRingArray<N,Q,K>* _mpk;
  public:
    IBEPublic() : _mpk{new GRingArray<N,Q,K>{2}} {}
    IBEPublic( const IBEPublic& other ) : IBEPublic{other.pkRings()} {}
    IBEPublic( GRingArray<N,Q,K>* polys ) : _mpk{polys} {}
    IBEPublic& operator=( const IBEPublic& other );
    virtual ~IBEPublic() { delete _mpk; }
    GRingArray<N,Q,K>* pkRings() const { return _mpk; }
    IBECiphertext<N,Q,K>* encrypt( std::string id, std::string message );
  };

  template< size_t N, size_t Q, size_t K >
  IBECiphertext<N,Q,K>*
  IBEPublic<N,Q,K>::encrypt( std::string id, std::string message ) {
    fmpz_t q_z;
    fmpz_init_set_ui( q_z, Q );

    ring_t mu;
    fmpz_mod_poly_init2( mu, q_z, N );
    encode<N,Q,K>( mu, message );

    ring_t* y = new ring_t[1];
    fmpz_mod_poly_init2( *y, q_z, N );
    idToY( *y, id, N, K );

    ring_t s;
    fmpz_mod_poly_init2( s, q_z, N );
    _mpk->makeUniformPoly( s );

    IBECiphertext<N,Q,K>* ctxt = new IBECiphertext<N,Q,K>{_mpk};
    ctxt->rings()->dualRegevEncode( s );
    fmpz_mod_poly_mulmod( *y, *y, s, _mpk->FF() );
    fmpz_mod_poly_add( *y, *y, mu );
    ctxt->setYS( y );

    fmpz_mod_poly_clear( s );
    fmpz_clear( q_z );

    return ctxt;
  }

  /**
   * IBE Master Secret Key
   */

  template< size_t N, size_t Q, size_t K >
  class IBEMasterSecret : public IBEPublic<N,Q,K> {
  private:
    GRingArray<N,Q,K>* _msk;
  public:
    IBEMasterSecret();
    IBEMasterSecret( const IBEMasterSecret& other );
    IBEMasterSecret& operator=( const IBEMasterSecret& other );
    ~IBEMasterSecret() { delete _msk; }
    inline GRingArray<N,Q,K>* mskRings() const { return _msk; }
    inline IBEPublic<N,Q,K>* publicKey() const { 
      return new IBEPublic<N,Q,K>{new GRingArray<N,Q,K>{*this->pkRings()}}; }
    IBESecret<N,Q,K>* secretKeyForID( std::string ID ) const;
  };

  template< size_t N, size_t Q, size_t K >
  IBEMasterSecret<N,Q,K>::IBEMasterSecret() : IBEPublic<N,Q,K>{}, 
      _msk{new GRingArray<N,Q,K>{2}} {
    trapGen( this->pkRings(), _msk, 0 );
  }

  template< size_t N, size_t Q, size_t K >
  IBEMasterSecret<N,Q,K>::IBEMasterSecret( const IBEMasterSecret& other ) : 
      IBEPublic<N,Q,K>{other.pkRings()} {
    this->_msk = other.msk();
  }

  template< size_t N, size_t Q, size_t K >
  IBEMasterSecret<N,Q,K>&
  IBEMasterSecret<N,Q,K>::operator=( const IBEMasterSecret& other ) {
    this->_mpk = other.pkRings();
    this->_msk = other.mskRings();
  }

  template< size_t N, size_t Q, size_t K >
  IBESecret<N,Q,K>*
  IBEMasterSecret<N,Q,K>::secretKeyForID( std::string ID ) const {
    fmpz_t q_z;
    fmpz_init_set_ui( q_z, Q );
    
    IBESecret<N,Q,K>* userSecret = new IBESecret<N,Q,K>{ID};
    GRingArray<N,Q,K>* data = userSecret->rings();

    ring_t y;
    fmpz_mod_poly_init2( y, q_z, N );
    invertID( y, ID, data, 1 );
    
    data->mulmodGRingSection( 1, 2, *_msk, 0, *data, 0 );

    fmpz_mod_poly_clear( y );
    fmpz_clear( q_z );

    return userSecret;    
  }

} // namespace gring

#endif // _IBE_H_

