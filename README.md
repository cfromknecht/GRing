GRing
=====
GRing is a lightweight, efficient implementation of G-Trapdoors under the Ring-LWE assumption.  It includes operations for arithmetic, memory management, and trapdoor operations.  The library provides implementations for both Identity-Based Encryption and Attribute-Based Encryption that leverage these cryptographic primitives.
  
Setup
-----
```
git clone https://github.com/cfromknecht/GRing.git
```

Running
-------
The main file is located at src/demo/main.cpp. To execute:
```
make  
./install/bin/ring-demo.x
```
The current demo includes examples for both ABE and IBE schemes.

Parameters
----------
* **N:** default security parameter
* **Q:** integer modulus, must be an even power of 2
* **K:** log_2( Q ) s.t. Q = 2^K
* **SS:** Guassian smoothing parameter, in practice this value need not be above 4
* **L:** length of ABE IDs
  
Values for N, K, SS, and L can be found at the top of the main.cpp file.  Note that we use K to derive Q to ensure that Q is a power of 2, which helps to simplify the G-Trapdoor algorithms.


Generating Guassian Samples
----------------
``` cpp
long* samples = gring::readOrBuildSamples<N,Q,K>( filename, smoothing_parameter );
```
To provide efficiency, GRing uses a finite set from which it consumes independent Guassian samples.  This method will first attempt to load existing samples from ```filename```.  If the file does not exist, contains an insufficient number of samples, or uses a different ```smoothing_parameter```, GRing will build a new sample set and replace ```filename```.


ABE for Arithmetic Circuits
--------------------------
* **Creating a Master Secret Key**
``` cpp
auto msk = new gring::ABEMasterSecret<N,Q,K>{L, samples};
```
where ```L``` is the maximum length of IDs for the system and ```samples``` are generated as described above.

* **Retrieving a Public Key**
``` cpp
auto pk = msk->publicKey();
```

* **Creating Gates**
``` cpp
std::vector< size_t > addWeights{1, 0, 4, 5, 2};
auto gAdd = new gring::ABEAddGate<N,Q,K>{addWeights};
  
std::vector< size_t > mulWeight{1};
auto gMul = new gring::ABEMulGate<N,Q,K>{mulWeight};
  
gAdd->print()       // ABEAddGate: 1 0 4 5 2
gMul->print()       // ABEMulGate: 1
```
The weight vector for an ABEAddGate may be of arbitrary length but is checked at runtime to ensure that it matches the length of the incoming ID or GRingArray.  ABEMulGate requires the weight vector to be of length one but can handle any number of inputs.

* **Generating Secret Keys for Gates**
``` cpp
auto skG = msk->keyGen( g );
```
where ```g``` is either an ABEAddGate or ABEMulGate.

* **Creating IDs**
``` cpp
std::vector< size_t > idWeights{1, 0, Q-1, 5, 3}
auto id = gring::ABEId{ idWeights };
```
* **Encrypting Messages**
``` cpp
std::string msg = "encrypt me";
auto ctxt = pk->encrypt( id, msg );
```
where the length of ```id``` is equal to ```L``` specified in the creation of the Master Secret Key.

* **Decrypting Messages**
``` cpp
std::string msgPrime = skG->decrypt( ctxt );
```
Decryption is only attempted if the ID beloning to ```ctxt``` evaluates to zero for the given gate.

* **Example**
``` cpp
// define N, Q, K, L
  
long* samples = gring::readOrBuildSamples<N,Q,K>( filename, smoothing_parameter );
  
// setup
auto msk = new gring::ABEMasterSecret<N,Q,K>{L, samples};
auto pk = msk->publicKey();
  
// creating IDs
std::vector< size_t > idWeights{1, 0, Q-1, 5, 3}
auto id = gring::ABEId{ idWeights };
  
// encrypting
std::string msg = "encrypt me";
auto ctxt = pk->encrypt( id, msg );
  
// creating gates
std::vector< size_t > addWeights{1, 0, 4, 5, 2};
auto gAdd = new gring::ABEAddGate<N,Q,K>{addWeights};
  
// decrypting
auto skG = msk->keyGen( gAdd );
std::string msgPrime = skG->decrypt( ctxt );
  
delete [] samples;
```







