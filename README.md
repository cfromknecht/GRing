GRing
=====
GRing is a lightweight, efficient implementation of G-Trapdoors on polynomials under the Ring-LWE assumption.  It includes operations for arithmetic, memory management, and trapdoor operations.  The current demo shows how one can use GRing to accomplish a simple IBE scheme.
  
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
The current output should a blank message, so every entry in the printed polynomial should be small.

Parameters
----------
Values for N, Q, K, etc. can be found at the top of the main.cpp file.  Note for simplicity, we use K to derive Q to ensure that Q is a power of 2.

Notes
-----
The encoder from hashes to polynomials can only support ```K <= 16```.  Behavior is undefined for ```K > 16```.  
Future improvements to the API will help to increase readibility.
