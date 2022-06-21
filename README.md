## Description of the sage-code implementation

This sage code implements the SwiftEC and ElligatorSwift algorithms and contains scripts for generating the required precomputation parameters.

### Example runs

Compute and benchmark a single encoding-decoding test with curve P256
`sage enc_test.py -p P256`

Compute and benchmark a single projective X-only encoding-decoding test with curve P256
`sage Xenc_test.py -p P256`

For curves with a=0, you must call `enc_test_0.py` and `Xenc_test_0.py` instead
`sage enc_test_0.py -p secp256k1`
`sage Xenc_test_0.py -p secp256k1`

Compute and benchmark a single encoding-decoding test with a curve that requires the isogeny trick (never needed for a=0)
`sage encisog_test.py -p P25519`

Generate the precomputation parameters for a new curve (both for a=0 and a!=0)
`sage generate_parameters.py -p P256`
`sage generate_parameters.py -p secp256k1`

Test precomputation parameters (not implemented for a=0)
`sage param_test.py -p P256`

Script for exhaustively counting the number of encoding preimages (only run this for less-than-15-bit primes) (not implemented for a=0)
`sage preimage_test.py -p P15`

### Adding new curves

To add a new curve E: x^3 + a*x + b over a prime field F_p, give it a name and enter p,a,b (one per line in this order, in decimal form) to the file curves/PRIME_NAME

You must then run `sage generate_parameters -p PRIME_NAME`.

Example for curve E: x^3 + 2*x + 7 over F_13,
```bash
echo 13\n2\n7 > curves/MyCurve
sage generate_parameters -p MyCurve
sage enc_test -p MyCurve
```

Additionally, when using the isogeny trick, add the corresponding objects to isogenies.py

## Authors

1. **Jorge Chávez-Saab** <jchavez@computacion.cs.cinvestav.mx>,
2. **Francisco Rodríguez-Henríquez** <francisco@cs.cinvestav.mx>, 
3. **Mehdi Tibouchi** <mehdi.tibouchi@normalesup.org>, 
