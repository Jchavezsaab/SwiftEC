## Description of the sage-code implementation

This sage code implements the SwiftEC and ElligatorSwift algorithms and contains scripts for generating the required precomputation parameters.

### Example runs

Compute and benchmark a single encoding-decoding test with curve P256
`sage enc_test.py -p P256`

Compute and benchmark a single projective X-only encoding-decoding test with curve P256
`sage Xenc_test.py -p P256`

Compute and benchmark a single encoding-decoding test with a curve that requires the isogeny trick
`sage encisog_test.py -p P25519`

Generate the precomputation parameters for a new curve
`sage generate_parameters.py -p P256`

Test precomputation parameters
`sage param_test.py -p P256`

Script for exhaustively counting the number of encoding preimages (only run this for less-than-15-bit primes)
`sage preimage_test.py -p P15`

### Adding new curves

To add a new curve E: x^3 + A*x + B over a prime field F_p, give it a name and enter p,A,B (one per line in this order, in decimal form) to the file curves/PRIME_NAME

You must then run `sage generate_parameters -p PRIME_NAME`.

Example for curve E: x^3 + 2*x + 7 over F_13,
```bash
echo 13\n2\n7 > curves/MyCurve
sage generate_parameters -p MyCurve
sage enc_test -p MyCurve
```

## Authors

1. **Jorge Chávez-Saab** <jchavez@computacion.cs.cinvestav.mx>,
2. **Francisco Rodríguez-Henríquez** <francisco@cs.cinvestav.mx>, 
3. **Mehdi Tibouchi** <mehdi.tibouchi@normalesup.org>, 
