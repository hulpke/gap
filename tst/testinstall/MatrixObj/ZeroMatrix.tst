gap> START_TEST("ZeroMatrix.tst");
gap> ReadGapRoot("tst/testinstall/MatrixObj/testmatobj.g");

#
# IsGF2MatrixRep
#
gap> TestZeroMatrix(IsGF2MatrixRep, GF(2), 2, 3);
<a 2x3 matrix over GF2>
gap> TestZeroMatrix(IsGF2MatrixRep, GF(2), 2, 0);
Error, Assertion failure
gap> TestZeroMatrix(IsGF2MatrixRep, GF(2), 0, 3); # TODO
Error, Assertion failure

# test error handling
gap> TestZeroMatrix(IsGF2MatrixRep, GF(2), -1, 3);
Error, Assertion failure
gap> TestZeroMatrix(IsGF2MatrixRep, GF(2), 2, -1);
Error, ZERO_GF2VEC_2: <len> must be a non-negative small integer (not the inte\
ger -1)

# test error handling
gap> TestZeroMatrix(IsGF2MatrixRep, GF(3), 2, 3);
Error, Assertion failure

#
# Is8BitMatrixRep
#
gap> TestZeroMatrix(Is8BitMatrixRep, GF(3), 2, 3);
[ [ 0*Z(3), 0*Z(3), 0*Z(3) ], [ 0*Z(3), 0*Z(3), 0*Z(3) ] ]
gap> TestZeroMatrix(Is8BitMatrixRep, GF(3), 2, 0);
Error, Assertion failure
gap> TestZeroMatrix(Is8BitMatrixRep, GF(3), 0, 3);
Error, Assertion failure

#
gap> TestZeroMatrix(Is8BitMatrixRep, GF(251), 2, 3);
[ [ 0*Z(251), 0*Z(251), 0*Z(251) ], [ 0*Z(251), 0*Z(251), 0*Z(251) ] ]
gap> TestZeroMatrix(Is8BitMatrixRep, GF(251), 2, 0);
Error, Assertion failure
gap> TestZeroMatrix(Is8BitMatrixRep, GF(251), 0, 3);
Error, Assertion failure

# test error handling
gap> TestZeroMatrix(Is8BitMatrixRep, GF(2), 2, 3);
Error, Assertion failure
gap> TestZeroMatrix(Is8BitMatrixRep, GF(257), 2, 3);
Error, Assertion failure

#
# IsPlistMatrixRep
#
gap> TestZeroMatrix(IsPlistMatrixRep, GF(2), 2, 3);
<2x3-matrix over GF(2)>
gap> TestZeroMatrix(IsPlistMatrixRep, GF(2), 2, 0);
<2x0-matrix over GF(2)>
gap> TestZeroMatrix(IsPlistMatrixRep, GF(2), 0, 3);
<0x3-matrix over GF(2)>

#
gap> TestZeroMatrix(IsPlistMatrixRep, Integers, 2, 3);
<2x3-matrix over Integers>
gap> TestZeroMatrix(IsPlistMatrixRep, Integers, 2, 0);
<2x0-matrix over Integers>
gap> TestZeroMatrix(IsPlistMatrixRep, Integers, 0, 3);
<0x3-matrix over Integers>

#
gap> TestZeroMatrix(IsPlistMatrixRep, Rationals, 2, 3);
<2x3-matrix over Rationals>
gap> TestZeroMatrix(IsPlistMatrixRep, Rationals, 2, 0);
<2x0-matrix over Rationals>
gap> TestZeroMatrix(IsPlistMatrixRep, Rationals, 0, 3);
<0x3-matrix over Rationals>

#
gap> TestZeroMatrix(IsPlistMatrixRep, Integers mod 4, 2, 3);
<2x3-matrix over (Integers mod 4)>
gap> TestZeroMatrix(IsPlistMatrixRep, Integers mod 4, 2, 0);
<2x0-matrix over (Integers mod 4)>
gap> TestZeroMatrix(IsPlistMatrixRep, Integers mod 4, 0, 3);
<0x3-matrix over (Integers mod 4)>

#
# Test ZeroMatrix variant which "guesses" a suitable representation, i.e.:
#    ZeroMatrix( <R>, <m>, <n> )
#

#
gap> ZeroMatrix(Integers, 2, 3);
<2x3-matrix over Integers>
gap> ZeroMatrix(Integers, 0, 3);
<0x3-matrix over Integers>
gap> ZeroMatrix(Integers, 2, 0);
<2x0-matrix over Integers>

#
gap> ZeroMatrix(Integers mod 4, 2, 3);
<2x3-matrix over (Integers mod 4)>
gap> ZeroMatrix(Integers mod 4, 0, 3);
<0x3-matrix over (Integers mod 4)>
gap> ZeroMatrix(Integers mod 4, 2, 0);
<2x0-matrix over (Integers mod 4)>

#
gap> ZeroMatrix(GF(2), 2, 3);
<a 2x3 matrix over GF2>
gap> ZeroMatrix(GF(2), 0, 3);
<a 1x3 matrix over GF2>
gap> ZeroMatrix(GF(2), 2, 0);
[ <a GF2 vector of length 0>, <a GF2 vector of length 0> ]

#
gap> ZeroMatrix(GF(3), 2, 3);
[ [ 0*Z(3), 0*Z(3), 0*Z(3) ], [ 0*Z(3), 0*Z(3), 0*Z(3) ] ]
gap> ZeroMatrix(GF(3), 0, 3);
[ [ 0*Z(3), 0*Z(3), 0*Z(3) ] ]
gap> ZeroMatrix(GF(3), 2, 0);
[ [  ], [  ] ]

#
gap> ZeroMatrix(GF(4), 2, 3);
[ [ 0*Z(2), 0*Z(2), 0*Z(2) ], [ 0*Z(2), 0*Z(2), 0*Z(2) ] ]
gap> ZeroMatrix(GF(4), 0, 3);
[ [ 0*Z(2), 0*Z(2), 0*Z(2) ] ]
gap> ZeroMatrix(GF(4), 2, 0);
[ [  ], [  ] ]

#
gap> STOP_TEST("ZeroMatrix.tst");