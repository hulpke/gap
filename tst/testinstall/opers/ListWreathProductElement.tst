#
gap> START_TEST("ListWreathProductElement.tst");

#
# Perm Wreath Product In Imprimitive Action
#

#
gap> K := SymmetricGroup(3);;
gap> H := SymmetricGroup(5);;
gap> G := WreathProduct(K, H);;

#
gap> list := [(1,2), (1,2,3), (), (1,2,3), (), (1,2)(3,4)];;
gap> x := WreathProductElementList(G, list);;
gap> list = ListWreathProductElement(G, x);
true

#
gap> list := [(1,2), (1,2,3), (), (1,2,3), (), (1,2,3)];;
gap> x := WreathProductElementList(G, list);;
gap> list = ListWreathProductElement(G, x);
true

#
gap> x := (1,5,9,12,13)(2,4,8,11,15,3,6,7,10,14);;
gap> list := ListWreathProductElement(G, x);;
gap> x = WreathProductElementList(G, list);
true

# Top Component fails
gap> x := (1,5);;
gap> x ^ Projection(G);
fail
gap> ListWreathProductElement(G, x);
fail

# Base Component fails
gap> x := (1,2,3,4,5);;
gap> x ^ Projection(G);
()
gap> ListWreathProductElement(G, x);
fail

#
# Perm Wreath Product In Product Action
#

#
gap> K := SymmetricGroup(3);;
gap> H := SymmetricGroup(5);;
gap> G := WreathProductProductAction(K, H);;

#
gap> list := [(1,2), (1,2,3), (), (1,2,3), (), (1,2)(3,4)];;
gap> x := WreathProductElementList(G, list);;
gap> list = ListWreathProductElement(G, x);
true

#
gap> list := [(1,2), (1,2,3), (), (1,2,3), (), (1,2,3)];;
gap> x := WreathProductElementList(G, list);;
gap> list = ListWreathProductElement(G, x);
true

#
gap> x :=
> (  1, 25, 24, 12, 17,  5)(  2,  7, 22, 21, 18, 14)(  3, 16, 23)(  4, 19, 27, 15, 11,  8)(  6, 10, 26)(  9, 13, 20)
> ( 28,106, 78,174, 44, 86, 55,187, 51, 93, 71,167)( 29, 88, 76,183, 45, 95, 56,169, 49,102, 72,176)
> ( 30, 97, 77,165, 43,104, 57,178, 50, 84, 70,185)( 31,100, 81,177, 38, 89, 58,181, 54, 96, 65,170)
> ( 32, 82, 79,186, 39, 98, 59,163, 52,105, 66,179)( 33, 91, 80,168, 37,107, 60,172, 53, 87, 64,188)
> ( 34,103, 75,180, 41, 83, 61,184, 48, 99, 68,164)( 35, 85, 73,189, 42, 92, 62,166, 46,108, 69,173)
> ( 36, 94, 74,171, 40,101, 63,175, 47, 90, 67,182)(109,160,240,201,125,140,217,214,132,147,233,194)
> (110,142,238,210,126,149,218,196,130,156,234,203)(111,151,239,192,124,158,219,205,131,138,232,212)
> (112,154,243,204,119,143,220,208,135,150,227,197)(113,136,241,213,120,152,221,190,133,159,228,206)
> (114,145,242,195,118,161,222,199,134,141,226,215)(115,157,237,207,122,137,223,211,129,153,230,191)
> (116,139,235,216,123,146,224,193,127,162,231,200)(117,148,236,198,121,155,225,202,128,144,229,209);;
gap> list := ListWreathProductElement(G, x);;
gap> x = WreathProductElementList(G, list);
true

# Top Projection fails
gap> x := (4,16);;
gap> x ^ Projection(G);
fail
gap> ListWreathProductElement(G, x);
fail

# Base Decomposition fails
gap> x := (1,2,3,4,5);;
gap> x ^ Projection(G);
()
gap> ListWreathProductElement(G, x);
fail

#
# Matrix Wreath Product
#

#
gap> K := GL(3,3);;
gap> H := SymmetricGroup(5);;
gap> G := WreathProduct(K, H);;

#
gap> list := [K.1, K.2, One(K), K.2, One(K), (1,2)(3,4)];;
gap> x := WreathProductElementList(G, list);;
gap> list = ListWreathProductElement(G, x);
true

#
gap> list := [K.1, K.2, One(K), K.2, One(K), (1,2,3)];;
gap> x := WreathProductElementList(G, list);;
gap> list = ListWreathProductElement(G, x);
true

#
gap> x := [
> [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, Z(3)^0, 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3), Z(3)^0, Z(3)^0 ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ Z(3)^0, Z(3)^0, Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), Z(3)^0, Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, Z(3)^0, Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ] ];;
gap> list := ListWreathProductElement(G, x);;
gap> x = WreathProductElementList(G, list);
true
gap> ConvertToMatrixRep(x);;
gap> list := ListWreathProductElement(G, x);;
gap> x = WreathProductElementList(G, list);
true

# Top Projection fails
gap> x := [
> [ Z(3), 0*Z(3), Z(3)^0, Z(3), Z(3)^0, Z(3)^0, Z(3)^0, Z(3), Z(3), Z(3)^0, 0*Z(3), 0*Z(3), Z(3)^0, Z(3)^0, Z(3)^0 ],
> [ Z(3), Z(3), Z(3)^0, 0*Z(3), Z(3)^0, Z(3)^0, Z(3)^0, Z(3)^0, Z(3)^0, Z(3), Z(3)^0, 0*Z(3), Z(3)^0, Z(3), Z(3) ],
> [ 0*Z(3), Z(3)^0, 0*Z(3), Z(3)^0, Z(3)^0, Z(3)^0, Z(3), 0*Z(3), Z(3)^0, Z(3), Z(3)^0, Z(3)^0, 0*Z(3), 0*Z(3), Z(3)^0 ],
> [ Z(3), 0*Z(3), Z(3), 0*Z(3), Z(3), Z(3), Z(3)^0, Z(3), Z(3), Z(3), Z(3), Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ],
> [ Z(3)^0, Z(3), Z(3)^0, 0*Z(3), 0*Z(3), Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, Z(3), 0*Z(3), Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3), Z(3)^0, 0*Z(3), Z(3)^0, Z(3), 0*Z(3), Z(3), Z(3)^0, 0*Z(3) ],
> [ Z(3)^0, 0*Z(3), Z(3)^0, Z(3)^0, 0*Z(3), Z(3), Z(3)^0, Z(3)^0, 0*Z(3), Z(3), 0*Z(3), 0*Z(3), Z(3), Z(3), Z(3)^0 ],
> [ Z(3)^0, Z(3)^0, Z(3)^0, 0*Z(3), Z(3)^0, Z(3), Z(3)^0, 0*Z(3), 0*Z(3), Z(3), Z(3)^0, Z(3)^0, Z(3), Z(3), 0*Z(3) ],
> [ Z(3), Z(3)^0, Z(3)^0, Z(3), Z(3)^0, 0*Z(3), Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ],
> [ Z(3), 0*Z(3), Z(3)^0, Z(3)^0, Z(3)^0, Z(3)^0, 0*Z(3), 0*Z(3), Z(3), Z(3), Z(3), Z(3)^0, Z(3)^0, Z(3)^0, 0*Z(3) ],
> [ 0*Z(3), Z(3)^0, Z(3)^0, Z(3)^0, 0*Z(3), 0*Z(3), Z(3)^0, Z(3), Z(3), 0*Z(3), 0*Z(3), Z(3), Z(3), Z(3)^0, Z(3) ],
> [ 0*Z(3), Z(3)^0, 0*Z(3), Z(3), Z(3)^0, Z(3)^0, Z(3)^0, 0*Z(3), Z(3), Z(3)^0, Z(3)^0, 0*Z(3), Z(3)^0, Z(3)^0, Z(3) ],
> [ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), Z(3), 0*Z(3), Z(3), 0*Z(3), Z(3), Z(3), 0*Z(3), Z(3)^0, Z(3)^0 ],
> [ 0*Z(3), 0*Z(3), Z(3), Z(3)^0, Z(3)^0, Z(3)^0, 0*Z(3), Z(3)^0, Z(3), Z(3), 0*Z(3), Z(3)^0, Z(3)^0, Z(3)^0, Z(3) ],
> [ 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), Z(3), 0*Z(3), Z(3)^0, Z(3), Z(3)^0, 0*Z(3), Z(3), Z(3)^0 ] ];;
gap> x ^ Projection(G);
fail
gap> ListWreathProductElement(G, x);
fail
gap> ConvertToMatrixRep(x);;
gap> x ^ Projection(G);
fail
gap> ListWreathProductElement(G, x);
fail

# Base Decomposition fails
gap> x := [
> [ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3), Z(3)^0, Z(3)^0, Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3), Z(3), 0*Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3), Z(3), 0*Z(3), Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, Z(3)^0, Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), Z(3)^0 ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), Z(3) ],
> [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, Z(3)^0 ] ];;
gap> x ^ Projection(G);
()
gap> ListWreathProductElement(G, x);
fail
gap> ConvertToMatrixRep(x);;
gap> x ^ Projection(G);
()
gap> ListWreathProductElement(G, x);
fail

# Bugfix, Projection checked only diagonal entries
gap> x :=
> [ [  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  Z(3)^0,  Z(3)^0,  0*Z(3),  0*Z(3),  0*Z(3) ],
>   [  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  Z(3)^0,  0*Z(3),    Z(3),  0*Z(3),  0*Z(3),  0*Z(3) ],
>   [  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),    Z(3),  Z(3)^0,  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3) ],
>   [  0*Z(3),  Z(3)^0,  Z(3)^0,  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3) ],
>   [  Z(3)^0,  0*Z(3),    Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3) ],
>   [  0*Z(3),  0*Z(3),  Z(3)^0,  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3) ],
>   [  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  Z(3)^0,  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3) ],
>   [  0*Z(3),  0*Z(3),  0*Z(3),  Z(3)^0,  Z(3)^0,  Z(3)^0,  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3) ],
>   [  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  Z(3)^0,    Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3) ],
>   [  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),    Z(3),  0*Z(3),  0*Z(3) ],
>   [  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  Z(3)^0,  Z(3)^0,  0*Z(3) ],
>   [  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),    Z(3),  Z(3)^0 ],
>   [  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),    Z(3),  Z(3)^0,    Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3) ],
>   [  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),    Z(3),  Z(3)^0,  Z(3)^0,  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3) ],
>   [  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),    Z(3),  0*Z(3),  Z(3)^0,  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3),  0*Z(3) ] ];;
gap> x ^ Projection(G);
(1,4,5,3,2)
gap> ListWreathProductElement(G, x) <> fail;
true
gap> ConvertToMatrixRep(x);;
gap> x ^ Projection(G);
(1,4,5,3,2)
gap> ListWreathProductElement(G, x) <> fail;
true

#
# Generic Wreath Product
#

#
gap> K := DihedralGroup(12);;
gap> H := SymmetricGroup(5);;
gap> G := WreathProduct(K, H);;

#
gap> list := [K.1, K.2, One(K), K.2, One(K), (1,2)(3,4)];;
gap> x := WreathProductElementList(G, list);;
gap> list = ListWreathProductElement(G, x);
true

#
gap> x := G.1 * G.5 * G.2 * G.5 ^ ((G.4 * G.5) ^ 2) * G.2;;
gap> list := ListWreathProductElement(G, x);;
gap> x = WreathProductElementList(G, list);
true

#
# Generic Wreath Product : Bugfix for immutable lists
#

#
gap> K := FreeGroup("x", "y");;
gap> x := K.1;;
gap> y := K.2;;
gap> H := SymmetricGroup(3);;
gap> W := WreathProduct(K, H);;
gap> l := [x*y, x, y, (1,2,3)];;
gap> MakeImmutable(l);;
gap> w := WreathProductElementList(W, l);;
gap> l = [x*y, x, y, (1,2,3)];
true

#
gap> STOP_TEST("ListWreathProductElement.tst");
