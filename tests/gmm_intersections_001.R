library( MixtureFitting )
P = matrix( c( 0.0443473274465592, 0.130138383023397,  0.825514289530041,
               2.21235358083254,   2.00352959127296,   2.03247868665627,
               0.27455506306116,   0.0275176004414257, 0.106694771531626 ),
            ncol = 3 )

test1 = all.equal( gmm_intersections( P[c(1,3),] ),
                   c( 1.66800776372449, 2.33295709358655 ) )
test2 = all.equal( gmm_intersections( P[c(2,3),] ),
                   c() )
test3 = all.equal( gmm_intersections( P[c(1,1),] ),
                   NaN )

if( test1 != TRUE ) {
    stop( test1 )
}
if( test2 != TRUE ) {
    stop( test2 )
}
if( test3 != TRUE ) {
    stop( test3 )
}
