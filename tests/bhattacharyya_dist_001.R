library( MixtureFitting )
orig1 = 0.125
dist1 = bhattacharyya_dist( 1, 2, 1, 1 )
orig2 = 1.95606150021424
dist2 = bhattacharyya_dist( 0, 0, 1, 0.01 )
if( all.equal( c( orig1, orig2 ), c( dist1, dist2 ) ) != TRUE ) {
    stop( all.equal( c( orig1, orig2 ), c( dist1, dist2 ) )
}
