library( MixtureFitting )

set.seed(42)
p = c( 0.5, 0.5, 1.5, 10, 3, 1 )

x = c( rnorm(2000 * p[1], p[3], p[5]), rnorm(2000 * p[2], p[4], p[6]) )

init = gmm_init_vector_quantile( x, 2 )
if( !all( abs( init[3:4] - c( 2.7, 9.5 ) ) < 0.1 ) ) {
    stop( 1,   init[3:4] - c( 2.7, 9.5 ) )
}

init = gmm_init_vector_quantile( x, 2, c( numeric(1000)+1, numeric(1000) ) )
if( !all( abs( init[3:4] - c( 0.2, 2.7 ) ) < 0.1 ) ) {
    stop( 2,   init[3:4] - c( 0.2, 2.7 ) )
}

init = gmm_init_vector_quantile( x, 2, c( numeric(1000), numeric(1000)+1 ) )
if( !all( abs( init[3:4] - c( 9.5, 10.4 ) ) < 0.1 ) ) {
    stop( 2,   init[3:4] - c( 9.5, 10.4 ) )
}
