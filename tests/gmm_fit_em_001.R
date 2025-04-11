library( MixtureFitting )

set.seed(42)
p = c( 0.5, 0.5, 1.5, 10, 3, 1 )

x = c( rnorm(2000 * p[1], p[3], p[5]), rnorm(2000 * p[2], p[4], p[6]) )
init = gmm_init_vector( x, 2 )

gf = gmm_fit_em( x, init, implementation = "C" )
if( !all( abs( gf$p - c( 0.5, 0.5, 1.5, 10, 3, 1 ) ) < 0.1 ) ) {
    stop( 1, gf$p - p )
}

gf = gmm_fit_em( x, init, implementation = "R" )
if( !all( abs( gf$p - c( 0.5, 0.5, 1.5, 10, 3, 1 ) ) < 0.1 ) ) {
    stop( 2, gf$p - p )
}

gf = gmm_fit_em( x, init, x * 0 + 0.5, implementation = "C" )
if( !all( abs( gf$p - c( 0.5, 0.5, 1.5, 10, 3, 1 ) ) < 0.1 ) ) {
    stop( 3, gf$p - p )
}

gf = gmm_fit_em( x, init, x * 0 + 0.5, implementation = "R" )
if( !all( abs( gf$p - c( 0.5, 0.5, 1.5, 10, 3, 1 ) ) < 0.1 ) ) {
    stop( 4, gf$p - p )
}
