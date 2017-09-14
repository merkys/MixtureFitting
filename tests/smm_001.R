library( MixtureFitting )
load( "outputs/smm_001.RData" )
set.seed( 42 )
x = rt( 100, 10 )
sinit = smm_init_vector( x, 1 )
sf = smm_fit_em_CWL04( x, sinit, debug = TRUE )
if( all.equal( sinit, output_init ) != TRUE ) {
    stop( all.equal( sinit, output_init ) )
}
if( all.equal( sf$p, output_parameters ) != TRUE ) {
    stop( all.equal( sf$p, output_parameters ) )
}
