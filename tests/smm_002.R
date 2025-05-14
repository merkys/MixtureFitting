library( MixtureFitting )
load( "outputs/smm_002.RData" )
set.seed( 42 )
x = rt( 100, 10 )
sinit = smm_init_vector( x, 3 )
sfjt = smm_fit_em( x, sinit, polyroot.solution = "jenkins_taub" )
sfnr = smm_fit_em( x, sinit, polyroot.solution = "newton_raphson" )
if( all.equal( sinit, output_init ) != TRUE ) {
    stop( all.equal( sinit, output_init ) )
}
if( all.equal( sfjt$p, output_parameters_jt ) != TRUE ) {
    stop( all.equal( sfjt$p, output_parameters_jt ) )
}
if( all.equal( sfjt$p, output_parameters_nr ) != TRUE ) {
    stop( all.equal( sfjt$p, output_parameters_nr ) )
}
