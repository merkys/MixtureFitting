load( "inputs/2-2-cryptand-dihedrals-trimmed.RData" )
load( "outputs/2-2-cryptand-dihedrals.RData" )
library( MixtureFitting )
init = list()
parameters = list()
ll = list()
for( i in 1:10 ) {
    vinit = vmm_init_vector( i )
    vf = vmm_fit_em( dihedrals, vinit )
    init[[i]] = vinit
    parameters[[i]] = vf$p
    ll[[i]] = llvmm( dihedrals, vf$p )
}
if( all.equal( init, output_init ) != TRUE ) {
    stop( all.equal( init, output_init ) )
}
if( all.equal( parameters, output_parameters ) != TRUE ) {
    stop( all.equal( parameters, output_parameters ) )
}
if( all.equal( ll, output_ll ) != TRUE ) {
    stop( all.equal( ll, output_ll ) )
}
