library( datasets )
library( MixtureFitting )

x = faithful[,1]
# Resulting values taken from Lin et al. (2007), Statistica Sinica 17, 909-92, Table 3
p = c( 0.3487, 1 - 0.3487, 1.7267, 4.8026, 0.3801, 0.6857, 5.8026, -3.4951 )

sf = snmm_fit_em( x, snmm_init_vector_kmeans( x, 2 ) )

if( !all( abs( sf - p ) < 0.1 ) ) {
    stop( sf )
}
