library( MixtureFitting )

x = read.csv("inputs/enzyme.csv", header = TRUE)[,1]

# Resulting values taken from Lin et al. (2007), Statistica Sinica 17, 909-92, Table 1
p = c( 0.6240, 1 - 0.6240, 0.0949, 0.7802, 0.1331^2, 0.7150^2, 3.2780, 6.6684 )

sf = snmm_init_vector_kmeans( x, 2 )
for (i in 1:5000) {
    sf = snmm_fit_em(x, sf)
}

if( !all( abs( sf - p ) < 0.1 ) ) {
    stop( sf )
}
