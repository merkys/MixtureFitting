dgmm <- function( x, p )
{
    if( length( p[is.na(p)] ) > 0 ) {
        return( rep( NaN, times = length( x ) ) )
    }
    buffer = numeric( length(x) )
    ret = .C( "dgmm",
              as.double(x),
              as.integer( length(x) ),
              as.double(p),
              as.integer( length(p) ),
              retvec = numeric( length(x) ))$retvec
    buffer[1:length(x)] <- ret[1:length(x)]
    return( buffer )
}

dvmm <- function( x, p )
{
    if( length( p[is.na(p)] ) > 0 ) {
        return( rep( NaN, times = length( x ) ) )
    }
    buffer = numeric( length(x) )
    ret = .C( "dvmm",
              as.double(x),
              as.integer( length(x) ),
              as.double(p),
              as.integer( length(p) ),
              retvec = numeric( length(x) ))$retvec
    buffer[1:length(x)] <- ret[1:length(x)]
    return( buffer )
}

dcmm <- function( x, p )
{
    if( length( p[is.na(p)] ) > 0 ) {
        return( rep( NaN, times = length( x ) ) )
    }
    buffer = numeric( length(x) )
    ret = .C( "dcmm",
              as.double(x),
              as.integer( length(x) ),
              as.double(p),
              as.integer( length(p) ),
              retvec = numeric( length(x) ))$retvec
    buffer[1:length(x)] <- ret[1:length(x)]
    return( buffer )
}

llgmm <- function( x, p )
{
    if( length( p[is.na(p)] ) > 0 ) {
        return( NaN )
    }
    ret = .C( "llgmm",
              as.double(x),
              as.integer( length(x) ),
              as.double(p),
              as.integer( length(p) ),
              retvec = numeric(1) )$retvec
    return( ret )
}

llvmm <- function( x, p )
{
    if( length( p[is.na(p)] ) > 0 ) {
        return( NaN )
    }
    ret = .C( "llvmm",
              as.double(x),
              as.integer( length(x) ),
              as.double(p),
              as.integer( length(p) ),
              retvec = numeric(1) )$retvec
    return( ret )
}

llcmm <- function( x, p )
{
    if( length( p[is.na(p)] ) > 0 ) {
        return( NaN )
    }
    ret = .C( "llcmm",
              as.double(x),
              as.integer( length(x) ),
              as.double(p),
              as.integer( length(p) ),
              retvec = numeric(1) )$retvec
    return( ret )
}

gmm_fit_em <- function( x, p, epsilon = c( 0.000001, 0.000001, 0.000001 ),
                        debug = FALSE )
{
    ret = .C( "gmm_fit_em",
              as.double(x),
              as.integer( length(x) ),
              as.double(p),
              as.integer( length(p) ),
              as.double( epsilon ),
              as.integer( debug ),
              retvec = numeric( length(p) ),
              steps = integer(1) )
    l = list( p = ret$retvec, steps = ret$steps )
    return( l )
}

vmm_fit_em <- function( x, p,
                        epsilon = c( 0.000001, 0.000001, 0.000001 ),
                        debug = FALSE )
{
    l = vmm_fit_em_by_diff( x, p, epsilon, debug )
    return( l )
}

vmm_fit_em_by_diff <- function( x, p,
                                epsilon = c( 0.000001, 0.000001, 0.000001 ),
                                debug = FALSE )
{
    debugflag = 0
    if( debug == TRUE ) {
        debugflag = 1
    }
    ret = .C( "vmm_fit_em_by_diff",
              as.double(x),
              as.integer( length(x) ),
              as.double(p),
              as.integer( length(p) ),
              as.double( epsilon ),
              as.integer( debug ),
              retvec = numeric( length(p) ),
              steps = integer(1) )
    l = list( p = ret$retvec, steps = ret$steps )
    return( l )
}

vmm_fit_em_by_ll <- function( x, p,
                              epsilon = .Machine$double.eps,
                              debug = FALSE )
{
    debugflag = 0
    if( debug == TRUE ) {
        debugflag = 1
    }
    ret = .C( "vmm_fit_em_by_ll",
              as.double(x),
              as.integer( length(x) ),
              as.double(p),
              as.integer( length(p) ),
              as.double( epsilon ),
              as.integer( debug ),
              retvec = numeric( length(p) ),
              steps = integer(1) )
    l = list( p = ret$retvec, steps = ret$steps )
    return( l )
}

cmm_fit_em <- function( x, p, epsilon = c( 0.000001, 0.000001, 0.000001 ),
                        iter.cauchy = 20, debug = FALSE )
{
    debugflag = 0
    if( debug == TRUE ) {
        debugflag = 1
    }
    ret = .C( "cmm_fit_em",
              as.double(x),
              as.integer( length(x) ),
              as.double(p),
              as.integer( length(p) ),
              as.double( epsilon ),
              as.integer( iter.cauchy ),
              as.integer( debug ),
              retvec = numeric( length(p) ),
              steps = integer(1) )
    l = list( p = ret$retvec, steps = ret$steps )
    return( l )
}

gmm_init_vector <- function( x, m )
{
    ret = .C( "gmm_init_vector",
              as.double(x),
              as.integer( length(x) ),
              as.integer(m),
              retvec = numeric( 3*m ) )
    return( ret$retvec )
}

cmm_init_vector <- function( x, m )
{
    ret = .C( "cmm_init_vector",
              as.double(x),
              as.integer( length(x) ),
              as.integer(m),
              retvec = numeric( 3*m ) )
    return( ret$retvec )
}

vmm_init_vector <- function( m )
{
    ret = .C( "vmm_init_vector",
              as.integer(m),
              retvec = numeric( 3*m ) )
    return( ret$retvec )
}

polyroot_NR <- function( p, init = 0, epsilon = 1e-6, debug = FALSE )
{
    ret = .C( "polyroot_NR",
              as.double(p),
              as.integer( length(p) ),
              as.double(init),
              as.double(epsilon),
              as.integer( debug ),
              retvec = numeric(1) )
    return( ret$retvec )
}

#=========================================================================
# R counterparts of functions, rewritten in C
#=========================================================================

# Returns distribution for a Gaussian Mixture Model at point x
# p - vector of 3*m parameters, where m is size of mixture in peaks.
# p[1:m] -- mixture proportions
# p[(m+1):(2*m)] -- means of the peaks
# p[(2*m+1):(3*m)] -- dispersions of the peaks
dgmm_R <- function( x, p, normalise_proportions = FALSE,
                    restrict_sigmas = FALSE )
{
    m     = length(p)/3
    A     = p[1:m]
    mu    = p[(m+1):(2*m)]
    sigma = p[(2*m+1):(3*m)]
    if( normalise_proportions == TRUE ) {
        A = A/sum(A)
    }
    sum = numeric( length(x) )
    for( i in 1:m ) {
        if( sigma[i] > 0 | restrict_sigmas == FALSE ) {
            sum = sum + A[i] * dnorm( x, mu[i], sigma[i] )
        }
    }
    return( sum )
}

dvmm_R <- function( x, p )
{
    m  = length(p)/3
    A  = p[1:m]
    mu = p[(m+1):(2*m)]
    k  = p[(2*m+1):(3*m)]
    sum = 0
    for( i in 1:m ) {
        sum = sum + A[i] * exp( k[i] * cos( deg2rad( x - mu[i] ) ) ) /
                           ( 2 * pi * besselI( k[i], 0 ) )
    }
    return( sum )
}

dcmm_R <- function (x, p)
{
    m = length(p)/3
    A = p[1:m]
    c = p[(m+1):(2*m)]
    s = p[(2*m+1):(3*m)]
    sum = numeric( length( x ) ) * 0
    for( i in 1:m ) {
        sum = sum + A[i] * dcauchy( x, c[i], s[i] )
    }
    return( sum )
}

# Log-likelihood of data x given Gaussian Mixture Model, described by
# parameter vector p.
llgmm_R <- function (x, p)
{
    n     = length(x)
    m     = length(p)/3
    A     = p[1:m]
    mu    = p[(m+1):(2*m)]
    sigma = p[(2*m+1):(3*m)]
    if( length(p[is.na(p)]) > 0 ) {
        return( NaN )
    }
    if( length(A[A>=0]) < length(A) ||
        length(sigma[sigma>=0]) < length(sigma) ) {
        return( -Inf )
    }
    diff = matrix(data = 0, nrow = n, ncol = m )
    ref_peak = vector( "numeric", n )
    for (i in 1:n) {
        diff[i,1:m] = ( x[i] - mu[1:m] ) ^ 2
        ref_peak[i] = which.min( diff[i,1:m] )
    }
    sum = sum( log( A[ref_peak] / (sqrt(2*pi) * sigma[ref_peak]) ) )
    for (i in 1:n) {
        rp = ref_peak[i]
        sum = sum - diff[i,rp] / ( 2 * sigma[rp]^2 )
        expsum = sum( exp( log( (A[1:m]*sigma[rp])/(A[rp]*sigma[1:m]) )
                            - diff[i,1:m] / (2*sigma[1:m]^2)
                            + diff[i,rp] / (2*sigma[rp]^2) ) )
        sum = sum + log( expsum )
    }
    return( sum )
}

llvmm_R <- function( x, p )
{
    n  = length(x)
    m  = length(p)/3
    A  = p[1:m]/sum(p[1:m])
    mu = p[(m+1):(2*m)]
    k  = p[(2*m+1):(3*m)]
    y  = vector( "numeric", n ) * 0
    for (i in 1:m) {
        y = y + dvmm( x, c( A[i], mu[i], k[i] ) )
    }
    return( sum( log( y ) ) )
}

llcmm_R <- function( x, p )
{
    n = length(x)
    m = length(p)/3
    A = p[1:m]/sum(p[1:m])
    c = p[(m+1):(2*m)]
    s = p[(2*m+1):(3*m)]
    y = numeric( n ) * 0
    for (i in 1:m) {
        y = y + dcmm( x, c( A[i], c[i], s[i] ) )
    }
    return( sum( log( y ) ) )
}

gmm_fit_em_R <- function( x, p, epsilon = c( 0.000001, 0.000001, 0.000001 ),
                          collect.history = FALSE, unif.component = FALSE,
                          convergence = abs_convergence )
{
    m     = length(p)/3
    A     = p[1:m]
    mu    = p[(m+1):(2*m)]
    sigma = p[(2*m+1):(3*m)]
    prev_A     = rep( Inf, m )
    prev_mu    = rep( Inf, m )
    prev_sigma = rep( Inf, m )
    steps = 0
    history = list()
    if( collect.history == TRUE ) {
        history[[1]] = p
    }
    while( steps == 0 ||
           !convergence( c( A, mu, sigma ), 
                         c( prev_A, prev_mu, prev_sigma ), epsilon ) ) {
        prev_A     = A
        prev_mu    = mu
        prev_sigma = sigma
        q = vector( "numeric", length( x ) ) * 0
        for( j in 1:m ) {
            q = q + A[j] * dnorm( x, mu[j], sigma[j] )
        }
        if( unif.component ) {
            # Allows additional component with uniform distribution for
            # the modelling of outliers as suggested in:
            # Cousineau, D. & Chartier, S.
            # Outliers detection and treatment: a review
            # International Journal of Psychological Research,
            # 2010, 3, 58-67
            # http://revistas.usb.edu.co/index.php/IJPR/article/view/844

            q = q + ( 1 - sum( A ) ) * dunif( x, min(x), max(x) )
        }
        for( j in 1:m ) {
            h = A[j] * dnorm( x, mu[j], sigma[j] ) / q
            A[j]     = sum( h ) / length( x )
            mu[j]    = sum( h * x ) / sum( h )
            sigma[j] = sqrt( sum( h * ( x - mu[j] ) ^ 2 ) / sum( h ) )
        }
        steps   = steps + 1
        if( collect.history == TRUE ) {
            history[[steps+1]] = c( A, mu, sigma )
        }
        if( length( A[    is.na(A)] )   +
            length( mu[   is.na(mu) ] ) +
            length( sigma[is.na(sigma)] ) > 0 ) {
            break
        }
    }
    p[1:m] = A
    p[(m+1):(2*m)] = mu
    p[(2*m+1):(3*m)] = sigma
    l = list( p = p, steps = steps )
    if( collect.history == TRUE ) {
        l$history = history
    }
    return( l )
}

# Estimate von Mises Mixture parameters using expectation maximisation.
# Implemented according to Banerjee et al., Expectation Maximization for
# Clustering on Hyperspheres, manuscript, 2003.
vmm_fit_em_R <- function( x, p, epsilon = c( 0.000001, 0.000001, 0.000001 ),
                          debug = FALSE )
{
    l = vmm_fit_em_by_diff_R( x, p, epsilon, debug )
    return( l )
}

vmm_fit_em_by_diff_R <- function( x, p,
                                  epsilon = c( 0.000001, 0.000001, 0.000001 ),
                                  debug = FALSE )
{
    m  = length(p)/3
    A  = p[1:m]
    mu = p[(m+1):(2*m)]
    k  = p[(2*m+1):(3*m)]
    d_A  = c( Inf )
    d_mu = c( Inf )
    d_k  = c( Inf )
    steps = 0
    while( length( d_A[ d_A  > epsilon[1]] ) > 0 ||
           length( d_mu[d_mu > epsilon[2]] ) > 0 ||
           length( d_k[ d_k  > epsilon[3]] ) > 0 ) {
        prev_A  = A
        prev_mu = mu
        prev_k  = k
        q = vector( "numeric", length( x ) ) * 0
        for( j in 1:m ) {
            q = q + dvmm( x, c( A[j], mu[j], k[j] ) )
        }
        for( j in 1:m ) {
            h = dvmm( x, c( A[j], mu[j], k[j] ) ) / q
            A[j]  = sum( h ) / length( x )
            mu[j] = rad2deg( atan2( sum( sin( deg2rad(x) ) * h ),
                                    sum( cos( deg2rad(x) ) * h ) ) )
            Rbar  = sqrt( sum( sin( deg2rad(x) ) * h )^2 +
                          sum( cos( deg2rad(x) ) * h )^2 ) / sum( h )
            k[j]  = ( 2 * Rbar - Rbar ^ 3 ) / ( 1 - Rbar ^ 2 )
            if( debug == TRUE ) {
                cat( A[j], " ", mu[j], " ", k[j], " " )
            }
        }
        if( debug == TRUE ) {
            cat( "\n" )
        }
        d_A  = abs( A  - prev_A  )
        d_mu = abs( mu - prev_mu )
        d_k  = abs( k  - prev_k  )
        steps = steps + 1
        if( length( d_A[ is.na(d_A)] )  +
            length( d_mu[is.na(d_mu)] ) +
            length( d_k[ is.na(d_k)] )  > 0 ) {
            break
        }
    }
    p[1:m] = A
    p[(m+1):(2*m)] = mu
    p[(2*m+1):(3*m)] = k
    return( list( p = p, steps = steps ) )
}

vmm_fit_em_by_ll_R <- function( x, p, epsilon = .Machine$double.eps,
                                debug = FALSE )
{
    m  = length(p)/3
    A  = p[1:m]
    mu = p[(m+1):(2*m)]
    k  = p[(2*m+1):(3*m)]
    prev_llog = llvmm( x, p )
    d_llog = Inf
    steps = 0
    while( !is.na(d_llog) && d_llog > epsilon ) {
        q = vector( "numeric", length( x ) ) * 0
        for( j in 1:m ) {
            q = q + dvmm( x, c( A[j], mu[j], k[j] ) )
        }
        for( j in 1:m ) {
            h = dvmm( x, c( A[j], mu[j], k[j] ) ) / q
            A[j]  = sum( h ) / length( x )
            mu[j] = rad2deg( atan2( sum( sin( deg2rad(x) ) * h ),
                                    sum( cos( deg2rad(x) ) * h ) ) )
            Rbar  = sqrt( sum( sin( deg2rad(x) ) * h )^2 +
                          sum( cos( deg2rad(x) ) * h )^2 ) / sum( h )
            k[j]  = ( 2 * Rbar - Rbar ^ 3 ) / ( 1 - Rbar ^ 2 )
            if( debug == TRUE ) {
                cat( A[j], " ", mu[j], " ", k[j], " " )
            }
        }
        if( debug == TRUE ) {
            cat( "\n" )
        }
        llog = llvmm( x, c( A, mu, k ) )
        d_llog = abs( llog - prev_llog )
        prev_llog = llog
        steps = steps + 1
        if( length( A[ is.na(A) ] ) +
            length( mu[is.na(mu)] ) +
            length( k[ is.na(k) ] ) > 0 ) {
            break
        }
    }
    p[1:m] = A
    p[(m+1):(2*m)] = mu
    p[(2*m+1):(3*m)] = k
    return( list( p = p, steps = steps ) )
}

# Estimate Cauchy Mixture parameters using expectation maximisation.
# Estimation of individual component's parameters is implemented according
# to Ferenc Nahy, Parameter Estimation of the Cauchy Distribution in
# Information Theory Approach, Journal of Universal Computer Science, 2006.
cmm_fit_em_R <- function( x, p, epsilon = c( 0.000001, 0.000001, 0.000001 ),
                          collect.history = FALSE,
                          unif.component = FALSE,
                          convergence = abs_convergence )
{
    m = length(p)/3
    A = p[1:m]
    c = p[(m+1):(2*m)]
    s = p[(2*m+1):(3*m)]
    prev_A = rep( Inf, m )
    prev_c = rep( Inf, m )
    prev_s = rep( Inf, m )
    steps = 0
    history = list()
    if( collect.history == TRUE ) {
        history[[1]] = p
    }
    while( steps == 0 ||
           !convergence( c( A, c, s ), 
                         c( prev_A, prev_c, prev_s ), epsilon ) ) {
        prev_A = A
        prev_c = c
        prev_s = s
        q = numeric( length( x ) ) * 0
        for( j in 1:m ) {
            q = q + A[j] * dcauchy( x, c[j], s[j] )
        }
        if( unif.component ) {
            # Allows additional component with uniform distribution for
            # the modelling of outliers as suggested in:
            # Cousineau, D. & Chartier, S.
            # Outliers detection and treatment: a review
            # International Journal of Psychological Research,
            # 2010, 3, 58-67
            # http://revistas.usb.edu.co/index.php/IJPR/article/view/844

            q = q + ( 1 - sum( A ) ) * dunif( x, min(x), max(x) )
        }
        for( j in 1:m ) {
            h = A[j] * dcauchy( x, c[j], s[j] ) / q
            A[j] = sum( h ) / length( x )
            cauchy_steps = 0
            prev_cj = Inf
            prev_sj = Inf
            while( cauchy_steps == 0 ||
                   !convergence( c( c[j], s[j] ),
                                 c( prev_cj, prev_sj ),
                                 epsilon[2:3] ) ) {
                prev_cj = c[j]
                prev_sj = s[j]
                e0k  = sum( h / (1+((x-c[j])/s[j])^2) ) / sum( h )
                e1k  = sum( h * ((x-c[j])/s[j])/(1+((x-c[j])/s[j])^2) ) /
                       sum( h )
                c[j] = c[j] + s[j] * e1k / e0k
                s[j] = s[j] * sqrt( 1/e0k - 1 )
                cauchy_steps = cauchy_steps + 1
            }
        }
        steps = steps + 1
        if( collect.history == TRUE ) {
            history[[steps+1]] = c( A, c, s )
        }
        if( length( A[is.na(A)] ) +
            length( c[is.na(c)] ) +
            length( s[is.na(s)] ) > 0 ) {
            break
        }
    }
    p[1:m] = A
    p[(m+1):(2*m)] = c
    p[(2*m+1):(3*m)] = s
    l = list( p = p, steps = steps )
    if( collect.history == TRUE ) {
        l$history = history
    }
    return( l )
}

gmm_init_vector_R <- function( x, m ) {
    start = numeric( 3 * m )
    start[1:m]           = 1/m
    start[(m+1):(2*m)]   = min(x) + (1:m)*(max(x)-min(x))/(m+1)
    start[(2*m+1):(3*m)] = (max(x)-min(x))/(m+1)/6
    return( start )
}

vmm_init_vector_R <- function( m ) {
    start = numeric( 3 * m )
    start[1:m]           = 1/m
    start[(m+1):(2*m)]   = 360 / m * seq( 0, m-1, 1 )
    start[(2*m+1):(3*m)] = (m/(12*180))^2
    return( start )
}

cmm_init_vector_R <- function( x, m ) {
    start = numeric( 3 * m )
    start[1:m]           = 1/m
    start[(m+1):(2*m)]   = min(x) + (1:m)*(max(x)-min(x))/(m+1)
    start[(2*m+1):(3*m)] = 1
    return( start )
}

# Finds one real polynomial root using Newton--Raphson method, implemented
# according to Wikipedia:
# https://en.wikipedia.org/w/index.php?title=Newton%27s_method&oldid=710342140
polyroot_NR_R <- function( p, init = 0, epsilon = 1e-6, debug = FALSE )
{
    x = init
    x_prev = Inf
    steps = 0

    n = length(p)

    d = p[2:n] * (1:(n-1))
    while( abs( x - x_prev ) > epsilon ) {
        x_prev = x
        powers = x^(0:(n-1))
        x = x - sum(p * powers) / sum(d * powers[1:(n-1)])
        steps = steps + 1
    }

    if( debug ) {
        cat( "Convergence reached after", steps, "iteration(s)\n" )
    }

    return( x )
}

#=========================================================================
# Functions, that are not yet rewritten in C
#=========================================================================

ds <- function( x, c, s, ni )
{
    return( dt( ( x - c ) / s, ni ) / s )
}

dsmm <- function( x, p )
{
    m = length( p ) / 4
    A  = p[1:m]
    c  = p[(m+1):(2*m)]
    s  = p[(2*m+1):(3*m)]
    ni = p[(3*m+1):(4*m)]

    ret = numeric( length( x ) )
    for( i in 1:m ) {
        ret = ret + A[i] * ds( x, c[i], s[i], ni[i] )
    }
    return( ret )
}

# Generates random sample of size n from Gaussian Mixture Model.
# GMM is parametrised using p vector, as described in dgmm.
rgmm <- function (n,p)
{
    m     = length(p)/3
    A     = p[1:m]/sum(p[1:m])
    mu    = p[(m+1):(2*m)]
    sigma = p[(2*m+1):(3*m)]
    x = vector( "numeric", n )
    prev = 0
    for( i in 1:m ) {
        quant = round(n*A[i])
        last  = prev + quant
        if( i == m ) {
            last = n
        }
        x[(prev+1):(last)] = rnorm( length((prev+1):(last)),
                                    mu[i], sigma[i] )
        prev = last
    }
    return( x )
}

# Generates random sample of size n from von Mises Mixture Model.
# vMM is parametrised using p vector.
# Accepts and returns angles in degrees.
# Simulation of random sampling is implemented according to
# Best & Fisher, Efficient Simulation of the von Mises Distribution,
# Journal of the RSS, Series C, 1979, 28, 152-157.
rvmm <- function (n,p)
{
    m  = length(p)/3
    A  = p[1:m]
    mu = deg2rad( p[(m+1):(2*m)] )
    k  = p[(2*m+1):(3*m)]
    x = vector( "numeric", n )
    prev = 0
    for( i in 1:m ) {
        quant = 0
        if( i == m ) {
            quant = n - prev
        } else {
            quant = round(n*A[i])
        }
        last  = prev + quant
        tau = 1 + sqrt( 1 + 4 * k[i]^2 )
        rho = ( tau - sqrt( 2 * tau ) ) / ( 2 * k[i] )
        r   = ( 1 + rho^2 ) / ( 2 * rho )
        c   = vector( "numeric", quant ) * NaN;
        f   = vector( "numeric", quant ) * NaN;
        while( length( c[is.na(c)] ) > 0 ) {
            na_count = length( c[is.na(c)] )
            z  = cos( pi * runif( na_count, 0, 1 ) )
            f[is.na(c)] = ( 1 + r * z ) / ( r + z )
            cn = k[i] * ( r - f[is.na(c)] )
            cn[ log( cn / runif( na_count, 0, 1 ) ) + 1 - cn < 0 ] = NaN
            c[is.na(c)] = cn
        }
        x[(prev+1):(last)] = sign( runif( quant, 0, 1 ) - 0.5 ) *
                             acos( f ) + mu[i]
        prev = last
    }
    return( rad2deg( x ) )
}

rcmm <- function (n,p)
{
    m     = length(p)/3
    A     = p[1:m]/sum(p[1:m])
    mu    = p[(m+1):(2*m)]
    sigma = p[(2*m+1):(3*m)]
    x = numeric( n )
    prev = 0
    for( i in 1:m ) {
        quant = round(n*A[i])
        last  = prev + quant
        if( i == m ) {
            last = n
        }
        x[(prev+1):(last)] = rcauchy( quant, mu[i], sigma[i] )
        prev = last
    }
    return( x )
}

llgmm_conservative <- function (x, p)
{
    n     = length(x)
    m     = length(p)/3
    A     = p[1:m]
    mu    = p[(m+1):(2*m)]
    sigma = p[(2*m+1):(3*m)]
    sum = 0
    for (i in 1:n) {
        sum = sum + log( sum( A / ( sqrt(2*pi) * sigma ) *
                         exp( -(mu-x[i])^2/(2*sigma^2) ) ) )
    }
    return( sum )
}

llsmm <- function( x, p )
{
    n = length(x)
    m = length(p)/4
    A  = p[1:m]/sum(p[1:m])
    c  = p[(m+1):(2*m)]
    s  = p[(2*m+1):(3*m)]
    ni = p[(3*m+1):(4*m)]
    y = numeric( n ) * 0
    for (i in 1:m) {
        y = y + dsmm( x, c( A[i], c[i], s[i], ni[i] ) )
    }
    return( sum( log( y ) ) )
}

llgmm_opposite <- function( x, p )
{
    return( -llgmm( x, p ) )
}

llvmm_opposite <- function( x, p )
{
    return( -llvmm( x, p ) )
}

# Calculate Bayesian Information Criterion (BIC) for any type of mixture
# model. Log-likelihood function has to be provided.
bic <- function( x, p, llf )
{
    return( -2 * llf( x, p ) + (length( p ) - 1) * log( length( x ) ) )
}

# Calculate posterior probability of given number of peaks in
# Gaussian Mixture Model
gmm_size_probability <- function (x, n, method = "SANN")
{
    p = vector( "numeric", n * 3 )
    l = vector( "numeric", n * 3 )
    u = vector( "numeric", n * 3 )
    for (i in 1:n) {
        p[i]     = 1
        p[n+i]   = min(x)+(max(x)-min(x))/n*i-(max(x)-min(x))/n/2
        p[2*n+i] = 1
    }
    f = optim( p,
               llgmm_opposite,
               hessian = TRUE,
               method  = method,
               x = x )
    return( f )
}

# Fit Gaussian Mixture Model to binned data (histogram).
# Lower bounds for mixture proportions and dispersions are fixed in order
# to avoid getting NaNs.
gmm_size_probability_nls <- function (x, n, bins = 100, trace = FALSE)
{
    lower = min( x )
    upper = max( x )
    p = vector( "numeric", n * 3 )
    l = vector( "numeric", n * 3 )
    binsize = (upper-lower)/bins
    for (i in 1:n) {
        p[i]     = 1/n
        p[n+i]   = lower + (upper-lower)/(n+1)*i
        p[2*n+i] = 1
        l[i]     = 0.001
        l[n+i]   = -Inf
        l[2*n+i] = 0.1
    }
    y = vector( "numeric", bins )
    for (i in 1:bins) {
        y[i] = length(x[x >= lower+(i-1)*binsize &
                        x < lower+i*binsize])/length(x)
    }
    leastsq = nls( y ~ dgmm( lower + seq( 0, bins - 1, 1 ) * binsize,
                             theta,
                             normalise_proportions = FALSE),
                   start = list( theta = p ),
                   trace = trace,
                   control = list( warnOnly = TRUE ),
                   algorithm = "port",
                   lower = l )
    par = coef( leastsq )
    prob = factorial( n ) * ( 4*pi )^n * exp( -sum(resid(leastsq)^2)/2 ) /
           (((lower-upper) * 1 * max(par[(2*n+1):(3*n)]))^n *
            sqrt(det(solve(vcov(leastsq)))))
    return( list( p = prob, par = par, residual = sum(resid(leastsq)^2),
                  hessian = solve(vcov(leastsq)), vcov = vcov(leastsq) ) )
}

gmm_fit_kmeans <- function(x, n)
{
    p = vector( "numeric", 3*n )
    km = kmeans( x, n )
    for( i in 1:n ) {
        p[i]     = length(   x[km$cluster==i] )/length(x)
        p[n+i]   = mean(     x[km$cluster==i] )
        p[2*n+i] = sqrt(var( x[km$cluster==1] ))
    }
    return( p )
}

ssd_gradient <- function(x, y, p)
{
    n     = length(x)
    m     = length(p)/3
    A     = p[1:m]
    mu    = p[(m+1):(2*m)]
    sigma = p[(2*m+1):(3*m)]

    grad  = vector( "numeric", length(p) )
    for( i in 1:m ) {
        grad[i]     = 0
        grad[i+m]   = 0
        grad[i+2*m] = 0
        for( k in 1:n ) {
            grad[i]     = grad[i] - 2 * exp( -( x[k] - mu[i] )^2 / (2 * sigma[i]^2) ) / ( (2*pi)^0.5 * sigma[i] ) *
                          ( y[k] - sum( A / ( (2*pi)^0.5 * sigma ) * exp( -( x[k] - mu )^2 / ( 2 * sigma^2 ) ) ) )
            grad[i+m]   = grad[i+m] - 2 * A[i] * ( x[k] - mu[i] ) * exp( -( x[k] - mu[i] )^2 / ( 2 * sigma[i]^2 ) ) *
                          ( y[k] - sum( A / ( (2*pi)^0.5 * sigma ) * exp( -( x[k] - mu )^2 / ( 2 * sigma^2 ) ) ) )
            grad[i+2*m] = grad[i+2*m] + 2 * A[i] * exp( -( x[k] - mu[i] )^2 / ( 2 * sigma[i]^2 ) ) / ( sigma[i]^2 * (2*pi)^0.5 ) *
                          ( 1 - ( x[k] - mu[i] )^2 / sigma[i]^2 ) *
                          ( y[k] - sum( A / ( (2*pi)^0.5 * sigma ) * exp( -( x[k] - mu )^2 / ( 2 * sigma^2 ) ) ) )
        }
    }
    return( grad )
}

pssd_gradient <- function(x, y, p)
{
    grad = ssd_gradient( x, y, p )
    n     = length(x)
    m     = length(p)/3
    A     = p[1:m]
    mu    = p[(m+1):(2*m)]
    sigma = p[(2*m+1):(3*m)]
    for( i in 1:m ) {
        grad[i] = grad[i] + 2 * ( sum( A ) - 1 )
        if( A[i] <= 0 ) {
            grad[i] = grad[i] - exp( -sum( A[A<=0] ) )
        }
        if( mu[i] < min(x) ) {
            grad[i+m] = grad[i+m] - 2 * (min(x) - mu[i])
        }
        if( mu[i] > max(x) ) {
            grad[i+m] = grad[i+m] - 2 * (max(x) - mu[i])
        }
        if( sigma[i] <= 0 ) {
            grad[i+2*m] = grad[i+2*m] - exp( -sum( sigma[sigma<=0] ) )
        }
    }
    return( grad )
}

gradient_descent <- function( gradfn, start, gamma = 0.1, ...,
                              epsilon = 0.01 )
{
    a = start
    while( TRUE ) {
        grad = gradfn( a, ... )
        prev_a = a
        a = a - gamma * grad
        if( sqrt(sum((a-prev_a)^2)) <= epsilon ) {
            break
        }
    }
    return( a )
}

ssd <- function( x, y, p )
{
    return( sum( ( y - dgmm( x, p ) )^2 ) )
}

pssd <- function( x, y, p )
{
    m     = length(p)/3
    A     = p[1:m]
    mu    = p[(m+1):(2*m)]
    sigma = p[(2*m+1):(3*m)]
    sum   = ssd( x, y, c( A[sigma>0], mu[sigma>0], sigma[sigma>0] ) )
    sum   = sum + sum( exp( -sigma[sigma<=0] ) - 1 )
    sum   = sum + sum( exp( -A[A<=0] ) - 1 )
    sum   = sum + ( sum( A ) - 1 )^2
    sum   = sum + sum( ( mu[mu<min(x)] - min(x) )^2 ) +
                  sum( ( mu[mu>max(x)] - max(x) )^2 )
    return( sum )
}

simplex <- function( fn, start, ..., epsilon = 0.000001, alpha = 1,
                     gamma = 2, rho = 0.5, delta = 0.5, trace = F )
{
    A = start
    while( TRUE ) {
        v = vector( "numeric", length( A ) )
        for (i in 1:length( A )) {
            v[i] = fn( A[[i]], ... )
        }
        A  = A[sort( v, index.return = T )$ix]
        v  = sort( v )
        maxdiff = 0
        for (i in 2:length( A )) {
            diff = sqrt( sum( (as.vector(A[[i]])-as.vector(A[[1]]))^2 ) )
            if( diff > maxdiff ) {
                maxdiff = diff;
            }
        }
        if( maxdiff / max( 1, sqrt( sum(as.vector(A[[1]])^2) ) ) <= epsilon ) {
            break;
        }
        x0 = vector( "numeric", length( A[[1]] ) ) * 0  # gravity center
        for (i in 1:(length( A )-1)) {
            x0 = x0 + A[[i]] / (length( A )-1)
        }
        xr = x0 + alpha * ( x0 - A[[length( A )]] )     # reflected point
        fr = fn( xr, ... )

        if( fr < v[1] ) {       # reflected point is the best point so far
            xe = x0 + gamma * ( x0 - A[[length( A )]] ) # expansion
            fe = fn( xe, ... )
            if( fe < fr ) {
                A[[length( A )]] = xe
                if( trace == T ) {
                    cat( "Expanding towards ", xe, " (", fe, ")\n" )
                }
            } else {
                A[[length( A )]] = xr
                if( trace == T ) {
                    cat( "Reflecting towards ", xr, " (", fr, ")\n" )
                }
            }
        } else if( fr >= v[length(v)-1] ) {
            xc = x0 + rho * ( x0 - A[[length( A )]] )   # contraction
            fc = fn( xc, ... )
            if( fc < v[length( A )] ) {
                A[[length( A )]] = xc
                if( trace == T ) {
                    cat( "Contracting towards ", xc, " (", fc, ")\n" )
                }
            } else {
                for (i in 2:length( A )) {  # reduction
                    A[[i]] = A[[1]] + delta * ( A[[i]] - A[[1]] )
                }
                if( trace == T ) {
                    cat( "Reducing\n" )
                }
            }
        } else if( fr < v[length(v)-1] ) {
            A[[length( A )]] = xr       # reflection
            if( trace == T ) {
                cat( "Reflecting towards ", xr, " (", fr, ")\n" )
            }
        } else {
            for (i in 2:length( A )) {  # reduction
                A[[i]] = A[[1]] + delta * ( A[[i]] - A[[1]] )
            }
            if( trace == T ) {
                cat( "Reducing\n" )
            }
        }
        
    }
    v = vector( "numeric", length( A ) )
    for (i in 1:length( A )) {
        v[i] = fn( A[[i]], ... )
    }
    A  = A[sort( v, index.return = T )$ix]
    v  = sort( v )
    return( list( best = A[[1]], score = v[1] ) )
}

rsimplex_start <- function(seed, n, lower, upper)
{
    set.seed( seed )
    l   = list()
    for( i in 1:(length(lower) * n + 1) ) {
        v = vector( "numeric", length(lower) * n ) * 0
        for( j in 1:length(lower) ) {
            v[((j-1)*n+1):(j*n)] = runif( n, lower[j], upper[j] )
        }
        v[1:n] = v[1:n] / sum( v[1:n] )
        l[[i]] = v
    }
    return( l )
}

gmm_fit_hwhm <- function( x, y, n )
{
    a = y
    p = vector( "numeric", 3 * n )
    for( i in 1:n ) {
        maxid = which.max(a)
        mu = x[maxid]
        aleft  = a[1:maxid]
        aright = a[maxid:length(a)]
        xleft  = x[1:maxid]
        xright = x[maxid:length(x)]
        sigma = min( c( x[maxid] - head(xleft[aleft>a[maxid]/2],1),
                        head(xright[aright<a[maxid]/2],1) - x[maxid] ) ) /
                sqrt( 2 * log( 2 ) )
        if( sigma == 0 ) {
            p = c( p[1:(i-1)], p[(n+1):(n+i-1)], p[(2*n+1):(2*n+i-1)] )
            break
        }
        A = a[maxid] / dnorm( mu, mu, sigma )
        a = a - A * dnorm( x, mu, sigma )
        p[i]     = A
        p[i+n]   = mu
        p[i+2*n] = sigma
    }
    return( p )
}

gmm_fit_hwhm_spline_derivatives <- function( x, y )
{
    a = y
    diff = y[2:length(y)]-y[1:(length(y)-1)]
    seq = 2:(length(y)-1)
    peak_positions = seq[diff[1:(length(diff)-1)] > 0 &
                         diff[2:length(diff)] < 0]
    sorted_peak_order = rev( sort( y[peak_positions],
                                   index.return = TRUE )$ix )
    peak_positions = peak_positions[sorted_peak_order]
    n = length( peak_positions )
    p = vector( "numeric", 3 * n )
    for( i in 1:n ) {
        maxid = peak_positions[i]
        mu = x[maxid]
        aleft  = a[1:maxid]
        aright = a[maxid:length(a)]
        xleft  = x[1:maxid]
        xright = x[maxid:length(x)]
        sigma = min( c( x[maxid] - head(xleft[aleft>a[maxid]/2],1),
                        head(xright[aright<a[maxid]/2],1) - x[maxid] ) ) /
                sqrt( 2 * log( 2 ) )
        if( sigma == 0 ) {
            p = c( p[1:(i-1)], p[(n+1):(n+i-1)], p[(2*n+1):(2*n+i-1)] )
            break
        }
        A = a[maxid] / dnorm( mu, mu, sigma )
        a = a - A * dnorm( x, mu, sigma )
        p[i]     = A
        p[i+n]   = mu
        p[i+2*n] = sigma
    }
    return( p )
}

cmm_fit_hwhm_spline_derivatives <- function( x, y )
{
    a = y
    diff = y[2:length(y)]-y[1:(length(y)-1)]
    seq = 2:(length(y)-1)
    peak_positions = seq[diff[1:(length(diff)-1)] > 0 &
                         diff[2:length(diff)] < 0]
    sorted_peak_order = rev( sort( y[peak_positions],
                                   index.return = TRUE )$ix )
    peak_positions = peak_positions[sorted_peak_order]
    n = length( peak_positions )
    p = vector( "numeric", 3 * n )
    for( i in 1:n ) {
        maxid = peak_positions[i]
        mu = x[maxid]
        aleft  = a[1:maxid]
        aright = a[maxid:length(a)]
        xleft  = x[1:maxid]
        xright = x[maxid:length(x)]
        sigma = min( c( x[maxid] - head(xleft[aleft>a[maxid]/2],1),
                        head(xright[aright<a[maxid]/2],1) - x[maxid] ) )
        if( sigma == 0 ) {
            p = c( p[1:(i-1)], p[(n+1):(n+i-1)], p[(2*n+1):(2*n+i-1)] )
            break
        }
        A = a[maxid] / dcauchy( mu, mu, sigma )
        a = a - A * dcauchy( x, mu, sigma )
        p[i]     = A
        p[i+n]   = mu
        p[i+2*n] = sigma
    }
    return( p )
}

smm_fit_em <- function( x, p, ... )
{
    return( smm_fit_em_GNL08( x, p, ... ) )
}

# Fits the distribution of observations with t-distribution (Student's
# distribution) mixture model. Estimation of component parameters is
# implemented according to Fig. 2 in:
# Aeschliman, C.; Park, J. & Kak, A. C. A
# Novel Parameter Estimation Algorithm for the Multivariate t-Distribution
# and Its Application to Computer Vision
# European Conference on Computer Vision 2010, 2010
# https://engineering.purdue.edu/RVL/Publications/Aeschliman2010ANovel.pdf
smm_fit_em_APK10 <- function( x, p, epsilon = c( 1e-6, 1e-6, 1e-6, 1e-6 ),
                              collect.history = FALSE, debug = FALSE )
{
    m  = length(p)/4
    A  = p[1:m]
    c  = p[(m+1):(2*m)]
    s  = p[(2*m+1):(3*m)]
    ni = p[(3*m+1):(4*m)]
    d_A  = c( Inf )
    d_c  = c( Inf )
    d_s  = c( Inf )
    d_ni = c( Inf )
    steps = 0
    history = list()
    if( collect.history == TRUE ) {
        history[[1]] = p
    }
    while( length( d_A[ d_A  > epsilon[1]] ) > 0 ||
           length( d_c[ d_c  > epsilon[2]] ) > 0 ||
           length( d_s[ d_s  > epsilon[3]] ) > 0 ||
           length( d_ni[d_ni > epsilon[4]] ) > 0 ) {
        prev_A  = A
        prev_c  = c
        prev_s  = s
        prev_ni = ni
        q = dsmm( x, c( A, c, s, ni ) )
        for( j in 1:m ) {
            h = A[j] * ds( x, c[j], s[j], ni[j] )

            w = h / q
            A[j] = sum( w ) / length( x )

            ord = order( x )
            xo = x[ord]
            ho = w[ord] / sum( w )

            c[j] = wmedian( xo, ho ) + 1e-6

            z = log( ( x - c[j] )^2 )

            zbar = sum( z * w ) / sum( w )

            b = sum( ( z - zbar )^2 * w ) / sum( w ) - trigamma( 0.5 )

            ni[j] = ( 1 + sqrt( 1 + 4 * b ) ) / b
            s[j] = exp( zbar - log( ni[j] ) + digamma( ni[j] / 2 ) -
                        digamma( 0.5 ) )
            if( debug == TRUE ) {
                cat( A[j], " ", c[j], " ", s[j], " ", ni[j], " " )
            }
        }
        if( debug == TRUE ) {
            cat( "\n" )
        }
        d_A  = abs( A  - prev_A  )
        d_c  = abs( c  - prev_c  )
        d_s  = abs( s  - prev_s  )
        d_ni = abs( ni - prev_ni )
        steps = steps + 1
        if( collect.history == TRUE ) {
            history[[steps+1]] = c( A, c, s, ni )
        }
        if( length( d_A[ is.na(d_A) ] ) +
            length( d_c[ is.na(d_c) ] ) +
            length( d_s[ is.na(d_s) ] ) +
            length( d_ni[is.na(d_ni)] ) > 0 ) {
            break
        }
    }
    l = list( p = c( A, c, s, ni ), steps = steps )
    if( collect.history == TRUE ) {
        l$history = history
    }
    return( l )
}

# Coeficients for digamma function approximation, that contains first
# eight non-zero members of asymptotic expression for digamma(x).
# Taken from Wikipedia (see "Computation and approximation"):
# https://en.wikipedia.org/w/index.php?title=Digamma_function&oldid=708779689
digamma_approx_coefs = c( 1/2, 1/12, 0, -1/120, 0, 1/252, 0,
                          -1/240, 0, 1/660, 0, -691/32760, 0, 1/12 )

smm_fit_em_GNL08 <- function( x, p, epsilon = c( 1e-6, 1e-6, 1e-6, 1e-6 ),
                              collect.history = FALSE, debug = FALSE,
                              min.sigma = 1e-256, min.ni = 1e-256,
                              max.df = 1000, max.steps = Inf,
                              polyroot.solution = 'jenkins_taub',
                              convergence = abs_convergence,
                              unif.component = FALSE )
{
    m  = length(p)/4
    A  = p[1:m]
    c  = p[(m+1):(2*m)]
    s  = p[(2*m+1):(3*m)]
    ni = p[(3*m+1):(4*m)]
    prev_A  = rep( Inf, m )
    prev_c  = rep( Inf, m )
    prev_s  = rep( Inf, m )
    prev_ni = rep( Inf, m )
    steps = 0
    history = list()
    if( collect.history ) {
        history[[1]] = p
    }
    while( steps < max.steps &&
           steps == 0 ||
           !convergence( c( A, c, s, ni ),
                         c( prev_A, prev_c, prev_s, prev_ni ), epsilon ) ) {
        prev_A  = A
        prev_c  = c
        prev_s  = s
        prev_ni = ni
        q = dsmm( x, c( A, c, s, ni ) )
        if( unif.component ) {
            # Allows additional component with uniform distribution for
            # the modelling of outliers as suggested in:
            # Cousineau, D. & Chartier, S.
            # Outliers detection and treatment: a review
            # International Journal of Psychological Research,
            # 2010, 3, 58-67
            # http://revistas.usb.edu.co/index.php/IJPR/article/view/844

            q = q + ( 1 - sum( A ) ) * dunif( x, min(x), max(x) )
        }
        for( j in 1:m ) {
            h = A[j] * ds( x, c[j], s[j], ni[j] )
            z = h / q

            u = ( ni[j] + 1 ) / ( ni[j] + ( ( x - c[j] ) / s[j] )^2 )

            A[j] = sum( z ) / length( x )
            c[j] = sum( z * u * x ) / sum( z * u )
            s[j] = sqrt( sum( z * u * ( x - c[j] )^2 ) / sum( z ) )

            # Solution of Eqn. 17 is implemented via digamma function
            # approximation using asymptotic expression of digamma(x).
            # Jenkins-Taub (implemented in R's polyroot() function) or
            # Newton-Raphson (implemented here) algorithm is used to find
            # the roots of the polynomial. For Jenkins-Taub, a positive
            # real root (should be single) is chosen as a solution.

            cl = length( digamma_approx_coefs )
            p = sum( rep( 2 / ( ni[j] + 1 ), cl )^(1:cl) *
                     digamma_approx_coefs )
            polynome = c( sum( z * ( log( u ) - u ) ) / sum( z ) - p + 1,
                          digamma_approx_coefs )

            roots = switch(
                polyroot.solution,
                jenkins_taub   = polyroot( polynome ),
                newton_raphson = polyroot_NR( polynome, init = 2/ni[j] ),
                NaN )

            ni[j] = 2 / switch(
                polyroot.solution,
                jenkins_taub   = Re(roots[abs(Im(roots)) < 1e-10 &
                                          Re(roots) > 1e-10]),
                newton_raphson = roots,
                NaN )

            if( ni[j] > max.df ) {
                ni[j] = max.df
            }

            if( debug ) {
                cat( A[j], " ", c[j], " ", s[j], " ", ni[j], " " )
            }
        }
        if( debug ) {
            cat( "\n" )
        }
        steps = steps + 1
        if( collect.history ) {
            history[[steps+1]] = c( A, c, s, ni )
        }
        if( length( A[ is.na(A) ] ) +
            length( c[ is.na(c) ] ) +
            length( s[ is.na(s) ] ) +
            length( ni[is.na(ni)] ) +
            length(  s[s  <= min.sigma] ) +
            length( ni[ni <= min.ni] ) > 0 ) {
            A  = A  * NaN
            c  = c  * NaN
            s  = s  * NaN
            ni = ni * NaN
            break
        }
    }

    l = list( p = c( A, c, s, ni ), steps = steps )
    if( collect.history ) {
        l$history = history
    }
    return( l )
}

# Greedy EM algorithm for Student's t-distribution mixture fitting.
# Implemented according to:
# Chen, S.; Wang, H. & Luo, B.
# Greedy EM Algorithm for Robust T-Mixture Modeling
# Third International Conference on Image and Graphics (ICIG'04),
# Institute of Electrical & Electronics Engineers (IEEE), 2004, 548--551
smm_fit_em_CWL04 <- function( x, p, collect.history = FALSE,
                              debug = FALSE, ... )
{
    bic_prev = Inf
    prev_p = p
    m = length(p) / 4
    run = TRUE
    history = list()

    while( run ) {
        if( debug ) {
            cat( "Starting EM with", m, "components\n" )
        }
        fit = smm_fit_em( x, p, ... )
        p = fit$p
        bic_now = bic( x, p, llsmm )
        if( debug ) {
            cat( "Achieving fit with BIC =", bic_now, "\n" )
        }
        if( bic_prev > bic_now ) {
            bic_prev = bic_now
            kldivs = numeric( m )
            for( i in 1:m ) {
                kldivs[i] = kldiv( x, p, i )
            }
            split = which.max( kldivs )
            if( debug ) {
                cat( "Splitting component", split, "\n" )
            }
            if( collect.history ) {
                history[[m]] = p
            }
            s = smm_split_component( p[0:3*m+split] )
            prev_p = p
            p = c( p[0*m+sort( c( 1:m, split ) )],
                   p[1*m+sort( c( 1:m, split ) )],
                   p[2*m+sort( c( 1:m, split ) )],
                   p[3*m+sort( c( 1:m, split ) )] )
            m = m + 1
            p[0*m+split+0:1] = s[1:2]
            p[1*m+split+0:1] = s[3:4]
            p[2*m+split+0:1] = s[5:6]
            p[3*m+split+0:1] = s[7:8]
        } else {
            run = FALSE
            p = prev_p
            if( debug ) {
                cat( "Stopping on convergence criterion\n" )
            }
        }
    }
    l = list( p = p )
    if( collect.history ) {
        l$history = history
    }
    return( l )
}

# Fits the distribution of observations with t-distribution (Student's
# distribution) mixture model. Implemented according to the Batch
# Approximation Algorithm, as given in Fig. 2 in:
# Aeschliman, C.; Park, J. & Kak, A. C. A
# Novel Parameter Estimation Algorithm for the Multivariate t-Distribution
# and Its Application to Computer Vision
# European Conference on Computer Vision 2010, 2010
# https://engineering.purdue.edu/RVL/Publications/Aeschliman2010ANovel.pdf
s_fit_primitive <- function( x )
{
    xbar = median( x )
    z = log( ( x - xbar )^2 )
    zbar = sum( z ) / length( x )
    b = sum( ( z - zbar )^2 ) / length( x ) - trigamma( 0.5 )
    ni = ( 1 + sqrt( 1 + 4 * b ) ) / b
    alpha = exp( zbar - log( ni ) + digamma( ni / 2 ) - digamma( 0.5 ) )
    return( c( xbar, alpha, ni ) )
}

mk_fit_images <- function( h, l, prefix = "img_" ) {
    maxstrlen = ceiling( log( length( l ) ) / log( 10 ) )
    for( i in 1:length( l ) ) {
        fname = paste( prefix,
                       sprintf( paste( "%0", maxstrlen, "d", sep = "" ), i ),
                       ".png", sep = "" )
        png( filename = fname )
        plot( h )
        lines( h$mids,
               sum( h$counts ) * dgmm( h$mids, l[[i]] )/sum( h$density ),
               col = "red", lwd = 2 )
        dev.off()
    }
}

gmm_init_vector_kmeans <- function( x, m ) {
    start = numeric( 3 * m )
    if( min(x) == max(x) ) {
        start[1:m]           = 1/m
        start[(m+1):(2*m)]   = min(x) + (1:m)*(max(x)-min(x))/(m+1)
        start[(2*m+1):(3*m)] = (max(x)-min(x))/(m+1)/6
    } else {
        k = kmeans( x, m )
        start[1:m]           = k$size / length( x )
        start[(m+1):(2*m)]   = k$centers
        start[(2*m+1):(3*m)] = sqrt( k$withinss / k$size )
    }
    return( start )
}

gmm_init_vector_quantile <- function( x, m ) {
    sorted = sort( x )
    start = numeric( 3 * m )
    start[1:m]           = 1/m
    start[(m+1):(2*m)]   = sorted[floor(length(x) / (m+1) * 1:m)]
    start[(2*m+1):(3*m)] = sqrt( sum( (x-mean(x))^2 ) / length(x) )
    return( start )
}

cmm_init_vector_kmeans <- function( x, m, iter.cauchy = 20 ) {
    start = numeric( 3 * m )
    if( min(x) == max(x) ) {
        start[1:m]           = 1/m
        start[(m+1):(2*m)]   = min(x) + (1:m)*(max(x)-min(x))/(m+1)
        start[(2*m+1):(3*m)] = 1
    } else {
        k = kmeans( x, m )
        start[1:m] = k$size / length( x )
        start[(m+1):(2*m)]   = 0
        start[(2*m+1):(3*m)] = 1
        for( n in 1:m ) {
            for( i in 1:iter.cauchy ) {
                u    = (x[k$cluster == n] - start[m+n]) / start[2*m+n]
                e0k  = sum( 1 / (1 + u^2) ) / k$size[n]
                e1k  = sum( u / (1 + u^2 ) ) / k$size[n]
                start[m+n] = start[m+n] + start[2*m+n] * e1k / e0k
                start[2*m+n] = start[2*m+n] * sqrt( 1/e0k - 1 )
            }
        }
    }
    return( start )
}

smm_init_vector <- function( x, m ) {
    start = numeric( 4 * m )
    start[1:m]           = 1/m
    start[(m+1):(2*m)]   = min(x) + (1:m)*(max(x)-min(x))/(m+1)
    start[(2*m+1):(3*m)] = 1
    start[(3*m+1):(4*m)] = 1
    return( start )
}

smm_init_vector_kmeans <- function( x, m ) {
    start = numeric( 4 * m )
    if( min(x) == max(x) ) {
        start[1:m]           = 1/m
        start[(m+1):(2*m)]   = min(x) + (1:m)*(max(x)-min(x))/(m+1)
        start[(2*m+1):(3*m)] = 1
        start[(3*m+1):(4*m)] = 1
    } else {
        k = kmeans( x, m )
        start[1:m] = k$size / length( x )
        for( i in 1:m ) {
            p = s_fit_primitive( x[k$cluster==i] )
            start[1*m+i] = p[1]
            start[2*m+i] = p[2]
            start[3*m+i] = p[3]
        }
    }
    return( start )
}

gmm_merge_components <- function( x, p, i, j ) {
    P = matrix( p, ncol = 3 )
    A = P[i,1] + P[j,1]

    # Performing an iteration of EM to find new mean and sd
    q = vector( "numeric", length( x ) ) * 0
    for( k in 1:nrow( P ) ) {
        q = q + P[k,1] * dnorm( x, P[k,2], P[k,3] )
    }
    h = ( P[i,1] * dnorm( x, P[i,2], P[i,3] ) +
          P[j,1] * dnorm( x, P[j,2], P[j,3] ) ) / q
    mu = sum( x * h ) / sum( h )
    sigma = sqrt( sum( h * ( x - mu ) ^ 2 ) / sum( h ) )

    P[i,] = c( A, mu, sigma )
    return( as.vector( P[setdiff( 1:nrow( P ), j ),] ) )
}

# Splits a component of Student's t-distribution mixture. Implemented
# according to Eqns. 30--36 of:
# Chen, S.-B. & Luo, B.
# Robust t-mixture modelling with SMEM algorithm
# Proceedings of 2004 International Conference on Machine Learning and
# Cybernetics (IEEE Cat. No.04EX826),
# Institute of Electrical & Electronics Engineers (IEEE), 2004, 6, 3689--3694
smm_split_component <- function( p, alpha = 0.5, beta = 0.5, u = 0.5 ) {
    A  = c( alpha, 1 - alpha ) * p[1]
    c = p[2] + c( -sqrt( A[2]/A[1] ), sqrt( A[1]/A[2] ) ) * u *
               sqrt( p[3] )
    s = p[1] * p[3] * (1 - (p[4] - 2) * u^2 / p[4]) *
        c( beta, 1 - beta ) / A
    return( c( A, c, s, p[4], p[4] ) )
}

plot_circular_hist <- function( x, breaks = 72, ball = 0.5, ... ) {
    xx = numeric( breaks * 3 + 1 )
    yy = numeric( breaks * 3 + 1 ) * 0
    xx[(1:breaks)*3] = 1:breaks * 2*pi / breaks
    for( i in 1:breaks ) {
        xx[((i-1)*3+1):((i-1)*3+2)] = (i-1) * 2*pi / breaks
        yy[((i-1)*3+2):(i*3)] =
            length( x[x >= 360 / breaks * (i-1) & x < 360 / breaks * i] )
    }
    yy = (yy / max(yy)) * (1-ball) + ball
    plot( yy * cos( xx ), yy * sin( xx ), type = "l", asp = 1,
          ann = FALSE, axes = FALSE, ... )
    lines( ball * cos( seq( 0, 2*pi, pi/180 ) ),
           ball * sin( seq( 0, 2*pi, pi/180 ) ) )
}

deg2rad <- function( x ) {
    return( x * pi / 180 )
}

rad2deg <- function( x ) {
    return( x * 180 / pi )
}

kmeans_circular <- function( x, centers, iter.max = 10 ) {
    centers = sort( deg2rad( centers ) )
    n = length( centers )
    x = deg2rad( x )
    for( i in 1:iter.max ) {
        cluster = numeric(n) * 0
        for( j in 2:(n-1) ) {
            cluster[x >= (centers[j] + centers[j-1])/2 &
                    x <  (centers[j+1] + centers[j])/2] = j
        }
        midpoint = (centers[n] - 2*pi + centers[1])/2
        if( midpoint < 0 ) {
            cluster[x < (centers[1] + centers[2])/2] = 1
            cluster[x >= midpoint + 2*pi] = 1
            cluster[x <  midpoint + 2*pi &
                    x >= (centers[n-1] + centers[n])/2] = n
            x[x >= midpoint + 2*pi] = x[x >= midpoint + 2*pi] - 2*pi
        } else {
            cluster[x >= (centers[n-1] + centers[n])/2] = n
            cluster[x <  midpoint] = n
            cluster[x < (centers[1] + centers[2])/2 & x >= midpoint] = 1
            x[x < midpoint] = x[x < midpoint] + 2*pi
        }
        for( j in 1:n ) {
            centers[j] = sum( x[cluster == j] ) /
                         length( cluster[cluster == j] )
        }
        centers[centers<0]    = centers[centers<0]    + 2*pi
        centers[centers>2*pi] = centers[centers>2*pi] - 2*pi
        x[x < 0]    = x[x < 0]    + 2*pi
        x[x > 2*pi] = x[x > 2*pi] - 2*pi
        centers = sort( centers )
    }
    return( rad2deg( centers ) )
}

# Weighted median function, implemented according to Wikipedia:
# https://en.wikipedia.org/w/index.php?title=Weighted_median&oldid=690896947
wmedian <- function( x, w, start = 1, end = length( x ) )
{
    # base case for single element
    if( start == end ) {
        return( x[start] )
    }

    # base case for two elements
    # make sure we return lower median
    if( end - start == 1 ) {
        if( w[start] >= w[end] ) {
            return( x[start] )
        } else {
            return( x[end] )
        }
    }

    # partition around center pivot
    q = round( ( start + end ) / 2 )

    w_left  = sum( w[start:(q-1)] )
    w_right = sum( w[(q+1):end] )

    if( w_left < 0.5 & w_right < 0.5 ) {
        return( x[q] )
    }

    if( w_left > w_right ) {
        w[q] = w[q] + w_right
        return( wmedian( x, w, start, q ) )
    } else {
        w[q] = w[q] + w_left
        return( wmedian( x, w, q, end ) )
    }
}

digamma_approx <- function( x )
{
    cl = length( digamma_approx_coefs )
    ret = numeric( length( x ) )

    for( i in 1:length(x) ) {
        ret[i] = log( x[i] ) -
                 sum( rep( 1/x[i], cl )^(1:cl) * digamma_approx_coefs )
    }

    return( ret )
}

# Kullback--Leibler divergence, using Dirac's delta function, implemented
# according to:
# Chen, S.; Wang, H. & Luo, B.
# Greedy EM Algorithm for Robust T-Mixture Modeling
# Third International Conference on Image and Graphics (ICIG'04),
# Institute of Electrical & Electronics Engineers (IEEE), 2004, 548-551
kldiv <- function( x, p, k )
{
    m  = length( p ) / 4
    A  = p[k]
    c  = p[m+k]
    s  = p[2*m+k]
    ni = p[3*m+k]
    z = A * ds( x, c, s, ni ) / dsmm( x, p )
    kld = 0
    for( i in unique(x) ) {
        pk = z[x==i] / sum( z )
        fk = ds( i, c, s, ni )
        kld = kld + pk * log( pk / fk )
    }
    return( kld )
}

bhattacharyya_dist <- function( mu1, mu2, sigma1, sigma2 )
{
    return( log( sum( c( sigma1, sigma2 )^2 /
                      c( sigma2, sigma1 )^2, 2 ) / 4 ) / 4 +
            ( mu1 - mu2 )^2 / ( 4 * ( sigma1^2 + sigma2^2 ) ) )
}

abs_convergence <- function( p_now, p_prev, epsilon = 1e-6 )
{
    if( length( epsilon ) > 1 && length( epsilon ) < length( p_now ) ) {
        n = length( p_now ) / length( epsilon )
        epsilon_now = numeric( 0 )
        for( i in length( epsilon ) ) {
            epsilon_now = c( epsilon_now, rep( epsilon[i], n ) )
        }
    }
    has_converged = all( abs( p_now - p_prev ) <= epsilon )
    if( is.na( has_converged ) ) {
        has_converged = TRUE
    }
    return( has_converged )
}

ratio_convergence <- function( p_now, p_prev, epsilon = 1e-6 )
{
    if( length( epsilon ) > 1 && length( epsilon ) < length( p_now ) ) {
        n = length( p_now ) / length( epsilon )
        epsilon_now = numeric( 0 )
        for( i in length( epsilon ) ) {
            epsilon_now = c( epsilon_now, rep( epsilon[i], n ) )
        }
    }
    has_converged = all( abs( p_now - p_prev ) / p_prev <= epsilon );
    if( is.na( has_converged ) ) {
        has_converged = TRUE
    }
    return( has_converged )
}

plot_density <- function( x, cuts = 400, main = '', model, density_f,
                          filename = '/dev/stdout',
                          width, height, obs_good = c(), obs_bad = c(),
                          scale_density = FALSE )
{
    png( filename, width = width, height = height )
    h = hist( x, cuts, main = main,
              xlim = c( min( c( x, obs_bad ) ), max( c( x, obs_bad ) ) ) )
    xmids = seq( min( c( x, obs_bad )  ),
                 max( c( x, obs_bad )  ),
                 h$mids[2] - h$mids[1] )
    density = do.call( density_f, list( xmids, model ) )
    if( scale_density == TRUE ) {
        density = density / max( density ) * max( h$counts )
    } else {
        density = length(x) / sum( h$density ) * density
    }
    lines( xmids, density, lwd = 2, col = 'green' )
    if( length( obs_good ) > 0 ) {
        rug( obs_good, lwd = 2, col = 'green' )
    }
    if( length( obs_bad ) > 0 ) {
        rug( obs_bad,  lwd = 2, col = 'red' )
    }
    dev.off()
}
