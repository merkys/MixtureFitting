#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>

double bessi0( double x ) {
    return bessel_i( x, 0, 1 );
}

double deg2rad( double x ) {
    return x * M_PI / 180;
}

double rad2deg( double x ) {
    return x * 180 / M_PI;
}

void dgmm( double *x, int *xlength,
           double *p, int *plength,
           double *ret )
{
    int m = *plength / 3;
    double mu[m];
    double factor[m];
    double divisor[m];
    int i;
    int j;
    double sqrtdblpi = sqrt( 2 * M_PI );
    for( i = 0; i < m; i++ ) {
        double A = p[i];
        mu[i]    = p[m+i];
        double sigma = p[2*m+i];
        factor[i]  = A / (sigma * sqrtdblpi);
        divisor[i] = 2 * sigma * sigma;
    }
    for( i = 0; i < *xlength; i++ ) {
        ret[i] = 0.0;
        for( j = 0; j < m; j++ ) {
            double diff = x[i] - mu[j];
            ret[i] = ret[i] + factor[j] * exp( -(diff * diff) / divisor[j] );
        }
    }
}

void llgmm( double *x, int *xlength,
            double *p, int *plength,
            double *ret )
{
    int m = *plength / 3;
    double * dgmms = calloc( *xlength, sizeof( double ) );
    if( !dgmms ) {
        error( "can not allocate memory" );
    }
    dgmm( x, xlength, p, plength, dgmms );
    int i;
    *ret = 0.0;
    for( i = 0; i < *xlength; i++ ) {
        *ret = *ret + log( dgmms[i] );
    }
    free( dgmms );
}

void gmm_fit_em( double *x, int *xlength,
                 double *p, int *plength,
                 double *epsilon,
                 int *debug,
                 double *ret,
                 int *steps )
{
    int m = *plength / 3;
    double A[m];
    double mu[m];
    double sigma[m];
    int i;
    int j;
    for( i = 0; i < m; i++ ) {
        A[i]  = p[i];
        mu[i] = p[m+i];
        sigma[i] = p[2*m+i];
    }
    double sqrtdblpi = sqrt( 2 * M_PI );
    double * q = calloc( *xlength, sizeof( double ) );
    if( !q ) {
        error( "can not allocate memory" );
    }
    int run = 1;
    *steps = 0;
    while( run == 1 ) {
        int is_nan = 0;
        run = 0;
        for( i = 0; i < *xlength; i++ ) {
            q[i] = 0.0;
            for( j = 0; j < m; j++ ) {
                double diff = x[i] - mu[j];
                q[i] = q[i] + A[j] / (sigma[j] * sqrtdblpi) *
                       exp( -(diff * diff) / (2 * sigma[j] * sigma[j]) );
            }
        }
        for( j = 0; j < m; j++ ) {
            double sumh = 0.0;
            double sumhdiff = 0.0;
            double sumhprod = 0.0;

            double factor = A[j] / (sigma[j] * sqrtdblpi);
            double twosqsigma = 2 * sigma[j] * sigma[j];

            for( i = 0; i < *xlength; i++ ) {
                double diff = x[i] - mu[j];
                double sqdiff = diff * diff;
                double weight = factor * exp( -sqdiff / twosqsigma ) / q[i];
                sumh = sumh + weight;
                sumhdiff = sumhdiff + weight * sqdiff;
                sumhprod = sumhprod + weight * x[i];
            }
            double prev_A = A[j];
            A[j] = sumh / *xlength;
            double prev_sigma = sigma[j];
            sigma[j] = sqrt( sumhdiff / sumh );
            double prev_mu = mu[j];
            mu[j] = sumhprod / sumh;
            if( *debug > 0 ) {
                printf( "%f %f %f ", A[j], mu[j], sigma[j] );
            }
            if( fabs( A[j]     - prev_A )     > epsilon[0] ||
                fabs( mu[j]    - prev_mu )    > epsilon[1] ||
                fabs( sigma[j] - prev_sigma ) > epsilon[2] ) {
                run = 1;
            }
            if( isnan( sigma[j] ) ) {
                is_nan = 1;
            }
        }
        if( is_nan == 1 ) {
            run = 0;
        }
        if( *debug > 0 ) {
            printf( "\n" );
        }
        *steps = *steps + 1;
    }
    for( i = 0; i < m; i++ ) {
        ret[i] = A[i];
        ret[m+i] = mu[i];
        ret[2*m+i] = sigma[i];
    }
    free( q );
}

void gmm_init_vector( double *x, int *xlength,
                      int *m,
                      double *ret )
{
    int i;
    double min = x[0];
    double max = x[0];
    for( i = 1; i < *xlength; i++ ) {
        if( x[i] < min ) { min = x[i]; }
        if( x[i] > max ) { max = x[i]; }
    }
    for( i = 0; i < *m; i++ ) {
        ret[i] = 1.0/(*m);
        ret[*m+i] = min + (i+1)*(max-min)/(*m+1);
        ret[2*(*m)+i] = (max-min)/(*m+1)/6;
    }
}

int main( int argc, char *argv[], char *env[] )
{
    return( 0 );
}

void dcmm( double *x, int *xlength,
           double *p, int *plength,
           double *ret )
{
    int m = *plength / 3;
    double mu[m];
    int i;
    int j;
    double sqrtdblpi = sqrt( 2 * M_PI );
    for( i = 0; i < *xlength; i++ ) {
        ret[i] = 0.0;
        for( j = 0; j < m; j++ ) {
            double normdiff = ( x[i] - p[j+m] ) / p[j+2*m];
            ret[i] = ret[i] + p[j] / ( M_PI * p[j+2*m] * ( 1 + normdiff*normdiff ) );
        }
    }
}

void llcmm( double *x, int *xlength,
            double *p, int *plength,
            double *ret )
{
    int m = *plength / 3;
    double * dcmms = calloc( *xlength, sizeof( double ) );
    if( !dcmms ) {
        error( "can not allocate memory" );
    }
    dcmm( x, xlength, p, plength, dcmms );
    int i;
    *ret = 0.0;
    for( i = 0; i < *xlength; i++ ) {
        *ret = *ret + log( dcmms[i] );
    }
    free( dcmms );
}

void cmm_fit_em( double *x, int *xlength,
                 double *p, int *plength,
                 double *epsilon,
                 int *itercauchy,
                 int *debug,
                 double *ret,
                 int *steps )
{
    int m = *plength / 3;
    double A[m];
    double c[m];
    double s[m];
    int i;
    int j;
    int k;
    for( i = 0; i < m; i++ ) {
        A[i] = p[i];
        c[i] = p[m+i];
        s[i] = p[2*m+i];
    }
    double * q = calloc( *xlength, sizeof( double ) );
    if( !q ) {
        error( "can not allocate memory" );
    }
    double * h = calloc( *xlength, sizeof( double ) );
    if( !h ) {
        error( "can not allocate memory" );
    }
    int run = 1;
    *steps = 0;
    while( run == 1 ) {
        int is_nan = 0;
        run = 0;
        for( i = 0; i < *xlength; i++ ) {
            q[i] = 0.0;
            for( j = 0; j < m; j++ ) {
                double diff = x[i] - c[j];
                q[i] = q[i] + A[j] * s[j] / ( M_PI * ( diff*diff + s[j]*s[j] ) );
            }
        }
        for( j = 0; j < m; j++ ) {
            double sumh = 0.0;

            double factor = A[j] * s[j] / M_PI;
            double sqs = s[j] * s[j];

            for( i = 0; i < *xlength; i++ ) {
                double diff = x[i] - c[j];
                h[i] = factor / ( diff*diff + sqs ) / q[i];
                sumh = sumh + h[i];
            }
            double prev_A = A[j];
            A[j] = sumh / *xlength;
            double prev_c = c[j];
            double prev_s = s[j];
            int is_converged = 0;
            k = 0;
            while( is_converged == 0 && k < *itercauchy ) {
                double hdiv  = 0.0;
                double hprod = 0.0;
                for( i = 0; i < *xlength; i++ ) {
                    double chi = (x[i] - c[j]) / s[j];
                    double sqchi = chi * chi;
                    hdiv  = hdiv  + h[i] / (1+sqchi);
                    hprod = hprod + h[i] * chi/(1+sqchi);
                }
                double e0k = hdiv  / sumh;
                double e1k = hprod / sumh;
                c[j] = c[j] + s[j] * e1k / e0k;
                s[j] = s[j] * sqrt( 1 / e0k - 1 );
                if( fabs( c[j] - prev_c ) < 1e-6 &&
                    fabs( s[j] - prev_s ) < 1e-6 ) {
                    is_converged = 1;
                }
                k++;
            }
            if( *debug > 0 ) {
                printf( "%f %f %f ", A[j], c[j], s[j] );
            }
            if( fabs( A[j] - prev_A ) > epsilon[0] ||
                fabs( c[j] - prev_c ) > epsilon[1] ||
                fabs( s[j] - prev_s ) > epsilon[2] ) {
                run = 1;
            }
            if( isnan( s[j] ) ) {
                is_nan = 1;
            }
        }
        if( is_nan == 1 ) {
            run = 0;
        }
        if( *debug > 0 ) {
            printf( "\n" );
        }
        *steps = *steps + 1;
    }
    for( i = 0; i < m; i++ ) {
        ret[i] = A[i];
        ret[m+i] = c[i];
        ret[2*m+i] = s[i];
    }
    free( q );
    free( h );
}

void cmm_init_vector( double *x, int *xlength,
                    int *m,
                    double *ret )
{
    int i;
    double min = x[0];
    double max = x[0];
    for( i = 1; i < *xlength; i++ ) {
        if( x[i] < min ) { min = x[i]; }
        if( x[i] > max ) { max = x[i]; }
    }
    for( i = 0; i < *m; i++ ) {
        ret[i] = 1.0/(*m);
        ret[*m+i] = min + (i+1)*(max-min)/(*m+1);
        ret[2*(*m)+i] = 1;
    }
}

void dvmm( double *x, int *xlength,
           double *p, int *plength,
           double *ret )
{
    int m = *plength / 3;
    double denom[m];
    int i;
    int j;
    for( j = 0; j < m; j++ ) {
        denom[j] = 2 * M_PI * bessi0( p[2*m+j] );
    }
    for( i = 0; i < *xlength; i++ ) {
        ret[i] = 0.0;
        for( j = 0; j < m; j++ ) {
            double diffcos = cos( deg2rad( x[i] - p[m+j] ) );
            ret[i] = ret[i] + p[j] * exp( p[2*m+j] * diffcos ) / denom[j];
        }
    }
}

void llvmm( double *x, int *xlength,
            double *p, int *plength,
            double *ret )
{
    double * dvmms = calloc( *xlength, sizeof( double ) );
    if( !dvmms ) {
        error( "can not allocate memory" );
    }
    dvmm( x, xlength, p, plength, dvmms );
    int i;
    *ret = 0.0;
    for( i = 0; i < *xlength; i++ ) {
        *ret = *ret + log( dvmms[i] );
    }
    free( dvmms );
}

void vmm_fit_em_by_diff( double *x, int *xlength,
                         double *p, int *plength,
                         double *epsilon,
                         int *debug,
                         double *ret,
                         int *steps )
{
    int m = *plength / 3;
    double A[m];
    double mu[m];
    double k[m];
    double denom[m];
    int i;
    int j;
    for( i = 0; i < m; i++ ) {
        A[i] = p[i];
        mu[i] = p[m+i];
        k[i] = p[2*m+i];
        denom[i] = 2.0 * M_PI * bessi0( k[i] );
    }
    double * q = calloc( *xlength, sizeof( double ) );
    if( !q ) {
        error( "can not allocate memory" );
    }
    double * h = calloc( *xlength, sizeof( double ) );
    if( !h ) {
        error( "can not allocate memory" );
    }
    int run = 1;
    *steps = 0;
    while( run == 1 ) {
        int is_nan = 0;
        run = 0;
        for( i = 0; i < *xlength; i++ ) {
            q[i] = 0.0;
            for( j = 0; j < m; j++ ) {
                double diffcos = cos( deg2rad( x[i] - mu[j] ) );
                q[i] = q[i] + A[j] * exp( k[j] * diffcos ) / denom[j];
            }
        }
        for( j = 0; j < m; j++ ) {
            double sumh = 0.0;
            double sumsin = 0.0;
            double sumcos = 0.0;
            for( i = 0; i < *xlength; i++ ) {
                double diffcos = cos( deg2rad( x[i] - mu[j] ) );
                h[i] = A[j] * exp( k[j] * diffcos ) / ( denom[j] * q[i] );
                sumh = sumh + h[i];
                sumsin = sumsin + sin( deg2rad(x[i]) ) * h[i];
                sumcos = sumcos + cos( deg2rad(x[i]) ) * h[i];
            }
            double prev_A = A[j];
            A[j] = sumh / *xlength;
            double prev_mu = mu[j];
            mu[j] = rad2deg( atan2( sumsin, sumcos ) );
            double Rbar = sqrt( sumsin*sumsin + sumcos*sumcos ) / sumh;
            double prev_k = k[j];
            k[j] = ( 2.0 * Rbar - Rbar*Rbar*Rbar ) / ( 1.0 - Rbar*Rbar );
            denom[j] = 2.0 * M_PI * bessi0( k[j] );
            if( *debug > 0 ) {
                printf( "%f %f %f ", A[j], mu[j], k[j] );
            }
            if( fabs( A[j]  - prev_A )  > epsilon[0] ||
                fabs( mu[j] - prev_mu ) > epsilon[1] ||
                fabs( k[j]  - prev_k )  > epsilon[2] ) {
                run = 1;
            }
            if( isnan( k[j] ) ) {
                is_nan = 1;
            }
        }
        if( is_nan == 1 ) {
            run = 0;
        }
        if( *debug > 0 ) {
            printf( "\n" );
        }
        *steps = *steps + 1;
    }
    for( i = 0; i < m; i++ ) {
        ret[i] = A[i];
        ret[m+i] = mu[i];
        ret[2*m+i] = k[i];
    }
    free( q );
    free( h );
}

void vmm_fit_em_by_ll( double *x, int *xlength,
                       double *p, int *plength,
                       double *epsilon,
                       int *debug,
                       double *ret,
                       int *steps )
{
    int m = *plength / 3;
    double A[m];
    double mu[m];
    double k[m];
    double denom[m];
    double prev_llog;
    llvmm( x, xlength, p, plength, &prev_llog );
    int i;
    int j;
    for( i = 0; i < m; i++ ) {
        A[i] = p[i];
        mu[i] = p[m+i];
        k[i] = p[2*m+i];
        denom[i] = 2.0 * M_PI * bessi0( k[i] );
    }
    double * q = calloc( *xlength, sizeof( double ) );
    if( !q ) {
        error( "can not allocate memory" );
    }
    double * h = calloc( *xlength, sizeof( double ) );
    if( !h ) {
        error( "can not allocate memory" );
    }
    int run = 1;
    *steps = 0;
    while( run == 1 ) {
        int is_nan = 0;
        for( i = 0; i < *xlength; i++ ) {
            q[i] = 0.0;
            for( j = 0; j < m; j++ ) {
                double diffcos = cos( deg2rad( x[i] - mu[j] ) );
                q[i] = q[i] + A[j] * exp( k[j] * diffcos ) / denom[j];
            }
        }
        double tp[3*m];
        for( j = 0; j < m; j++ ) {
            double sumh = 0.0;
            double sumsin = 0.0;
            double sumcos = 0.0;
            for( i = 0; i < *xlength; i++ ) {
                double diffcos = cos( deg2rad( x[i] - mu[j] ) );
                h[i] = A[j] * exp( k[j] * diffcos ) / ( denom[j] * q[i] );
                sumh = sumh + h[i];
                sumsin = sumsin + sin( deg2rad(x[i]) ) * h[i];
                sumcos = sumcos + cos( deg2rad(x[i]) ) * h[i];
            }
            tp[j]   = sumh / *xlength;
            tp[m+j] = rad2deg( atan2( sumsin, sumcos ) );
            double Rbar = sqrt( sumsin*sumsin + sumcos*sumcos ) / sumh;
            tp[2*m+j] = ( 2.0 * Rbar - Rbar*Rbar*Rbar ) / ( 1.0 - Rbar*Rbar );
            denom[j] = 2.0 * M_PI * bessi0( tp[2*m+j] );
            if( *debug > 0 ) {
                printf( "%f %f %f ", tp[j], tp[m+j], tp[2*m+j] );
            }
            if( isnan( tp[2*m+j] ) ) {
                is_nan = 1;
            }
        }
        double llog;
        llvmm( x, xlength, tp, plength, &llog );
        if( is_nan == 1 || llog - prev_llog < *epsilon ) {
            run = 0;
        } else {
            for( j = 0; j < m; j++ ) {
                A[j]  = tp[j];
                mu[j] = tp[m+j];
                k[j]  = tp[2*m+j];
            }
        }
        prev_llog = llog;
        if( *debug > 0 ) {
            printf( "\n" );
        }
        *steps = *steps + 1;
    }
    for( i = 0; i < m; i++ ) {
        ret[i] = A[i];
        ret[m+i] = mu[i];
        ret[2*m+i] = k[i];
    }
    free( q );
    free( h );
}

void vmm_init_vector( int *m, double *ret )
{
    int i;
    for( i = 0; i < *m; i++ ) {
        ret[i] = 1.0/(*m);
        ret[*m+i] = 360/(*m) * i;
        ret[2*(*m)+i] = (((double)*m)/(12*180))*(((double)*m)/(12*180));
    }
}

void polyroot_NR( double *p, int *plength,
                  double *init, double *epsilon,
                  int *debug, double *ret )
{
    double x = *init;
    int steps = 0;

    double d[*plength-1];
    int i;
    for( i = 0; i < *plength-1; i++ ) {
        d[i] = p[i+1] * (i+1);
    }

    int run = 1;
    while( run == 1 ) {
        double powers[*plength];
        powers[0] = 1;
        for( i = 1; i < *plength; i++ ) {
            powers[i] = powers[i-1] * x;
        }

        double numerator = 0.0;
        double denominator = 0.0;
        for( i = 0; i < *plength; i++ ) {
            numerator = numerator + p[i] * powers[i];
            if( i < *plength-1 ) {
                denominator = denominator + d[i] * powers[i];
            }
        }

        double diff = numerator / denominator;
        x = x - diff;
        steps++;
        if( fabs( diff ) < *epsilon ) {
            run = 0;
        }
    }

    if( *debug > 0 ) {
        printf( "Convergence reached after %u iteration(s)\n", steps );
    }

    *ret = x;
}
