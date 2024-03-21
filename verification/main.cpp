#include "cmath"
#include "erfinv.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>

using scalar = double;

template<typename erfinv_func_t>
void report( scalar x, erfinv_func_t erfinv, std::ofstream & fstream )
{
    scalar inv      = erfinv( x );
    scalar erf_inv  = std::erfl( inv );
    scalar rel_err  = std::max( std::abs( ( erf_inv - x ) / x ), std::numeric_limits<scalar>::epsilon() );
    scalar n_digits = -std::log10( rel_err );

    fstream << std::setprecision( 16 );
    fstream << x << " " << inv << " " << erf_inv << " " << n_digits << "\n";
}

template<typename erfinv_func_t>
void run_function_on_samples( const std::string & filename, erfinv_func_t erfinv )
{
    std::ofstream file_out( filename );
    int N_Samples = 10000;
    for( int i = 0; i < N_Samples; i++ )
    {
        scalar x = -1.0 + 2.0 * i / ( N_Samples - 1.0 );
        report( x, erfinv, file_out );
    }
    file_out.close();
}

int main()
{
    run_function_on_samples( "series_5.txt", erfinv::detail::erfinv_series<scalar, 5> );
    run_function_on_samples( "series_10.txt", erfinv::detail::erfinv_series<scalar, 10> );
    run_function_on_samples( "winitzki.txt", erfinv::detail::erfinv_winitzki<scalar> );
    run_function_on_samples( "full.txt", erfinv::erfinv<scalar> );
    run_function_on_samples( "erfinv_golang.txt", erfinv::detail::erfinv_golang_math<scalar> );

    using erfinv::detail::erfinv_series;
    using erfinv::detail::erfinv_winitzki;

    run_function_on_samples( "win_2halley.txt", []( scalar x ) {
        return erfinv::detail::erfinv_halley<scalar, 2>( x, erfinv_winitzki( x ) );
    } );

    run_function_on_samples( "series_10_1halley.txt", []( scalar x ) {
        return erfinv::detail::erfinv_halley<scalar, 1>( x, erfinv_series<scalar, 10>( x ) );
    } );

    run_function_on_samples( "series_10_2halley.txt", []( scalar x ) {
        return erfinv::detail::erfinv_halley<scalar, 2>( x, erfinv_series<scalar, 10>( x ) );
    } );
}