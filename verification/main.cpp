#include "cmath"
#include "erfinv.hpp"
#include "iostream"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>

template<typename erfinv_func_t>
void report( double x, erfinv_func_t erfinv, std::ofstream & fstream )
{
    auto inv      = erfinv( x );
    auto erf_inv  = std::erf( inv );
    auto rel_err  = std::max( std::abs( ( erf_inv - x ) / x ), std::numeric_limits<double>::epsilon() );
    auto n_digits = -std::log10( rel_err );

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
        double x = -1.0 + 2.0 * i / ( N_Samples - 1.0 );
        report( x, erfinv, file_out );
    }
    file_out.close();
}

int main()
{
    run_function_on_samples( "series_5.txt", erfinv::detail::erfinv_series<double, 5> );
    run_function_on_samples( "series_10.txt", erfinv::detail::erfinv_series<double, 10> );
    run_function_on_samples( "winitzki.txt", erfinv::detail::erfinv_winitzki<double> );
    run_function_on_samples( "full.txt", erfinv::erfinv<double> );

    using erfinv::detail::erfinv_series;
    using erfinv::detail::erfinv_winitzki;

    run_function_on_samples( "win_newton.txt", []( double x ) {
        return erfinv::detail::erfinv_newton<double, 1>( x, erfinv_winitzki( x ) );
    } );

    run_function_on_samples( "win_2halley.txt", []( double x ) {
        return erfinv::detail::erfinv_halley<double, 2>( x, erfinv_winitzki( x ) );
    } );

    run_function_on_samples( "series_10_1halley.txt", []( double x ) {
        return erfinv::detail::erfinv_halley<double, 1>( x, erfinv_series<double, 10>( x ) );
    } );

    run_function_on_samples( "series_10_2halley.txt", []( double x ) {
        return erfinv::detail::erfinv_halley<double, 2>( x, erfinv_series<double, 10>( x ) );
    } );
}