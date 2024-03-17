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
    auto inv     = erfinv( x );
    auto erf_inv = std::erf( inv );

    auto n_digits = -std::log10( std::abs( ( erf_inv - x ) / x ) + std::numeric_limits<double>::epsilon() );

    fstream << std::setprecision( 16 );
    fstream << x << " " << inv << " " << erf_inv << " " << n_digits << "\n";
}

int main()
{
    std::ofstream f_series_5( "series_5.txt" );
    std::ofstream f_series_10( "series_10.txt" );
    // std::ofstream f_series_15( "series_15.txt" );
    // std::ofstream f_series_20( "series_20.txt" );
    std::ofstream f_winitzki( "winitzki.txt" );
    std::ofstream f_full( "full.txt" );

    int N_Samples = 10000;
    for( int i = 0; i < N_Samples; i++ )
    {
        double x = -1.0 + 2.0 * i / ( N_Samples - 1.0 );
        report( x, erfinv::detail::erfinv_series<double, 5>, f_series_5 );
        report( x, erfinv::detail::erfinv_series<double, 10>, f_series_10 );
        // report( x, erfinv::detail::erfinv_series<double, 15>, f_series_15 );
        // report( x, erfinv::detail::erfinv_series<double, 10>, f_series_20 );
        report( x, erfinv::detail::erfinv_winitzki<double>, f_winitzki );
        report( x, erfinv::erfinv<double>, f_full );
    }

    f_series_5.close();
    f_series_10.close();
    // f_series_15.close();
    // f_series_20.close();
    f_winitzki.close();
    f_full.close();
}