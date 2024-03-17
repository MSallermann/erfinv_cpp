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

    std::cout << "erfinv(x)       " << inv << "\n";
    std::cout << "x               " << x << "\n";
    std::cout << "erff(erfinv(x)) " << erf_inv << "\n";
    std::cout << "n_digits        " << n_digits << "\n";
}

int main()
{
    std::ofstream f_series( "series.txt" );
    std::ofstream f_winitzki( "winitzki.txt" );
    std::ofstream f_full( "full.txt" );

    int N_Samples = 100;
    for( int i = 0; i < N_Samples; i++ )
    {
        double x = -1.0 + 2.0 * i / ( N_Samples - 1.0 );
        report( x, erfinv::detail::erfinv_series<double, 5>, f_series );
        report( x, erfinv::detail::erfinv_winitzki<double>, f_winitzki );
        report( x, erfinv::erfinv<double>, f_full );
    }

    f_series.close();
    f_winitzki.close();
    f_full.close();
}