#pragma once
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
namespace erfinv
{
namespace detail
{

constexpr double pi      = 3.14159265358979323846;
constexpr double sqrt_pi = 1.77245385090551602729;

template<typename T>
constexpr T constexpr_pow( T x, int n )
{
    T res = 1;
    for( int i = 0; i < n; i++ )
    {
        res *= x;
    }
    return res;
}

template<typename T>
constexpr T c_k( int k )
{
    if( k == 0 )
        return 1.0;

    T res = 0.0;
    for( int m = 0; m < k; m++ )
    {
        res += c_k<T>( m ) * c_k<T>( k - 1 - m ) / ( ( m + 1.0 ) * ( 2.0 * m + 1.0 ) );
    }
    return res;
}

template<typename T>
constexpr T pref_k( int k )
{
    return constexpr_pow( sqrt_pi / 2.0, 2 * k + 1 ) / ( 2.0 * k + 1.0 );
}

template<typename T, size_t... Indices>
constexpr auto create_taylor_coeff_array( std::index_sequence<Indices...> )
{
    // Use the index sequence to call computeValue for each index
    return std::array<T, sizeof...( Indices )>{ pref_k<T>( Indices ) * c_k<T>( Indices )... };
}

template<typename T, int n_terms = 5>
constexpr T erfinv_series( T x )
{
    constexpr auto coefficient_array = create_taylor_coeff_array<T>( std::make_index_sequence<n_terms>{} );

    T res     = 0;
    T x_pow_k = x;

    for( int k = 0; k < n_terms; k++ )
    {
        res += coefficient_array[k] * x_pow_k;
        x_pow_k *= x * x;
    }
    return res;
}

template<typename T>
constexpr T sign( T x )
{
    return ( x > 0 ) ? 1.0 : -1.0;
}

template<typename T>
T erfinv_winitzki( T x )
{
    constexpr T a = 8.0 * ( pi - 3.0 ) / ( 3.0 * pi * ( 4.0 - pi ) );
    const T x2    = x * x;

    const T t1   = 2.0 / ( pi * a ) + std::log1p( -x2 ) / 2.0;
    const T t1_2 = t1 * t1;

    const T t2 = std::log1p( -x2 ) / a;
    const T t3 = 2.0 / ( pi * a ) + std::log1p( -x2 ) / 2.0;

    return sign( x ) * std::sqrt( std::sqrt( t1_2 - t2 ) - t3 );
}

template<typename T>
constexpr T erf_deriv( T x )
{
    constexpr T pref = 2.0 / sqrt_pi;
    return pref * std::exp( -x * x );
}

template<typename T, int n_steps>
T erfinv_newton( T x, T y )
{
    // find the root of the function f(y) = x - erf( y )
    for( size_t i = 0; i < n_steps; i++ )
    {
        const auto f     = x - std::erf( y );
        const auto df_dy = -erf_deriv( y ) + std::numeric_limits<T>::epsilon(); // guard against divide by zero
        y -= f / df_dy;
    }
    return y;
}

// Compute the inverse of the error function erf.
// x has to be a real value on the interval [-1,1]
template<typename T>
constexpr T erfinv( T x )
{
    assert( x >= -1.0 );
    assert( x <= 1.0 );

    if( x == 1.0 )
        return std::numeric_limits<T>::infinity();

    if( x == -1.0 )
        return -std::numeric_limits<T>::infinity();

    T res{};

    if( x > 0.7 || x < -0.7 )
    {
        res = detail::erfinv_winitzki( x );
        res = detail::erfinv_newton<T, 3>( x, res );
        return res;
    }

    res = detail::erfinv_series<T, 10>( x );
    if( x < 0.2 && x > -0.2 )
        return res;

    return detail::erfinv_newton<T, 2>( x, res );
}

}; // namespace detail

using detail::erfinv;

} // namespace erfinv
