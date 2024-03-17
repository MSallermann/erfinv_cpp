#pragma once
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>

namespace erfinv
{
namespace detail
{

constexpr double pi = 3.14159265358979323846;

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

template<typename T, size_t... Indices>
constexpr auto create_ck_array( std::index_sequence<Indices...> )
{
    // Use the index sequence to call computeValue for each index
    return std::array<T, sizeof...( Indices )>{ c_k<T>( Indices )... };
}

template<typename T, int n_terms = 5>
T erfinv_series( T x )
{
    constexpr auto coefficient_array = create_ck_array<T>( std::make_index_sequence<n_terms>{} );

    T res     = 0;
    T x_pow_k = x;
    for( int k = 0; k < n_terms; k++ )
    {
        res += std::pow( std::sqrt( pi ) / 2.0 * x, 2 * k + 1 ) / ( 2.0 * k + 1.0 ) * coefficient_array[k];
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
T erf_deriv( T x )
{
    constexpr T pref = 2.0 / std::sqrt( pi );
    return ( pref * std::exp( -x * x ) );
}

template<typename T, int n_steps>
T erfinv_newton( T x, T y )
{
    // find the root of the function f(y) = x - erf( y )
    for( size_t i = 0; i < n_steps; i++ )
    {
        auto f = x - std::erf( y );
        auto df_dy
            = -erf_deriv( y ) + std::numeric_limits<T>::epsilon(); // guard agains divide by zero, by adding epsilon
        y -= f / df_dy;
    }
    return y;
}
}; // namespace detail

// Compute the inverse of the error function erf.
// x has to be a real value on the interval [-1,1]
template<typename T>
T erfinv( T x )
{
    assert( x >= -1.0 );
    assert( x <= 1.0 );

    if( x == 1.0 )
        return std::numeric_limits<T>::infinity();

    if( x == -1.0 )
        return -std::numeric_limits<T>::infinity();

    T initial_approx{};
    if( std::abs( x ) < 0.5 )
        initial_approx = detail::erfinv_series<T, 3>( x );
    else
        initial_approx = detail::erfinv_winitzki<T>( x );

    return detail::erfinv_newton<T, 3>( x, initial_approx );
}

} // namespace erfinv
