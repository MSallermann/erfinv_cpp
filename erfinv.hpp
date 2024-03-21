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
constexpr double pi_3_2  = 5.56832799683170784528;

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
    for( int i = 0; i < n_steps; i++ )
    {
        const T y2       = y * y;
        const T erfy     = std::erf( y );
        const T exp_m_y2 = std::exp( -y2 );
        const T f        = x - erfy;
        const T df_dy    = -2.0 / sqrt_pi * exp_m_y2;
        y -= f / df_dy;
    }
    return y;
}

template<typename T, int n_steps>
T erfinv_halley( T x, T y )
{
    // find the root of the function f(y) = x - erf( y )
    for( int i = 0; i < n_steps; i++ )
    {
        const T y2       = y * y;
        const T erfy     = std::erf( y );
        const T exp_m_y2 = std::exp( -y2 );
        const T f        = x - erfy;
        const T df_dy    = -2.0 / sqrt_pi * exp_m_y2;
        const T df2_d2y  = 4.0 * exp_m_y2 * y / sqrt_pi;

        y -= 2.0 * f * df_dy / ( 2 * df_dy * df_dy - f * df2_d2y );
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

    // for |x| < 0.65 we use the taylor series
    if( x < 0.65 && x > -0.65 )
    {
        res = detail::erfinv_series<T, 10>( x );

        // for |x| < 0.2 the taylor series alone is accurate enough
        if( x < 0.2 && x > -0.2 )
        {
            return res;
        }
        else // else we apply one iteration of halleys method
        {
            return detail::erfinv_halley<T, 1>( x, res );
        }
    }

    // fo |x| >= 0.65 we use winitzkis approximation with two iteration of halleys method
    res = detail::erfinv_winitzki( x );
    res = detail::erfinv_halley<T, 2>( x, res );
    return res;
}

}; // namespace detail

using detail::erfinv;

} // namespace erfinv
