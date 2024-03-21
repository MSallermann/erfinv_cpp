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

constexpr long double pi      = 3.1415926535897932384626433832795028841971693993751L;
constexpr long double sqrt_pi = 1.7724538509055160272981674833411451827975494561224L;

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

// Implementation adapted from https://github.com/lakshayg/erfinv same as used in golang math library
template<typename T>
T erfinv_golang_math( T x )
{
    if( x < -1 || x > 1 )
    {
        return std::numeric_limits<T>::quiet_NaN();
    }
    else if( x == 1.0 )
    {
        return std::numeric_limits<T>::infinity();
    }
    else if( x == -1.0 )
    {
        return -std::numeric_limits<T>::infinity();
    }

    const T LN2 = 6.931471805599453094172321214581e-1L;

    const T A0 = 1.1975323115670912564578e0L;
    const T A1 = 4.7072688112383978012285e1L;
    const T A2 = 6.9706266534389598238465e2L;
    const T A3 = 4.8548868893843886794648e3L;
    const T A4 = 1.6235862515167575384252e4L;
    const T A5 = 2.3782041382114385731252e4L;
    const T A6 = 1.1819493347062294404278e4L;
    const T A7 = 8.8709406962545514830200e2L;

    const T B0 = 1.0000000000000000000e0L;
    const T B1 = 4.2313330701600911252e1L;
    const T B2 = 6.8718700749205790830e2L;
    const T B3 = 5.3941960214247511077e3L;
    const T B4 = 2.1213794301586595867e4L;
    const T B5 = 3.9307895800092710610e4L;
    const T B6 = 2.8729085735721942674e4L;
    const T B7 = 5.2264952788528545610e3L;

    const T C0 = 1.42343711074968357734e0L;
    const T C1 = 4.63033784615654529590e0L;
    const T C2 = 5.76949722146069140550e0L;
    const T C3 = 3.64784832476320460504e0L;
    const T C4 = 1.27045825245236838258e0L;
    const T C5 = 2.41780725177450611770e-1L;
    const T C6 = 2.27238449892691845833e-2L;
    const T C7 = 7.74545014278341407640e-4L;

    const T D0 = 1.4142135623730950488016887e0L;
    const T D1 = 2.9036514445419946173133295e0L;
    const T D2 = 2.3707661626024532365971225e0L;
    const T D3 = 9.7547832001787427186894837e-1L;
    const T D4 = 2.0945065210512749128288442e-1L;
    const T D5 = 2.1494160384252876777097297e-2L;
    const T D6 = 7.7441459065157709165577218e-4L;
    const T D7 = 1.4859850019840355905497876e-9L;

    const T E0 = 6.65790464350110377720e0L;
    const T E1 = 5.46378491116411436990e0L;
    const T E2 = 1.78482653991729133580e0L;
    const T E3 = 2.96560571828504891230e-1L;
    const T E4 = 2.65321895265761230930e-2L;
    const T E5 = 1.24266094738807843860e-3L;
    const T E6 = 2.71155556874348757815e-5L;
    const T E7 = 2.01033439929228813265e-7L;

    const T F0 = 1.414213562373095048801689e0L;
    const T F1 = 8.482908416595164588112026e-1L;
    const T F2 = 1.936480946950659106176712e-1L;
    const T F3 = 2.103693768272068968719679e-2L;
    const T F4 = 1.112800997078859844711555e-3L;
    const T F5 = 2.611088405080593625138020e-5L;
    const T F6 = 2.010321207683943062279931e-7L;
    const T F7 = 2.891024605872965461538222e-15L;

    T abs_x = std::abs( x );

    T r, num, den;

    if( abs_x <= 0.85 )
    {
        r   = 0.180625 - 0.25 * x * x;
        num = ( ( ( ( ( ( ( A7 * r + A6 ) * r + A5 ) * r + A4 ) * r + A3 ) * r + A2 ) * r + A1 ) * r + A0 );
        den = ( ( ( ( ( ( ( B7 * r + B6 ) * r + B5 ) * r + B4 ) * r + B3 ) * r + B2 ) * r + B1 ) * r + B0 );
        return x * num / den;
    }

    r = std::sqrt( LN2 - std::log1p( -abs_x ) );
    if( r <= 5.0 )
    {
        r   = r - 1.6;
        num = ( ( ( ( ( ( ( C7 * r + C6 ) * r + C5 ) * r + C4 ) * r + C3 ) * r + C2 ) * r + C1 ) * r + C0 );
        den = ( ( ( ( ( ( ( D7 * r + D6 ) * r + D5 ) * r + D4 ) * r + D3 ) * r + D2 ) * r + D1 ) * r + D0 );
    }
    else
    {
        r   = r - 5.0;
        num = ( ( ( ( ( ( ( E7 * r + E6 ) * r + E5 ) * r + E4 ) * r + E3 ) * r + E2 ) * r + E1 ) * r + E0 );
        den = ( ( ( ( ( ( ( F7 * r + F6 ) * r + F5 ) * r + F4 ) * r + F3 ) * r + F2 ) * r + F1 ) * r + F0 );
    }

    return std::copysign<T>( num / den, x );
}

// Compute the inverse of the error function erf.
// x has to be a real value on the interval [-1,1]
template<typename T>
constexpr T erfinv( T x )
{
    if( x > 1.0 || x < -1.0 )
        return std::numeric_limits<T>::quiet_NaN();

    if( x == 1.0 )
        return std::numeric_limits<T>::infinity();

    if( x == -1.0 )
        return -std::numeric_limits<T>::infinity();

    T res{};

    // for |x| < 0.6 we use the taylor series
    if( x < 0.55 && x > -0.55 )
    {
        res = detail::erfinv_series<T, 10>( x );

        // for |x| < 0.15 the taylor series alone is accurate enough
        if( x < 0.125 && x > -0.125 )
        {
            return res;
        }
        else // else we apply one iteration of halleys method
        {
            return detail::erfinv_halley<T, 1>( x, res );
        }
    }

    // fo |x| >= 0.6 we use winitzkis approximation with two iteration of halleys method
    res = detail::erfinv_winitzki( x );
    res = detail::erfinv_halley<T, 2>( x, res );
    return res;
}

}; // namespace detail

using detail::erfinv;

} // namespace erfinv
