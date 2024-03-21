#include "erfinv.hpp"
#include "nanobench.h"

// NOLINTNEXTLINE

int main()
{
    using scalar = double;

    scalar x       = 0.5;
    scalar erfinvx = erfinv::erfinv( x );

    auto bench = ankerl::nanobench::Bench().minEpochIterations( 1000000 );

    bench.run( "erfinv_winitzki", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_winitzki( x ) );
    } );

    bench.run( "erfinv_series_1", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<scalar, 1>( x ) );
    } );

    bench.run( "erfinv_series_2", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<scalar, 2>( x ) );
    } );

    bench.run( "erfinv_series_3", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<scalar, 3>( x ) );
    } );

    bench.run( "erfinv_series_4", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<scalar, 4>( x ) );
    } );

    bench.run( "erfinv_series_5", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<scalar, 5>( x ) );
    } );

    bench.run( "erfinv_series_10", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<scalar, 10>( x ) );
    } );

    bench.run( "erfinv_newton_1", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_newton<scalar, 1>( x, erfinvx ) );
    } );

    bench.run( "erfinv_newton_2", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_newton<scalar, 2>( x, erfinvx ) );
    } );

    bench.run( "erfinv_newton_3", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_newton<scalar, 3>( x, erfinvx ) );
    } );

    bench.run( "erfinv_halley_1", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_halley<scalar, 1>( x, erfinvx ) );
    } );

    bench.run( "erfinv_halley_2", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_halley<scalar, 2>( x, erfinvx ) );
    } );

    bench.run( "erfinv_halley_3", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_halley<scalar, 3>( x, erfinvx ) );
    } );

    bench.run( "erfinv", [&]() { ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::erfinv( x ) ); } );

    bench.run( "erfinv_golang", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_golang_math( x ) );
    } );

    scalar erferfinvx = 0;
    bench.run( "erfl", [&]() { ankerl::nanobench::doNotOptimizeAway( erferfinvx = std::erfl( erfinvx ) ); } );

    // bench.render( ankerl::nanobench::templates::csv(), std::cout );
    // bench.render( ankerl::nanobench::templates::htmlBoxplot(), std::cout );
}