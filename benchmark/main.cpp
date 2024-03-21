#include "erfinv.hpp"
#include "nanobench.h"
#include <iostream>

// NOLINTNEXTLINE

int main()
{
    double x       = 0.5;
    double erfinvx = erfinv::erfinv( x );

    auto bench = ankerl::nanobench::Bench().minEpochIterations( 1000000 );

    bench.run( "erfinv_winitzki", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_winitzki( x ) );
    } );

    bench.run( "erfinv_series_1", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<double, 1>( x ) );
    } );

    bench.run( "erfinv_series_2", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<double, 2>( x ) );
    } );

    bench.run( "erfinv_series_3", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<double, 3>( x ) );
    } );

    bench.run( "erfinv_series_4", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<double, 4>( x ) );
    } );

    bench.run( "erfinv_series_5", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<double, 5>( x ) );
    } );

    bench.run( "erfinv_series_10", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<double, 10>( x ) );
    } );

    bench.run( "erfinv_newton_1", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_newton<double, 1>( x, erfinvx ) );
    } );

    bench.run( "erfinv_newton_2", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_newton<double, 2>( x, erfinvx ) );
    } );

    bench.run( "erfinv_newton_3", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_newton<double, 3>( x, erfinvx ) );
    } );

    bench.run( "erfinv_halley_1", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_halley<double, 1>( x, erfinvx ) );
    } );

    bench.run( "erfinv_halley_2", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_halley<double, 2>( x, erfinvx ) );
    } );

    bench.run( "erfinv_halley_3", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_halley<double, 3>( x, erfinvx ) );
    } );

    bench.run( "erfinv", [&]() { ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::erfinv( x ) ); } );

    double erferfinvx = 0;
    bench.run( "erf", [&]() { ankerl::nanobench::doNotOptimizeAway( erferfinvx = std::erf( erfinvx ) ); } );

    // bench.render( ankerl::nanobench::templates::csv(), std::cout );
    // bench.render( ankerl::nanobench::templates::htmlBoxplot(), std::cout );
}