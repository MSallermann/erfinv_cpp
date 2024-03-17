#include "erfinv.hpp"
#include "nanobench.h"

// NOLINTNEXTLINE

int main()
{
    double x       = 0.5;
    double erfinvx = erfinv::erfinv( x );

    ankerl::nanobench::Bench().run( "erfinv_winitzki", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_winitzki( x ) );
    } );

    ankerl::nanobench::Bench().run( "erfinv_series_1", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<double, 1>( x ) );
    } );

    ankerl::nanobench::Bench().run( "erfinv_series_2", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<double, 2>( x ) );
    } );

    ankerl::nanobench::Bench().run( "erfinv_series_3", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<double, 2>( x ) );
    } );

    ankerl::nanobench::Bench().run( "erfinv_series_4", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<double, 2>( x ) );
    } );

    ankerl::nanobench::Bench().run( "erfinv_series_5", [&]() {
        ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::detail::erfinv_series<double, 2>( x ) );
    } );

    ankerl::nanobench::Bench().run(
        "erfinv", [&]() { ankerl::nanobench::doNotOptimizeAway( erfinvx = erfinv::erfinv( x ) ); } );

    double erferfinvx = 0;
    ankerl::nanobench::Bench().run(
        "erf", [&]() { ankerl::nanobench::doNotOptimizeAway( erferfinvx = std::erf( erfinvx ) ); } );
}