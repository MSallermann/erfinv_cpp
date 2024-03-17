#include "cmath"
#include "erfinv.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE( "Test erfinv", "[erfinv]" )
{
    using namespace Catch::Matchers;

    int N_Samples = 100;
    for( int i = 0; i < N_Samples; i++ )
    {
        double x = -1.0 + 2.0 * i / ( N_Samples - 1.0 );

        auto inv     = erfinv::erfinv( x );
        auto erf_inv = std::erf( inv );

        REQUIRE_THAT( x, WithinRel( erf_inv ) );
    }
}