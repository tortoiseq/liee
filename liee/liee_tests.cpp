#include "my_util.hpp"

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
//#include <boost/test/unit_test_log.hpp>
//#include <boost/filesystem/fstream.hpp>
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( my_utils )
BOOST_AUTO_TEST_CASE( cerf_test )
{
    dcmplx result = liee::cerf( dcmplx(0.0, 0.0), 100 );
    BOOST_CHECK_CLOSE( result.real(), 0, 1e-6);
    BOOST_CHECK_CLOSE( result.imag(), 0, 1e-6);

    result = liee::cerf( dcmplx(0.0, 1.0), 100 );
    BOOST_CHECK_CLOSE( result.real(), 0, 1e-6);
    BOOST_CHECK_CLOSE( result.imag(), 1.65042575880, 1e-6);

    result = liee::cerf( dcmplx(1.0, 1.0), 100 );
    BOOST_CHECK_CLOSE( result.real(), 1.31615128170, 1e-6);
    BOOST_CHECK_CLOSE( result.imag(), 0.190453469238, 1e-6);
}
BOOST_AUTO_TEST_SUITE_END()
