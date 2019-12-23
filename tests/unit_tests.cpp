#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "../src/vars.h"
#include "../src/derivs.h"
//~ #include "../src/fftw_deriv.h"
#include "../src/ints.h"
#include "../src/masks.h"

#include "../src/functions.h"
#include "../src/compNACV.h"
#include "../src/linsol.h"
#include "../src/state_objs.h"
#include "../src/RK45.h"
#include "../src/io.h"

#include "headers/catch.hpp"


TEST_CASE( "density verification", "[density]" ) {
	// check that having a zero resource semaphore times out
	auto nd0 = FileRead1D_dcomp( "tests/verify/149nd.txt" );
	auto nd1 = FileRead1D_dcomp( "ShinMetiu/output_verify/149nd.txt" );
    REQUIRE(nd0 == nd1);
}
	
TEST_CASE( "tdpes verification", "[tdpes]" ) {
	// check that having a zero resource semaphore times out
	auto tdpes0 = FileRead1D_dcomp( "tests/verify/149TDPES.txt" );
	auto tdpes1 = FileRead1D_dcomp( "ShinMetiu/output_verify/149TDPES.txt" );
    REQUIRE(tdpes0 == tdpes1);
}
	
TEST_CASE( "C1 verification", "[c1]" ) {
	// check that having a zero resource semaphore times out
	auto c10 = FileRead1D_dcomp( "tests/verify/149C1_sq.txt" );
	auto c11 = FileRead1D_dcomp( "ShinMetiu/output_verify/149C1_sq.txt" );
    REQUIRE(c10 == c11);
}
