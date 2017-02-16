#include <fstream>
#include <cppunit/TestCase.h>
#include <cppunit/TestFixture.h>
#include <cppunit/ui/text/TextTestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/XmlOutputter.h>
#include "../Kernel.h"
#include "../../extern/eigen/include/Eigen/Core"
using namespace CppUnit;
using namespace std;

class KernelTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(KernelTest);
    CPPUNIT_TEST(testPoly6_isZeroifVectorLengthGreaterH);
    CPPUNIT_TEST(testPoly6_isZeroPoint0826179682_ifVectorLengthSmallerH);
    CPPUNIT_TEST(testSpiky_isZeroifVectorLengthGreaterH);
    CPPUNIT_TEST(testSpiky_isZeroPoint0024446199_ifVectorLengthSmallerH);
    CPPUNIT_TEST(testCubic_isZeroifVectorLengthDividedByHGreater2);
    CPPUNIT_TEST(testCubic_isZeroPoint4791666666666ifVectorLengthDividedByHGreater0Smaller1);
    CPPUNIT_TEST(testCubic_isZeroPoint02083333333ifVectorLengthDividedByHGreater1Smaller2);
    CPPUNIT_TEST_SUITE_END();

public:
    void setUp(void);
    void tearDown(void);

protected:
    void testPoly6_isZeroifVectorLengthGreaterH(void);
    void testPoly6_isZeroPoint0826179682_ifVectorLengthSmallerH(void);
    void testSpiky_isZeroifVectorLengthGreaterH(void);
    void testSpiky_isZeroPoint0024446199_ifVectorLengthSmallerH(void);
    void testCubic_isZeroifVectorLengthDividedByHGreater2(void);
	void testCubic_isZeroPoint4791666666666ifVectorLengthDividedByHGreater0Smaller1(void);
	void testCubic_isZeroPoint02083333333ifVectorLengthDividedByHGreater1Smaller2(void);
private:

    Kernel *mTestObj;
	
	
};


void KernelTest::testPoly6_isZeroifVectorLengthGreaterH(void)
{
    CPPUNIT_ASSERT_EQUAL(0.0, mTestObj->w_poly6(Eigen::Vector3d(1.0, 1.0, 1.0), 1.0));
}

void KernelTest::testPoly6_isZeroPoint0826179682_ifVectorLengthSmallerH(void)
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0826179682, mTestObj->w_poly6(Eigen::Vector3d(1.0, 0.0, 0.0), 2.0), 0.000001);
}

void KernelTest::testSpiky_isZeroifVectorLengthGreaterH(void)
{
    CPPUNIT_ASSERT_EQUAL(0.0, mTestObj->w_spiky(Eigen::Vector3d(1.0, 1.0, 1.0), 1.0));
}

void KernelTest::testSpiky_isZeroPoint0024446199_ifVectorLengthSmallerH(void)
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0024446199, mTestObj->w_spiky(Eigen::Vector3d(2.0, 2.0, 1.0), 5.0), 0.000001);
}

void KernelTest::testCubic_isZeroifVectorLengthDividedByHGreater2(void)
{
    CPPUNIT_ASSERT_EQUAL(0.0, mTestObj->w_cubicSpline(Eigen::Vector3d(2.0, 2.0, 1.0), 1.0));
}

void KernelTest::testCubic_isZeroPoint4791666666666ifVectorLengthDividedByHGreater0Smaller1(void)
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.4791666666666, mTestObj->w_cubicSpline(Eigen::Vector3d(0.0, 1.0, 0.0), 2.0), 0.000001);
}

void KernelTest::testCubic_isZeroPoint02083333333ifVectorLengthDividedByHGreater1Smaller2(void)
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333, mTestObj->w_cubicSpline(Eigen::Vector3d(2.0, 1.0, 2.0), 2.0), 0.000001);
}

void KernelTest::setUp(void)
{
    mTestObj = new Kernel();
}

void KernelTest::tearDown(void)
{
    delete mTestObj;
}

//-----------------------------------------------------------------------------

CPPUNIT_TEST_SUITE_REGISTRATION( KernelTest );

int main(int argc, char* argv[])
{
    // informs test-listener about testresults
    CPPUNIT_NS::TestResult testresult;

    // register listener for collecting the test-results
    CPPUNIT_NS::TestResultCollector collectedresults;
    testresult.addListener (&collectedresults);

    // register listener for per-test progress output
    CPPUNIT_NS::BriefTestProgressListener progress;
    testresult.addListener (&progress);

    // insert test-suite at test-runner by registry
    CPPUNIT_NS::TestRunner testrunner;
    testrunner.addTest (CPPUNIT_NS::TestFactoryRegistry::getRegistry().makeTest ());
    testrunner.run(testresult);

    // output results in compiler-format
    CPPUNIT_NS::CompilerOutputter compileroutputter(&collectedresults, std::cerr);
    compileroutputter.write ();

    // Output XML for Jenkins CPPunit plugin
    ofstream xmlFileOut("cppTestBasicMathResults.xml");
    XmlOutputter xmlOut(&collectedresults, xmlFileOut);
    xmlOut.write();

    // return 0 if tests were successful
    return collectedresults.wasSuccessful() ? 0 : 1;
}
