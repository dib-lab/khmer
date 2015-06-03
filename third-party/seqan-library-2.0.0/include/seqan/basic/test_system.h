// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Alternative, more flexible test system, based on the ideas and mechanisms
// used by google-test system.
//
// See http://code.google.com/p/googletest/ for the better implementation. ;)
// ==========================================================================

// TODO(holtgrew): Add fixture support.
// TODO(holtgrew): Add DDDoc documentation.
// TODO(holtgrew): Port all tests over to new system?

#ifndef INCLUDE_SEQAN_BASIC_TEST_SYSTEM_H_
#define INCLUDE_SEQAN_BASIC_TEST_SYSTEM_H_

#ifdef PLATFORM_WINDOWS
#include <typeinfo>
#endif  // #ifdef PLATFORM_WINDOWS

#include <seqan/basic/fundamental_tags.h>

#include <memory>
#include <string>

#ifdef PLATFORM_GCC
#include <cxxabi.h>
#endif  // #ifdef PLATFORM_GCC

namespace seqan {

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class Test
// --------------------------------------------------------------------------

// Base class for tests in the test system.

class Test
{
public:
    virtual void setUp() {}
    virtual void runTest() = 0;
    virtual void tearDown() {}

    virtual ~Test() {}
};

// --------------------------------------------------------------------------
// Class TestDescription_
// --------------------------------------------------------------------------

// Stores information about a test.

class TestDescription_
{
public:
    std::string testCaseName;
    std::string testName;
    std::string typeName;
    std::SEQAN_AUTO_PTR_NAME<Test> instance;

    TestDescription_(char const * testCaseName, char const * testName,
                     char const * typeName, Test * instance) :
            testCaseName(testCaseName), testName(testName), typeName(typeName), instance(instance)
    {}
};

// --------------------------------------------------------------------------
// Class TestSystem
// --------------------------------------------------------------------------

// Registry for the test system.
//
// A main function using this system looks as follows:
//
// int main(int argc, char const ** argv)
// {
//     seqan::TestSystem::init(argc, argv);
//     return seqan::TestSystem::runAll();
// }

class TestSystem
{
public:
    std::vector<TestDescription_ *> testDescriptions;

    TestSystem()
    {}

    static void init(int argc, char const ** argv)
    {
        (void)argc;
        seqan::ClassTest::beginTestSuite("tests", argv[0]);
    }

    static TestSystem * getInstance()
    {
        static TestSystem instance;
        return &instance;
    }

    void registerTest(TestDescription_ * description)
    {
        testDescriptions.push_back(description);
    }

    static int runAll()
    {
        TestSystem & instance = *getInstance();
        typedef std::vector<TestDescription_ *>::iterator TIt;
        for (TIt it = instance.testDescriptions.begin(); it != instance.testDescriptions.end(); ++it)
        {
            std::string testName = (*it)->testCaseName;
            testName += "_";
            testName += (*it)->testName;
            if (!(*it)->typeName.empty())
            {
                testName += " type parameter ";
                testName += (*it)->typeName;
            }
            seqan::ClassTest::beginTest(testName.c_str());
            try {
                (*it)->instance->setUp();
                (*it)->instance->runTest();
                (*it)->instance->tearDown();
            } catch(seqan::ClassTest::AssertionFailedException e) {
                /* Swallow exception, go on with next test. */
                (void) e;  /* Get rid of unused variable warning. */
            } catch (seqan::Exception const & e) {
                std::cerr << "Unexpected exception of type "
                            << toCString(seqan::Demangler< seqan::Exception>(e))
                            << "; message: " << e.what() << "\n";
                seqan::ClassTest::StaticData::thisTestOk() = false;
                seqan::ClassTest::StaticData::errorCount() += 1;
            } catch (...) {
                std::cerr << "Unexpected exception of unknown type\n";
                seqan::ClassTest::StaticData::thisTestOk() = false;
                seqan::ClassTest::StaticData::errorCount() += 1;
            }
            seqan::ClassTest::endTest();
        }
        return seqan::ClassTest::endTestSuite();
    }
};

// --------------------------------------------------------------------------
// Class TestCaseFactory_
// --------------------------------------------------------------------------

// Helper class for creating tests.

template <typename TTest>
class TestCaseFactory_
{
public:
    static TestDescription_ * make(char const * testCaseName, char const * testName,
                                  char const * typeName = "")
    {
        TestSystem * testSystem = TestSystem::getInstance();
        TestDescription_ * desc = new TestDescription_(testCaseName, testName, typeName, new TTest);
        testSystem->testDescriptions.push_back(desc);
        return desc;
    }
};

// --------------------------------------------------------------------------
// Class TypedTestFactory_
// --------------------------------------------------------------------------

// Helper class for creating typed tests.

template <template <typename> class TTestCase, typename TTagList>
class TypedTestFactory_;

template <template <typename> class TTestCase>
class TypedTestFactory_<TTestCase, void>
{
public:
    static bool make(char const *, char const *)
    {
        return true;
    }
};

template <template <typename> class TTestCase, typename TType, typename TSubList>
class TypedTestFactory_<TTestCase, TagList<TType, TSubList> >
{
public:
    // TODO(esiragusa): use Demangler.
    template <typename T>
    static std::string getTypeName()
    {
        const char* const name = typeid(T).name();
#ifdef PLATFORM_GCC
        int status = 0;
        char* const readableName = abi::__cxa_demangle(name, 0, 0, &status);
        std::string nameString(status == 0 ? readableName : name);
        free(readableName);
        return nameString;
#else  // #ifdef PLATFORM_GCC
        return name;
#endif  // #ifdef PLATFORM_GCC
    }

    static bool make(char const * testCaseName, char const * testName)
    {
        TestSystem * testSystem = TestSystem::getInstance();
        std::string typeName = getTypeName<TType>();
        TestDescription_ * desc = new TestDescription_(testCaseName, testName, typeName.c_str(), new TTestCase<TType>);
        testSystem->testDescriptions.push_back(desc);
        return TypedTestFactory_<TTestCase, TSubList>::make(testCaseName, testName);
    }
};

// ==========================================================================
// Metafunctions
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

// ==========================================================================
// Macros
// ==========================================================================

// --------------------------------------------------------------------------
// Helper Macro SEQAN_TEST_NAME_()
// --------------------------------------------------------------------------

// Helper macro for getting the name of a test.

#define SEQAN_TEST_NAME_(testCaseName, testName) \
    testCaseName ## __ ## testName

// --------------------------------------------------------------------------
// Macro SEQAN_TEST()
// --------------------------------------------------------------------------

// Macro for defining a test.

#define SEQAN_TEST(testCaseName, testName)                                    \
    class SEQAN_TEST_NAME_(testCaseName, testName) : public seqan::Test     \
    {                                                                         \
    public:                                                                   \
        SEQAN_TEST_NAME_(testCaseName, testName)() {}                         \
                                                                              \
        virtual void runTest();                                               \
                                                                              \
        static seqan::TestDescription_ * description;                       \
    };                                                                        \
                                                                              \
    seqan::TestDescription_ *                                               \
    SEQAN_TEST_NAME_(testCaseName, testName)::description =                   \
    seqan::TestCaseFactory_<SEQAN_TEST_NAME_(testCaseName, testName)>::make(\
            SEQAN_MKSTRING(testCaseName),                                     \
            SEQAN_MKSTRING(testName));                                        \
                                                                              \
    void SEQAN_TEST_NAME_(testCaseName, testName)::runTest()

// --------------------------------------------------------------------------
// Macro SEQAN_TEST_F()
// --------------------------------------------------------------------------

// Macro for defining a test with a fixture.

#define SEQAN_TEST_F(testCaseName, testName)                                  \
    class SEQAN_TEST_NAME_(testCaseName, testName) : public testCaseName      \
    {                                                                         \
    public:                                                                   \
        SEQAN_TEST_NAME_(testCaseName, testName)() {}                         \
                                                                              \
        virtual void runTest();                                               \
                                                                              \
        static seqan::TestDescription_ * description;                       \
    };                                                                        \
                                                                              \
    seqan::TestDescription_ *                                               \
    SEQAN_TEST_NAME_(testCaseName, testName)::description =                   \
    seqan::TestCaseFactory_<SEQAN_TEST_NAME_(testCaseName, testName)>::make(\
            SEQAN_MKSTRING(testCaseName),                                     \
            SEQAN_MKSTRING(testName));                                        \
                                                                              \
    void SEQAN_TEST_NAME_(testCaseName, testName)::runTest()

// --------------------------------------------------------------------------
// Helper Macro SEQAN_TYPED_TEST_CASE_TYPES_NAME_()
// --------------------------------------------------------------------------

// Helper macro for getting the name of a test case.

#define SEQAN_TYPED_TEST_CASE_TYPES_NAME_(testCaseName, types)  \
    SEQAN_TYPED_TEST_CASE_TYPES_ ## testCaseName ## _

// --------------------------------------------------------------------------
// Macro SEQAN_TYPED_TEST_CASE()
// --------------------------------------------------------------------------

// Helper macro for fixing the types to run with a typed test case.

#define SEQAN_TYPED_TEST_CASE(testCaseName, types)                      \
    typedef types SEQAN_TYPED_TEST_CASE_TYPES_NAME_(testCaseName, types)

// --------------------------------------------------------------------------
// Macro SEQAN_TYPED_TEST()
// --------------------------------------------------------------------------

// Define a typed test, i.e. one that is run for all types in a list.

#define SEQAN_TYPED_TEST(testCaseName, testName)    \
    template <typename SEQAN_TParam> \
    class SEQAN_TEST_NAME_(testCaseName, testName) : public testCaseName<SEQAN_TParam> \
    {                                                                   \
    public:                                                             \
        SEQAN_TEST_NAME_(testCaseName, testName)(){}                    \
                                                                        \
        virtual void runTest();                                         \
                                                                        \
        typedef testCaseName<SEQAN_TParam> TestFixture;                 \
    };                                                                  \
                                                                        \
    bool SEQAN_ ## testCaseName ## __ ## testName ## _registered_ =     \
            seqan::TypedTestFactory_<SEQAN_TEST_NAME_(testCaseName, testName), \
                               SEQAN_TYPED_TEST_CASE_TYPES_NAME_(testCaseName, types) \
                               >::make(SEQAN_MKSTRING(testCaseName), SEQAN_MKSTRING(testName)); \
                                                                        \
    template <typename SEQAN_TParam>                                    \
    void SEQAN_TEST_NAME_(testCaseName, testName)<SEQAN_TParam>::runTest()

// --------------------------------------------------------------------------
// Macro SEQAN_TEST_EXCEPTION()
// --------------------------------------------------------------------------

#define SEQAN_TEST_EXCEPTION(_exception_type, command)                              \
    do                                                                              \
    {                                                                               \
        bool caughtException = false;                                               \
        try                                                                         \
        {                                                                           \
            command;                                                                \
        }                                                                           \
        catch(_exception_type& ex)                                                  \
        {                                                                           \
            caughtException = true;                                                 \
        }                                                                           \
        catch(...)                                                                  \
        {                                                                           \
            SEQAN_FAIL("Wrong exception thrown: %s", #_exception_type);             \
        }                                                                           \
        if (!caughtException)                                                       \
            SEQAN_FAIL("No exception thrown!");                                     \
    } while(false)

#define SEQAN_TEST_EXCEPTION_WITH_MESSAGE(_exception_type, command, _message)       \
    do                                                                              \
    {                                                                               \
        bool caughtException = false;                                               \
        try                                                                         \
        {                                                                           \
            command;                                                                \
        }                                                                           \
        catch(_exception_type& ex)                                                  \
        {                                                                           \
            if(std::string(ex.what()) != _message)                                  \
                SEQAN_FAIL("Got correct exception but wrong message: '%s' != '%s'", \
                           ex.what(), _message);                                    \
            caughtException = true;                                                 \
        }                                                                           \
        catch(...)                                                                  \
        {                                                                           \
            SEQAN_FAIL("Wrong exception thrown: %s", #_exception_type);             \
        }                                                                           \
        if (!caughtException)                                                       \
            SEQAN_FAIL("No exception thrown!");                                     \
    } while(false)


}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BASIC_TEST_SYSTEM_H_
