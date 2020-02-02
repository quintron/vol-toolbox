#include <gtest/gtest.h>
#include <voltlbx/hello.h>

TEST(Hello, TestAdd)
{    
    ASSERT_EQ(voltlbx::add(3,7), 10);
}

TEST(Hello, TestEigen)
{
    voltlbx::test_eigen();
    ASSERT_TRUE(true);
}
