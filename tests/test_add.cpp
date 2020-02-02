#include <gtest/gtest.h>
#include <voltlbx/hello.h>

TEST(Hello, TestAdd)
{    
    ASSERT_EQ(voltlbx::add(3,7), 10);
}
