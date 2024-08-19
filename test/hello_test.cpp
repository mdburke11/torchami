#include <gtest/gtest.h>

#include "../src/tami_base_src/tami_base.hpp"

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
  // Expect two strings not to be equal.
  EXPECT_STRNE("hello", "world");
  // Expect equality.
  EXPECT_EQ(7 * 6, 42);
}

TEST(HelloTest, second){

std::cout<<"In test"<<std::endl;
ASSERT_NEAR(0.0,0.0,1e-13);

}