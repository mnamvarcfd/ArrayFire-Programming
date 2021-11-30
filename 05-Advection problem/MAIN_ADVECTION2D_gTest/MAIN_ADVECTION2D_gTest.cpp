#include "gtest/gtest.h"

int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);

    //::testing::GTEST_FLAG(filter) = "-Validation.*";

    int result = RUN_ALL_TESTS();

    return result;
}
