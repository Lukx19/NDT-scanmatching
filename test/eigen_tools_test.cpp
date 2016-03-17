// Bring in my package's API, which is what I'm testing
#include <eigen_tools.h>
// Bring in gtest
#include <gtest/gtest.h>

using namespace eigt;

double EPSILON = 0.001;

TEST(EigenTools, createTransFromPosesSimple)
{
  try {
    pose2d_t<double> first_pose;
    first_pose << 0, 0, 0;
    pose2d_t<double> second_pose;
    second_pose << 10, 10, 0;
    transform2d_t<double> res_trans;
    res_trans.matrix() << 1,0,10,0,1,10,0,0,1;
    auto calc = transBtwPoses(first_pose,second_pose).matrix();
    EXPECT_TRUE ((res_trans.matrix() - calc).cwiseAbs().sum() < EPSILON);
    //ADD_FAILURE() <<res_trans.matrix()<< "\n\n"<< calc;
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

TEST(EigenTools, createTransFromPosesAdvance)
{
  try {
    pose2d_t<double> first_pose;
    first_pose << 10, 10, 0;
    pose2d_t<double> second_pose;
    second_pose << 20, 20, 0;
    transform2d_t<double> res_trans;
    res_trans.matrix() << 1,0,10,0,1,10,0,0,1;
    auto calc = transBtwPoses(first_pose,second_pose).matrix();
    EXPECT_TRUE ((res_trans.matrix() - calc).cwiseAbs().sum() < EPSILON);
    //ADD_FAILURE() <<res_trans.matrix()<< "\n\n"<< calc;
  } catch (...) {
    ADD_FAILURE() << "Uncaught exception";
  }
}

// Run all the tests that were declared with TEST()
int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
