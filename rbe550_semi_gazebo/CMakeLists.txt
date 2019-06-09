cmake_minimum_required(VERSION 2.8.3)
project(rbe550_semi_gazebo)

find_package(catkin REQUIRED COMPONENTS
  rospy
  std_msgs
  geometry_msgs
)

catkin_package()

find_package(Boost REQUIRED COMPONENTS system)
find_package(gazebo REQUIRED)

# install(DIRECTORY launch
#   DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION})
#
# install(DIRECTORY worlds
#   DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION})
