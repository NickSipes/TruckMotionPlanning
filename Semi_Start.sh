#!/bin/bash

# Starts ROS and Gazebo with the Semi-Truck model for RBE-550 Project.

cd ~/catkin_ws/
catkin_make
source devel/setup.bash

roslaunch rbe550_semi_gazebo rbe550_world.launch


