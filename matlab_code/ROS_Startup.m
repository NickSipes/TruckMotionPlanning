%% ROS Setup
% This script sets up the environment and variables necessary in the
% workspace to allow other scripts and functions to interface with the ROS
% network and Gazebo simulation.

% Y if MATLAB and ROS/Gazebo in the same Linux Environment, N if not
environ_same = "N";

% Set up Environment Variables
if environ_same ~= "Y" && ~exist('ip_ros','var')
    ip_ros = input('Enter the inet addr from the virtual machine','s');
    ip_matlab = input('Enter the IPv4 Address from the MATLAB machine','s');
    setenv('ROS_MASTER_URI',strcat('http://',ip_ros,':11311'));
    setenv('ROS_IP',ip_matlab);
end

% Start ROS Node for MATLAB
rosinit(ip_ros)

%% Variables

% Global Variables for joint state and semitruck position data

global joint_states joint_states_time semitruck_position

joint_states_time  = zeros(2,1);
semitruck_position = zeros(3,1);

% Truck Parameters
L1 = 5.3273;   % Tractor Wheelbase (m)
L2 = 9.46755;  % Trailer Wheelbase (m)

%% Publishers

% Joint Position Publishers

left_steer_pub = rospublisher('/semitruck/left_steer_position_controller/command');
right_steer_pub = rospublisher('/semitruck/right_steer_position_controller/command');
drive1_pub = rospublisher('/semitruck/drive1_velocity_controller/command');    % Needs changed when changed to velcity controller
drive2_pub = rospublisher('/semitruck/drive2_velocity_controller/command');    % Needs changed when changed to velcity controller

left_steer_msg = rosmessage(left_steer_pub);
right_steer_msg = rosmessage(right_steer_pub);
drive1_msg = rosmessage(drive1_pub);
drive2_msg = rosmessage(drive2_pub);

% Publishers/Messages Combined for passing to functions

semi_pubs = [left_steer_pub; right_steer_pub; drive1_pub; drive2_pub];
semi_msgs = [left_steer_msg; right_steer_msg; drive1_msg; drive2_msg];

%% Subscribers

joint_states_sub = rossubscriber('/semitruck/joint_states',@semitruck_state_callback);

joint_data = receive(joint_states_sub);

semitruck_position_sub = rossubscriber('/gazebo/model_states',@semitruck_position_callback);

model_data = receive(semitruck_position_sub);
