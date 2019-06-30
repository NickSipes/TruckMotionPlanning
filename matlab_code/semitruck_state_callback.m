function semitruck_state_callback(~,message)
    
    global joint_states joint_states_time
    
    % Store Time Associated with the States
    joint_states_time(1,1) = message.Header.Stamp.Sec;
    joint_states_time(2,1) = message.Header.Stamp.Nsec;
    
    % Store Current Joint Position and Velocity
    joint_states(:,1) = message.Position;
    joint_states(:,2) = message.Velocity;
    
end