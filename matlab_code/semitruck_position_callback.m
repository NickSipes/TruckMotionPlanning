function semitruck_position_callback(~,message)
    
    global semitruck_position

    % Store Current semitruck position
    semitruck_position(1,1) = message.Pose(2).Position.X;
    semitruck_position(2,1) = message.Pose(2).Position.Y;
    semitruck_position(3,1) = message.Pose(2).Position.Z;
    
end