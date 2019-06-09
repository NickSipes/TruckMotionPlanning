% All of these closed-form solutions hold true with the following
% conditions:
% Constant velocity of the the tractor
% Constant steering angle for the tractor

% Inital Conditions
x_inital = 0;
y_inital = 0;
inital_tractor_orientation = 0;
inital_trailer_to_tractor_angle = 0;

%% Kinematic parameters of a tractor-trailer robot

% Independent Variables
steering_angle = 15; % degrees
tractor_velocity = 1; % m/s

x = x_inital;
y = y_inital;
theta = inital_tractor_orientation;
alpha = inital_trailer_to_tractor_angle;

% Main execution loop
% Edit steering angle and tractor verlociy at different times to see the
% motion
for time = 1:500
    
    if time < 5
        steering_angle = 0;
        tractor_velocity = 10;
    else
       steering_angle = -1 
       tractor_velocity = -2;
    end
    
    output = computeKinematicResults(steering_angle, tractor_velocity, x, y,...
        theta, alpha, time);

    data(time,:) = real(output);
    alpha = real(output(1,1));
    theta = real(output(1,2));
    x = real(output(1,3));
    y = real(output(1,4));
end

% Plot the tractor motion
comet(data(:,3),data(:,4));


% Kinematic Equations
function results = computeKinematicResults(steering_angle, tractor_velocity, previous_x, previous_y,...
    previous_tractor_orientation, previous_trailer_to_tractor_angle, time)

    % Physical Constants of tractor-trailer robot
    tractor_length = 1;
    trailer_length = 3;
    
    % Simplified case if the steering angle is 0
    if steering_angle == 0
        
        trailer_to_tractor_angle = previous_trailer_to_tractor_angle;
        tractor_x_position = previous_x + sind(previous_tractor_orientation);
        tractor_y_position = previous_y - cosd(previous_tractor_orientation);
        tractor_orientation = previous_tractor_orientation;
        
        results(1,1) = trailer_to_tractor_angle;
        results(1,2) = tractor_orientation;
        results(1,3) = tractor_x_position;
        results(1,4) = tractor_y_position;
        
    else
        % Constants of intergration
        C_x = previous_x;
        C_y = previous_y;
        C_tractor_orientation = previous_tractor_orientation;
        epsilon = (tand(steering_angle))^2 * (trailer_length^2 - tractor_length^2);
        C_trailer_to_tractor_angle = 2 * ((atan(epsilon^-0.5*(-tand(0.5 * previous_trailer_to_tractor_angle)...
                * tand(steering_angle) * trailer_length - tractor_length)) * tractor_length ^ 2 * trailer_length)...
                / (sqrt(epsilon) * tractor_velocity^2 * tand(steering_angle)));

        % Alpha in the diagram
        trailer_to_tractor_angle = 2 * atand((tand(steering_angle) * trailer_length)^-1 * ...
            (tand((-tractor_velocity * sqrt(epsilon) * (time * tractor_length +...
            (C_trailer_to_tractor_angle * tractor_velocity * tand(steering_angle))))...
            / (2 * tractor_length^2 * trailer_length^2)) * sqrt(epsilon) - tractor_length));
        results(1,1) = trailer_to_tractor_angle;

        % Theta in the diagram
        tractor_orientation = C_tractor_orientation + (tractor_velocity / tractor_length)...
                                * tand(steering_angle);
        results(1,2) = tractor_orientation;

        % X position of the tractor
        tractor_x_position = C_x + (tractor_length/tand(steering_angle)) * ...
                            (sind(tractor_orientation + (tractor_velocity * time /tractor_length) * tand(steering_angle)) ...
                             - sind(tractor_orientation));
        results(1,3) = tractor_x_position;

        % Y position of the tractor
        tractor_y_position = C_y - (tractor_length/tand(steering_angle)) * ...
                            (cosd(tractor_orientation + (tractor_velocity * time /tractor_length) * tand(steering_angle)) ...
                            - cosd(tractor_orientation));
        results(1,4) = tractor_y_position;
    end
end
