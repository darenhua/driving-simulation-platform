function next_state = vehicle_model(servo_motor, drive_motor,state)
    driveL = drive_motor;
    driveR = drive_motor;
    [v, steering_rate, steering_angle] = steering_system(servo_motor, driveL, driveR);
    next_state = kinematic_model(state, v, steering_rate, steering_angle);
end            
