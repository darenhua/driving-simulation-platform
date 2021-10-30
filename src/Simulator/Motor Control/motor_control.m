function [drive, servo] = motor_control(centerIdx, drive_gain, steering_gain)
% This function decides how the virtual motors should activate depending on
% the center of the camera waveform.
% ========================================================================
% INPUTS:
% - The predicted center of the camera waveform
% - The gain/multiplier on the drive motor strength
% - The gain/multiplier on the servo (steering) motor strength
% OUTPUTS:
% - Strength of drive motor [-1,1]
% - Strength of servo motor [-30, 30]

drive = drive_gain;
servo = (centerIdx - 64)*.4*steering_gain;
end
