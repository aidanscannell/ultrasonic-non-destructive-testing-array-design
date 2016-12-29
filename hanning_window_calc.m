function [window, input_signal,time,distance] = hanning_window_calc(velocity,frequency,wavelength,distance_to_wall)

%Number of cycles
N = 5; 

%Time
total_time = (2*distance_to_wall)/velocity;
duration_of_pulse = N/frequency;
max_time = total_time + 5*duration_of_pulse;
time_step = 1/(frequency*400); %time step for 10 points per cycle
time = [0:time_step:max_time];

%Distance
distance_step = wavelength/4; %distance step for 4 points per cycle
max_distance = 2*distance_to_wall + (N*wavelength);
distance = [0:distance_step:max_distance];

%Wave Propagation
nondispersive_propagation = 1; 
velocity_at_centre_frequency = velocity; 
centre_frequency = frequency;

%Hanning Windowed Function Inputs
time_at_centre_of_pulse = -(N/2)*wavelength/velocity; 
number_of_points = size(time,2);
peak_pos_fract = (N*wavelength/2)/(max_distance); %(N*5)/(number_of_points)
half_width_fract = (N/frequency)/(2*max_time);

%Create Input Signal Using Hanning Window Function
window = fn_hanning(number_of_points, peak_pos_fract, half_width_fract);
sine_wave = sin(2*pi*centre_frequency*time+(time_at_centre_of_pulse));
input_signal = window(:) .* sine_wave(:);
end