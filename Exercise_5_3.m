clear; %clear all variables from memory
close all; %close all windows
clc; %clear command window
tic

%------------------------------------------------------------------
%SETUP TRANSDUCER
%------------------------------------------------------------------
%Ultrasonic wave parameters
velocity = 1500;
frequency = 12.5e6;
wavelength = velocity/frequency;
k = 2*pi/wavelength; 
w=2*pi*frequency;
distance_to_wall = 20/1000;

%Transducer array parameters
transducer_elements = 64;
element_width = 0.05e-3;
spacing = 0.05e-3;
transducer_pitch = element_width + spacing;
transducer_width = transducer_pitch * (transducer_elements-1);
a = element_width;
sources_per_wavelength = transducer_elements / (transducer_width/wavelength);

%Setup Transducer Sources
source_x_positions = linspace((-transducer_width/2),(transducer_width/2),transducer_elements);

%------------------------------------------------------------------
%SETUP GRID
%------------------------------------------------------------------
grid_size = 22e-3;
grid_pts = 200;
x = transpose([-grid_size/2:grid_size/grid_pts:grid_size/2]);
z = [0:grid_size/grid_pts:grid_size];

%------------------------------------------------------------------
%SETUP MESHGRID
%------------------------------------------------------------------ 
[mx,mz,ms] = meshgrid(x,z,source_x_positions);

%------------------------------------------------------------------
%DEFINE FOCAL POINTS
%------------------------------------------------------------------
point_reflectors = [-2.25 -1.00 0.25 1.50 2.75 0.00 0.00 0.00 0.00 0.00; 5.00 7.50 10.00 12.50 15.00 5.00 7.50 10.00 12.50 15.00];
point_reflectors = point_reflectors/1000;

%------------------------------------------------------------------
%RETRIEVE INDEXES OF FOCAL POINTS
%------------------------------------------------------------------
focal_points_index = round(transpose(grid_pts * (point_reflectors)/grid_size));
focal_points_index(:,1) = focal_points_index(:,1) + grid_pts/2;
back_wall_index = round(transpose(grid_pts * (distance_to_wall)/grid_size));

%------------------------------------------------------------------
%CREATE MATRICES AND VARIABLES
%------------------------------------------------------------------
A = 1; %complex number representing the size and phase of the source
A_scatter = 0.01;
focal_points = transpose(point_reflectors); %Create Matrix of All Focal Points

%------------------------------------------------------------------
%SETUP HANNING WINDOW
%------------------------------------------------------------------
r = sqrt((mx-ms).^2 + mz.^2);
N = 5;
fft_points = 2048;
time_centre_pulse = (N/2)*wavelength/velocity;
pulse_time = (N*wavelength)/velocity;
max_time = 2*distance_to_wall/velocity + pulse_time;
%max_time = 2*r(end,grid_pts/2,1)/ velocity + pulse_time;
time = linspace(0,max_time,fft_points);
peak_pos_fract = (pulse_time./2)./max_time;
window = fn_hanning(fft_points, peak_pos_fract, peak_pos_fract);
input_signal = window.' .* sin(2.*pi.*frequency.*time - (time_centre_pulse));
 
%------------------------------------------------------------------
%CALCULATE FFT OF HANNING WINDOW
%------------------------------------------------------------------ 
fftpoints = 2.^nextpow2(length(time));
fourier_transform = fft(input_signal,fftpoints);
%fourier_transform(fftpoints/2:end) = [];
fourier_transform = fourier_transform(1:end/2);
frequency_sample = 1./(max(gradient(time)));
frequency_step = frequency_sample./fftpoints;
frequency_fourier = [0:frequency_step:frequency_step.*(length(fourier_transform)-1)];

%------------------------------------------------------------------
%CREATE 3D MATRICES FOR VARIABLES
%------------------------------------------------------------------
fourier_transform = repmat(fourier_transform, [transducer_elements 1 transducer_elements]);
fourier_transform = permute(fourier_transform, [3 1 2]);

k = 2.*pi.*frequency_fourier./velocity;
k_a = repmat(k, [transducer_elements 1]);
k_b = repmat(k, [transducer_elements 1 transducer_elements]);
k_b = permute(k_b, [3 1 2]);
k_a = permute(k_a, [3 1 2]);

H_txrx = zeros(transducer_elements,transducer_elements,length(frequency_fourier));
for ii = 1 : length(focal_points)
    %------------------------------------------------------------------
    %CALCULATE DISTANCES
    %------------------------------------------------------------------
    R_r = sqrt( (focal_points(ii,1) - source_x_positions(:)).^2 + (focal_points(ii,2)).^2);
    R_repmat = repmat(R_r, [1 transducer_elements]);
    R_permute = permute(R_repmat, [2 1]);
    d_txrx = R_permute + R_repmat;
    d_txrx = repmat(d_txrx, [1 1 size(frequency_fourier,2)]);

    %------------------------------------------------------------------
    %COMPLEX SPECTRUM OF PHASE SHIFTED SIGNAL
    %------------------------------------------------------------------
    G_txrx = fourier_transform .* exp(-1i .* k_b .* d_txrx);
    
    %------------------------------------------------------------------
    %CALCULATE AMPLITUDE
    %------------------------------------------------------------------
    A_txrx = A_scatter ./ sqrt(R_repmat .* R_permute);
    A_txrx = repmat(A_txrx, [1 1 length(fourier_transform)]);

    %------------------------------------------------------------------
    %CALCULATE ANGLE
    %------------------------------------------------------------------
    angle_points = atan( (focal_points(ii,1) - source_x_positions(:)) ./ focal_points(ii,2));
    angle_points = repmat(angle_points, [1 transducer_elements length(fourier_transform)]);
    angle_points(isnan(angle_points))=1;
    %angle_points = permute(angle_points, [2 1 3]);
    
    %------------------------------------------------------------------
    %CALCULATE DIRECTIVITY
    %------------------------------------------------------------------
    %p = pi * sinc(pi .* a .* k_a .* sin(angle_points));
    ptx = sin( 0.5 .* k_b .* a .* sin(angle_points))./( 0.5 .* k_b .* a .* sin(angle_points));
    %p = pi * sinc( a .* sin(angle_points).* k_a);
    %p1 = repmat(p, [transducer_elements 1 1]);
    prx = permute(ptx, [2 1 3]);
    
%     H_txr = p2 .* p1;
%     plot(H_txr(:,:,100))
%     
%     pause(0.5)
    
    %------------------------------------------------------------------
    %CALCULATE RESULTING SPECTRUM
    %------------------------------------------------------------------
    H_txrx_old =  A_txrx .* G_txrx;% .* ptx .* prx;
    H_txrx = H_txrx + H_txrx_old;
end

%------------------------------------------------------------------
%SAME AGAIN FOR BACK WALL
%------------------------------------------------------------------
r_wall = repmat(source_x_positions, [transducer_elements 1]);
r_wall = (r_wall - permute(r_wall, [2 1])) ./ 2;
r_wall = sqrt(r_wall.^2 + distance_to_wall.^2);
r_wall = 2 .* r_wall;
d_txrx_wall = repmat(r_wall, [1 1 length(fourier_transform)]);

A_txrx_wall = A ./ sqrt(d_txrx_wall);

G_txrx_wall = fourier_transform .* exp(-1i .* k_b .* d_txrx_wall);

%------------------------------------------------------------------
%DIRECTIVITIES
%------------------------------------------------------------------
angle_points_wall = acos( (distance_to_wall) ./ (d_txrx_wall ./ 2));
p_wall = sin( 0.5 .*  k_b .* a .* sin(angle_points_wall))./( 0.5 .* k_b .* a .* sin(angle_points_wall));
p_wall(isnan(p_wall))=1;
p_t_wall = permute(p_wall, [2 1 3]);

H_txrx_wall =  p_wall .* p_t_wall .* A_txrx_wall .* G_txrx_wall;

%------------------------------------------------------------------
%SUM SIGNALS FOR POINTS AND WALL
%------------------------------------------------------------------
H_txrx =   H_txrx + H_txrx_wall;

%------------------------------------------------------------------
%CALCULATE TIME DOMAIN SIGNAL
%------------------------------------------------------------------
time_domain = ifft(H_txrx,size(time,2),3);

trans = 1;
rec = 1;

plot_dat = time_domain(:,:,1:length(time));

plot_dat = plot_dat(trans,rec,:);
plot_dat = (permute(plot_dat, [3 1 2]));
hold off
plot(time(1:length(plot_dat)),(plot_dat));
hold on

title('Time Domain Signal for Transmitter 1 & Receiver 1')
xlabel('Time (s)'); 
ylabel('Amplitude');
caxis([0 150])
