clear; %clear all variables from memory
close all; %close all windows
clc; %clear command window

%------------------------------------------------------------------
%SETUP TRANSDUCER
%------------------------------------------------------------------
%Ultrasonic wave parameters
velocity = 1500;
frequency = 12.5e6;
wavelength = velocity/frequency;
k = 2*pi/wavelength; 
w=2*pi*frequency;

%Transducer array parameters
transducer_elements = 64;
N = transducer_elements;
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
grid_size = 20e-3;
grid_pts = 500;
x = transpose([-grid_size/2:grid_size/grid_pts:grid_size/2]);
z = [0:grid_size/grid_pts:grid_size];

%------------------------------------------------------------------
%SETUP MESHGRID
%------------------------------------------------------------------ 
[mx,mz,ms] = meshgrid(x,z,source_x_positions);

%------------------------------------------------------------------
%DEFINE FOCAL POINTS
%------------------------------------------------------------------
%for focal_points = [-2.25 -1.00 0.25 1.50 2.75 0.00 0.00 0.00 0.00 0.00; 5.00 7.50 10.00 12.50 15.00 5.00 7.50 10.00 12.50 15.00]/1000
for focal_points = [0;5]/1000
focal_points_index = round(transpose(grid_pts * (focal_points)/grid_size));
focal_points_index(1) = focal_points_index(1) + grid_pts/2;

%------------------------------------------------------------------
%CHECK SETUP
%------------------------------------------------------------------
near_field = transducer_width^2/(4*wavelength);
a_over_lambda = a/wavelength;
%------------------------------------------------------------------
%CREATE MATRICES AND VARIABLES
%-----------------------------------------------------------------
p = zeros([length(x)]); %prepare output matrix
A = 1; %complex number representing the size and phase of the source
t=0;

%------------------------------------------------------------------
%CALCULATE THE FIELD
%------------------------------------------------------------------
%Calculate r (distance from point source)
r = sqrt((mx-ms).^2 + mz.^2);

%Calculate p for all x and z
p = A .* r.^(-0.5) .* exp(1i.*(k.*r-w.*t));         

%Calculate Phi
angle = acos(mz./r);

%Calculate the Directivity Function
Df = sin(0.5*k*a*sin(angle))./(0.5*k*a*sin(angle));

%Calculate Time Delay
for ii = 1: length(source_x_positions)
    d0(ii) = sqrt(focal_points(2)^2 + focal_points(1)^2);
    dj(ii) = sqrt(focal_points(2)^2 + (focal_points(1) - (source_x_positions(ii)) )^2);
    tj(ii) = (dj(ii)-d0(ii))/velocity;
    B(ii) = exp(-1i*w*tj(ii));
end

%Calculate Field P
for ii = 1 : length(source_x_positions)    
    P_xz(:,:,ii) = B(ii) .* Df(1,ii) .* p(:,:,ii); 
end
p_xz = sum(P_xz,3);
p_xz = abs(p_xz);

%Calculate Phi
angle_q = 0.4229;

%------------------------------------------------------------------
%CALCULATE MAIN LOBE SHARPNESS FACTOR, q
%------------------------------------------------------------------
q = (1/pi)*((asin(sin(angle_q) + (wavelength/(transducer_elements*spacing)))) ...
    - asin(sin(angle_q) - (wavelength/(transducer_elements*spacing))))


%------------------------------------------------------------------
%CALCULATE GRATING LOBE INTENSITY
%------------------------------------------------------------------
max_grating_lobe_I = max(max(p_xz(grid_pts/2:end,1:100)));
max_beam_I = max(max(p_xz(1:end,100:end)));

grating_lobe_ratio = max_grating_lobe_I / max_beam_I;

%------------------------------------------------------------------
%OPTIMISATION TABLE
%------------------------------------------------------------------
 table(a_over_lambda,near_field,q,grating_lobe_ratio)
 
%------------------------------------------------------------------
%PLOT FIELD
%------------------------------------------------------------------
figure(1)
imagesc(p_xz)
hold on
title('Focused Beam Profile')
xticklabels = -10:2:10;
xticks = linspace(1, size(p_xz, 2), numel(xticklabels));
yticklabels = 0:2:20;
yticks = linspace(1, size(p_xz, 2), numel(yticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
xlabel('Z Position (mm)'); 
ylabel('X Position (mm)');
scatter(focal_points_index(1),focal_points_index(2),30,'r')
pause(1)

end