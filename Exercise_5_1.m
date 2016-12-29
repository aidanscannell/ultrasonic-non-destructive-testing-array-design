clear; %clear all variables from memory
close all; %close all windows
clc; %clear command window

%------------------------------------------------------------------
%SETUP TRANSDUCER
%------------------------------------------------------------------
%Ultrasonic wave parameters
velocity = 1500;
frequency = 5e6;
wavelength = velocity/frequency;
k = 2*pi/wavelength;
w=2*pi*frequency;

%Transducer array parameters
transducer_elements = 20;
element_width = 0.09e-3;
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
%CREATE MATRICES AND VARIABLES
%-----------------------------------------------------------------
p = zeros([length(x)]); %prepare output matrix
A = 1; %complex number representing the size and phase of the source
t=0;

tic
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

%Calculate Field P
P_xz = Df .* p;

%------------------------------------------------------------------
%PLOT FIELD
%------------------------------------------------------------------
p_xz = sum(P_xz,3);
p_xz = abs(p_xz);

%Create Figure
figure(1)
imagesc(p_xz)
hold on
title('Focused Beam Profile From a Linear Array')
xticklabels = -10:2:10;
xticks = linspace(1, size(p_xz, 2), numel(xticklabels));
yticklabels = 0:2:20;
yticks = linspace(1, size(p_xz, 2), numel(yticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
xlabel('X Position (mm)'); 
ylabel('Z Position (mm)');
caxis([0 150])