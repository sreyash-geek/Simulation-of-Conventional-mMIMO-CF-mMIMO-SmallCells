clc;
clear
close all;

%% Simulation Parameters
wavelength = 10;                % wavelength in cm
squareSide = 10 * wavelength;   % side length of the simulation area

% Generate grid of sample points where normalized SNR is computed
x = 4:1:(squareSide-4);
y = 4:1:(squareSide-4);

% Define antenna locations on the perimeter of a square (as in original script)
pointsperdim = 10;  % Number of antennas per side in the square
antenna_locations = squareSide * [ ...
    linspace(0, 1, pointsperdim+1), ...                              % bottom edge
    1i * linspace(1/pointsperdim, 1, pointsperdim), ...                % right edge
    1i + linspace(1/pointsperdim, 1, pointsperdim), ...                % top edge
    1 + 1i * linspace(1/pointsperdim, 1 - 1/pointsperdim, pointsperdim-1)];  % left edge

% Set the original target location (used in static plots)
targetlocation = round(squareSide/4) + 1i*round(squareSide/2);

%% Preallocation for results 
numX = length(x);
numY = length(y);
numAnt = length(antenna_locations);
phaseshifts = zeros(numX, numY, numAnt);
distances = zeros(numX, numY, numAnt);

% Compute distances and phase-shifts for each antenna and grid point (target-independent)
for n = 1:numAnt
    for k = 1:numX
        for j = 1:numY
            samplePoint = x(k) + 1i * y(j);
            distances(k,j,n) = abs(samplePoint - antenna_locations(n));
            phaseshifts(k,j,n) = 2*pi * distances(k,j,n) / wavelength;
        end
    end
end

% For the original target location, compute phase term (target-dependent)
angletarget = zeros(numX, numY, numAnt);
for n = 1:numAnt
    angletarget(:,:,n) = 2*pi * abs(targetlocation - antenna_locations(n)) / wavelength;
end

% Compute distances from all antennas to the target (for normalization)
distances_target = abs(antenna_locations - targetlocation);

%% Combined Channel Gain - Normalized SNR
% Free-space propagation pathloss is assumed
channel_gain = abs(sum(exp(1i*(phaseshifts - angletarget)) ./ sqrt(4*pi*distances.^2), 3)).^2;
channel_gain_target = abs(sum(1 ./ sqrt(4*pi*distances_target.^2))).^2; % at target, phase diff = 0
normalization = channel_gain_target; % normalization factor

%% Insightful Plots
% Part a: 3D Surface Plot with antenna markers and target marker
figure;
surf(y/wavelength, x/wavelength, channel_gain/normalization);
colormap(flipud(parula));
colorbar;
shading interp;
hold on;
for n = 1:numAnt
    plot3(imag(antenna_locations(n))/wavelength, real(antenna_locations(n))/wavelength, 0, 'k*', 'MarkerSize', 5);
end
plot3(imag(targetlocation)/wavelength, real(targetlocation)/wavelength, channel_gain_target/normalization, ...
    'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
hold off;
xlabel('Distance (wavelengths)');
ylabel('Distance (wavelengths)');
zlabel('Normalized SNR');
set(gca, 'fontsize', 16);
xlim([0 squareSide/wavelength]);
ylim([0 squareSide/wavelength]);
view([-20 26]);

% Part b: 2D Overhead View
figure;
surf(y/wavelength, x/wavelength, channel_gain/normalization);
colormap(flipud(parula));
colorbar;
shading interp;
hold on;
for n = 1:numAnt
    plot3(imag(antenna_locations(n))/wavelength, real(antenna_locations(n))/wavelength, 0, 'k*', 'MarkerSize', 5);
end
plot3(imag(targetlocation)/wavelength, real(targetlocation)/wavelength, channel_gain_target/normalization, ...
    'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
hold off;
xlabel('Distance (wavelengths)');
ylabel('Distance (wavelengths)');
zlabel('Normalized SNR');
set(gca, 'fontsize', 16);
view(2);
xlim([0 squareSide/wavelength]);
ylim([0 squareSide/wavelength]);
axis equal;

%% 3D Phase Difference Surface Plot (for 1st antenna)
% Compute phase difference (difference between actual phase and target phase) for the first antenna
phase_diff = phaseshifts(:,:,1) - angletarget(:,:,1);
figure;
surf(y/wavelength, x/wavelength, phase_diff);
colormap(jet);
colorbar;
shading interp;
xlabel('Distance (wavelengths)');
ylabel('Distance (wavelengths)');
zlabel('Phase Difference (rad)');
title('Phase Difference Surface (Antenna 1)');
set(gca, 'fontsize', 16);
view([-30 30]);
rotate3d on;

%% Animated 3D Plot for a Moving Target
% Define a trajectory for the moving target along a horizontal line
numFrames = 20;
target_start = squareSide/4 + 1i*(squareSide/2);
target_end   = 3*squareSide/4 + 1i*(squareSide/2);
target_traj = linspace(target_start, target_end, numFrames);

figure;
hSurf = surf(y/wavelength, x/wavelength, zeros(numX, numY));
colormap(flipud(parula));
colorbar;
shading interp;
hold on;
% Plot antenna locations (static)
for n = 1:numAnt
    plot3(imag(antenna_locations(n))/wavelength, real(antenna_locations(n))/wavelength, 0, 'k*', 'MarkerSize', 5);
end
hTarget = plot3(0,0,0,'ko','MarkerSize',5,'MarkerFaceColor','k');
xlabel('Distance (wavelengths)');
ylabel('Distance (wavelengths)');
zlabel('Normalized SNR');
title('Animated Normalized SNR with Moving Target');
set(gca, 'fontsize', 16);
xlim([0 squareSide/wavelength]);
ylim([0 squareSide/wavelength]);
zlim([0 1.2]);
rotate3d on;

% Animation loop: update target location, recompute phase terms and channel gain
for frame = 1:numFrames
    new_target = target_traj(frame);
    % Recompute target-dependent phase shifts for new target
    new_angletarget = zeros(numX, numY, numAnt);
    for n = 1:numAnt
        new_angletarget(:,:,n) = 2*pi * abs(new_target - antenna_locations(n)) / wavelength;
    end
    % Recompute channel gain for new target
    channel_gain_new = abs(sum(exp(1i*(phaseshifts - new_angletarget)) ./ sqrt(4*pi*distances.^2), 3)).^2;
    % Normalization: channel gain at new target (phase diff zero)
    distances_target_new = abs(antenna_locations - new_target);
    channel_gain_target_new = abs(sum(1 ./ sqrt(4*pi*distances_target_new.^2)))^2;
    norm_SNR_new = channel_gain_new / channel_gain_target_new;
    
    % Update the surface plot
    set(hSurf, 'ZData', norm_SNR_new);
    % Update target marker
    set(hTarget, 'XData', imag(new_target)/wavelength, 'YData', real(new_target)/wavelength, ...
        'ZData', 1);  % place marker at normalized gain =1
    title(sprintf('Frame %d: Target = (%.1f, %.1f)', frame, real(new_target), imag(new_target)));
    
    drawnow;
    pause(1);  % pause for 0.2 sec between frames
end
