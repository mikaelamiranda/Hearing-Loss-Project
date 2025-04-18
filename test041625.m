clc;
clear;

%% ==== Load STL Pipe ====
stlFile = 'Pipe for matlab code.STL';
fv = stlread(stlFile);
vertices = fv.Points;
faces = fv.ConnectivityList;

% Plot STL pipe
figure;
trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), ...
    'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Tool Moving Through Curved Pipe (XZ Bend)');
lighting gouraud;
camlight;
hold on;

%% ==== Define Centerline ====
tool_radius = 2;
y_base = 11.623;  % fixed height of the pipe
R = 40;           % Radius of curvature
theta = linspace(0, deg2rad(140), 70)';  % arc from 0 to 140 degrees

% --- Straight Segment (X increasing, constant Y, Z = 0)
n_straight = 50;
x_straight = linspace(0, 90, n_straight)';
y_straight = ones(n_straight,1) * y_base;
z_straight = zeros(n_straight,1);

% --- Curved Segment (Bend in XZ plane, starting at x = 90)
x_center = 90;            % center of curve on X
z_center = 0 + R;         % center offset up in Z

x_arc = x_center + R * sin(theta);
z_arc = z_center - R * cos(theta);
y_arc = ones(size(x_arc)) * y_base;  % Y stays constant

% Combine centerline segments
x_path = [x_straight; x_arc];
y_path = [y_straight; y_arc];
z_path = [z_straight; z_arc];

centerline = [x_path, y_path, z_path];

% Plot centerline
plot3(x_path, y_path, z_path, 'k--', 'LineWidth', 2);

%% ==== Animate Tool Along Centerline ====
tool = sphereSurface(centerline(1,:), tool_radius);
h_tool = surf(tool.X, tool.Y, tool.Z, 'FaceColor', 'red', 'EdgeColor', 'none');

for i = 1:size(centerline,1)
    [tool.X, tool.Y, tool.Z] = sphereSurface(centerline(i,:), tool_radius);
    set(h_tool, 'XData', tool.X, 'YData', tool.Y, 'ZData', tool_
% Your main script content (loading STL, plotting, animating, etc.)

%% ==== Helper Function ====
function s = sphereSurface(center, radius)
    [X, Y, Z] = sphere(20);
    s.X = X * radius + center(1);
    s.Y = Y * radius + center(2);
    s.Z = Z * radius + center(3);
end
