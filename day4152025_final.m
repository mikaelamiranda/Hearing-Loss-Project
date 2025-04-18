% --- Load STL Pipe ---
stlFile = 'Pipe for matlab code.STL';
fv = stlread(stlFile);
vertices = fv.Points;
faces = fv.ConnectivityList;

% --- Toolbox-free PCA (using SVD) ---
centeredVertices = vertices - mean(vertices, 1);
[~, ~, V] = svd(centeredVertices, 'econ');
coeff = V;
pipeDir = coeff(:,1);           % Main axis of pipe
targetDir = [1; 0; 0];          % Desired direction: horizontal (X-axis)

% --- Find rotation matrix to align with X-axis ---
v = cross(pipeDir, targetDir);
s = norm(v);
if s ~= 0
    c = dot(pipeDir, targetDir);
    vx = [   0   -v(3)  v(2);
           v(3)   0   -v(1);
          -v(2)  v(1)   0 ];
    R = eye(3) + vx + vx^2 * ((1 - c)/(s^2));
else
    R = eye(3);
end

% --- Rotate STL ---
vertices_rot = (R * vertices')';

% --- New midpoint between selected points ---
P1 = [-2.05743, 13.2827, -1.08084];  % Updated P1
P2 = [2.06781, -13.4239, -1.0688];    % Updated P2
midpoint = (P1 + P2) / 2;             % Midpoint: [0.00519, -0.0706, -1.07482]
midpoint_rot = (R * midpoint')';

% --- Center the rotated STL ---
vertices_centered = vertices_rot - midpoint_rot;

% --- Create proper straight centerline along rotated pipe axis ---
minBound = min(vertices_centered);
maxBound = max(vertices_centered);
numPoints = 200;
x_centerline = linspace(minBound(1), maxBound(1), numPoints)';
y_centerline = zeros(size(x_centerline));
z_centerline = zeros(size(x_centerline));
centerline = [x_centerline, y_centerline, z_centerline];

% --- Plotting Section ---
figure;
trisurf(faces, vertices_centered(:,1), vertices_centered(:,2), vertices_centered(:,3), ...
    'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;

% 1. Main Pipe Centerline (Red)
plot3(centerline(:,1), centerline(:,2), centerline(:,3), 'r-', 'LineWidth', 3, 'DisplayName', 'Pipe Centerline');

% 2. Midpoint Line Along X-Axis (Green Dashed)
x_midline = linspace(-10, 10, 100)';  % Adjust range as needed
y_midline = zeros(size(x_midline));
z_midline = zeros(size(x_midline));
plot3(x_midline, y_midline, z_midline, 'g--', 'LineWidth', 2, 'DisplayName', 'Midpoint X-axis Line');

% 3. Original Midpoint Line (Magenta Dashed: Connects P1-P2)
P1_rot_centered = (R * P1')' - midpoint_rot;
P2_rot_centered = (R * P2')' - midpoint_rot;
plot3([P1_rot_centered(1), P2_rot_centered(1)], ...
      [P1_rot_centered(2), P2_rot_centered(2)], ...
      [P1_rot_centered(3), P2_rot_centered(3)], ...
      'm--', 'LineWidth', 2, 'DisplayName', 'Midpoint Line (P1-P2)');

% 4. Midpoint Marker (Yellow Dot)
plot3(0, 0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow', 'DisplayName', 'Midpoint');

axis equal; grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
legend show;
view(3);  % 3D view