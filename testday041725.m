function simulatePipeMotion(stlFilename)
fv = stlread('Pipe for matlab code.STL');

    % Load the STL file using MATLAB's built-in function
    vertices = fv.Points;
    faces = fv.ConnectivityList;

    % Extract an improved centerline
    centerline = improvedCenterlineExtraction(vertices, faces);
    
    % Define tool properties
    baseRadius = 0.475/2; % Convert diameter to radius in mm
    tipRadius = 0.35/2;   % Convert diameter to radius in mm
    toolLength = 5.5;     % Total length in mm
    
    % Number of simulation steps
    numSteps = 100;
    
    % Perform the simulation and visualization
    simulateToolMotion(vertices, faces, centerline, baseRadius, tipRadius, toolLength, numSteps);
end

function centerline = improvedCenterlineExtraction(vertices, faces)
    % Get bounding box
    minCoords = min(vertices);
    maxCoords = max(vertices);
    
    % Find the longest dimension (likely the main pipe axis)
    [~, mainAxis] = max(maxCoords - minCoords);
    secondAxis = mod(mainAxis, 3) + 1;
    thirdAxis = mod(secondAxis, 3) + 1;
    
    % Define number of slices for cross-sectional analysis
    numSlices = 100;
    
    % Define slice positions along the main axis
    slicePositions = linspace(minCoords(mainAxis), maxCoords(mainAxis), numSlices);
    
    % Initialize storage for cross-section centers
    crossSectionCenters = zeros(numSlices, 3);
    validSlices = false(numSlices, 1);
    
    % For each slice position, find vertices near the slice plane
    for i = 1:numSlices
        slicePos = slicePositions(i);
        
        % Find vertices within a small tolerance of the slice position
        tolerance = (maxCoords(mainAxis) - minCoords(mainAxis)) / (numSlices * 2);
        sliceIndices = find(abs(vertices(:, mainAxis) - slicePos) < tolerance);
        
        if length(sliceIndices) > 10  % Ensure we have enough points for a meaningful cross-section
            % Extract slice vertices
            sliceVertices = vertices(sliceIndices, :);
            
            % Find the boundary of the cross-section (convex hull in 2D)
            % Project vertices onto the plane perpendicular to the main axis
            projectedVertices = sliceVertices(:, [secondAxis, thirdAxis]);
            
            try
                % Try to compute convex hull
                k = convhull(projectedVertices);
                boundaryVertices = projectedVertices(k, :);
                
                % Estimate center as the centroid of the boundary polygon
                crossSectionCenters(i, secondAxis) = mean(boundaryVertices(:, 1));
                crossSectionCenters(i, thirdAxis) = mean(boundaryVertices(:, 2));
                crossSectionCenters(i, mainAxis) = slicePos;
                validSlices(i) = true;
            catch
                % Skip if convex hull fails (can happen with very few points)
                continue;
            end
        end
    end
    
    % Keep only valid slice centers
    centerline = crossSectionCenters(validSlices, :);
    
    % Ensure we have enough points for the centerline
    if size(centerline, 1) < 3
        error('Failed to extract enough valid cross-sections for centerline.');
    end
    
    % Sort points along the main axis to ensure proper ordering
    [~, sortIdx] = sort(centerline(:, mainAxis));
    centerline = centerline(sortIdx, :);
    
    % Apply smoothing to the centerline using moving average
    centerline = smoothCenterlineMovingAverage(centerline);
    
    % Increase the point density along the centerline for smoother motion
    centerline = densifyCenterline(centerline, 200);
end

function smoothedCenterline = smoothCenterlineMovingAverage(centerline)
    % Apply a moving average filter to smooth the centerline
    windowSize = 5;
    numPoints = size(centerline, 1);
    smoothedCenterline = centerline;
    
    % Skip smoothing if we don't have enough points
    if numPoints <= windowSize
        return;
    end
    
    % Apply moving average for each dimension
    for i = 1:numPoints
        startIdx = max(1, i - floor(windowSize/2));
        endIdx = min(numPoints, i + floor(windowSize/2));
        smoothedCenterline(i, :) = mean(centerline(startIdx:endIdx, :), 1);
    end
    
    % Apply the smoothing multiple times for a smoother result
    for j = 1:3
        for i = 1:numPoints
            startIdx = max(1, i - floor(windowSize/2));
            endIdx = min(numPoints, i + floor(windowSize/2));
            smoothedCenterline(i, :) = mean(smoothedCenterline(startIdx:endIdx, :), 1);
        end
    end
end

function denseCenterline = densifyCenterline(centerline, numDesiredPoints)
    % Create a denser set of points along the centerline using interpolation
    
    % Calculate cumulative distances
    diffs = diff(centerline);
    segmentLengths = sqrt(sum(diffs.^2, 2));
    cumDist = [0; cumsum(segmentLengths)];
    
    % Create new parameter values
    newDist = linspace(0, cumDist(end), numDesiredPoints);
    
    % Interpolate each coordinate
    denseCenterline = zeros(numDesiredPoints, 3);
    for dim = 1:3
        denseCenterline(:,dim) = interp1(cumDist, centerline(:,dim), newDist, 'pchip');
    end
end

function [cylinderVertices, cylinderFaces] = createTaperedCylinder(startPoint, endPoint, baseRadius, tipRadius, nSides)
    % Create a tapered cylinder from startPoint to endPoint
    
    % Calculate direction vector
    direction = endPoint - startPoint;
    length = norm(direction);
    if length < eps
        error('Start and end points are too close');
    end
    
    direction = direction / length;
    
    % Create orthogonal vectors for the cylinder cross-section
    if abs(direction(3)) < 0.9
        perpVector = cross(direction, [0, 0, 1]);
    else
        perpVector = cross(direction, [1, 0, 0]);
    end
    perpVector = perpVector / norm(perpVector);
    
    perpVector2 = cross(direction, perpVector);
    perpVector2 = perpVector2 / norm(perpVector2);
    
    % Generate vertices around the base (larger end)
    basePoints = zeros(nSides, 3);
    for i = 1:nSides
        angle = 2 * pi * (i-1) / nSides;
        basePoints(i, :) = startPoint + baseRadius * (cos(angle) * perpVector + sin(angle) * perpVector2);
    end
    
    % Generate vertices around the tip (smaller end)
    tipPoints = zeros(nSides, 3);
    for i = 1:nSides
        angle = 2 * pi * (i-1) / nSides;
        tipPoints(i, :) = endPoint + tipRadius * (cos(angle) * perpVector + sin(angle) * perpVector2);
    end
    
    % Combine vertices (base circle, tip circle, centers)
    cylinderVertices = [basePoints; tipPoints; startPoint; endPoint];
    
    % Create triangular faces for the cylinder sides
    cylinderFaces = zeros(2*nSides, 3);
    for i = 1:nSides
        next = mod(i, nSides) + 1;
        
        % First triangle of the quad
        cylinderFaces(2*i-1, :) = [i, next, i + nSides];
        
        % Second triangle of the quad
        cylinderFaces(2*i, :) = [next, next + nSides, i + nSides];
    end
    
    % Create faces for the base cap
    baseFaces = zeros(nSides, 3);
    for i = 1:nSides
        next = mod(i, nSides) + 1;
        baseFaces(i, :) = [2*nSides + 1, next, i];
    end
    
    % Create faces for the tip cap
    tipFaces = zeros(nSides, 3);
    for i = 1:nSides
        next = mod(i, nSides) + 1;
        tipFaces(i, :) = [2*nSides + 2, i + nSides, next + nSides];
    end
    
    % Combine all faces
    cylinderFaces = [cylinderFaces; baseFaces; tipFaces];
end

function simulateToolMotion(pipeVertices, pipeFaces, centerline, baseRadius, tipRadius, toolLength, numSteps)
    % Create figure for visualization
    figure('Position', [100, 100, 1200, 800]);
    
    % Display the pipe with transparency
    h_pipe = trisurf(pipeFaces, pipeVertices(:,1), pipeVertices(:,2), pipeVertices(:,3), ...
        'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold on;
    
    % Display the centerline
    plot3(centerline(:,1), centerline(:,2), centerline(:,3), 'r-', 'LineWidth', 1.5);
    
    % Set up the axes
    axis equal;
    grid on;
    xlabel('X (mm)');
    ylabel('Y (mm)');
    zlabel('Z (mm)');
    title('Simulation of Tapered Tool Moving Through Pipe');
    view(3);
    light('Position',[1 1 1],'Style','infinite');
    lighting gouraud;
    
    % Calculate parameterized positions along the centerline
    centerlinePoints = size(centerline, 1);
    
    % Create waypoints for the tool motion
    numWaypoints = max(numSteps, centerlinePoints);
    waypoints = interpolateWaypoints(centerline, numWaypoints);
    
    % Calculate tangent directions along the centerline
    tangents = calculateTangents(waypoints);
    
    % Initialize the tool
    nSides = 16; % Number of sides for the tool cylinder
    
    % Main simulation loop
    
    % Create initial tool
    halfLength = toolLength / 2;
    
    % Place the tool at the start with correct orientation
    startPt = waypoints(1,:) - halfLength * tangents(1,:);
    endPt = waypoints(1,:) + halfLength * tangents(1,:);
    [toolVertices, toolFaces] = createTaperedCylinder(startPt, endPt, baseRadius, tipRadius, nSides);
    
    % Initial tool visualization
    h_tool = trisurf(toolFaces, toolVertices(:,1), toolVertices(:,2), toolVertices(:,3), ...
        'FaceColor', [0.2, 0.6, 0.9], 'EdgeColor', 'none');
    
    % Store the collision state and tool path
    collisionHistory = false(numSteps, 1);
    toolPath = zeros(numSteps, 3);
    
    % Main animation loop
    for i = 1:numSteps
        % Current position index along the waypoints
        posIdx = min(i, size(waypoints, 1));
        
        % Current position on the centerline
        currentPos = waypoints(posIdx,:);
        toolPath(i,:) = currentPos;
        
        % Current direction (tangent to the centerline)
        currentTangent = tangents(posIdx,:);
        
        % Place the tool with its center at the current position
        startPt = currentPos - halfLength * currentTangent;
        endPt = currentPos + halfLength * currentTangent;
        
        % Create tapered tool
        [toolVertices, toolFaces] = createTaperedCylinder(startPt, endPt, baseRadius, tipRadius, nSides);
        
        % Check for collisions
        collisionDetected = checkCollision(toolVertices, pipeVertices, pipeFaces);
        collisionHistory(i) = collisionDetected;
        
        % Update tool color based on collision state
        if collisionDetected
            toolColor = [0.9, 0.2, 0.2]; % Red for collision
        else
            toolColor = [0.2, 0.6, 0.9]; % Blue for no collision
        end
        
        % Update tool visualization
        set(h_tool, 'Vertices', toolVertices);
        set(h_tool, 'Faces', toolFaces);
        set(h_tool, 'FaceColor', toolColor);
        
        % Add current position marker
        plot3(currentPos(1), currentPos(2), currentPos(3), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 4);
        
        % Update display
        drawnow;
        pause(0.05);
    end
    
    % Display collision statistics
    collisionCount = sum(collisionHistory);
    disp(['Simulation complete: ' num2str(collisionCount) ' collisions detected out of ' num2str(numSteps) ' steps.']);
    
    % Plot path with color based on collision
    figure;
    for i = 1:numSteps-1
        if collisionHistory(i)
            line(toolPath(i:i+1,1), toolPath(i:i+1,2), toolPath(i:i+1,3), 'Color', 'r', 'LineWidth', 2);
        else
            line(toolPath(i:i+1,1), toolPath(i:i+1,2), toolPath(i:i+1,3), 'Color', 'g', 'LineWidth', 2);
        end
    end
    title('Tool Path Colored by Collision Status (Red=Collision, Green=No Collision)');
    axis equal; grid on;
    xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
end

function waypoints = interpolateWaypoints(centerline, numPoints)
    % Calculate cumulative distances along the centerline
    segments = diff(centerline);
    segmentLengths = sqrt(sum(segments.^2, 2));
    cumDist = [0; cumsum(segmentLengths)];
    
    % Create new distance values for interpolation
    newDist = linspace(0, cumDist(end), numPoints);
    
    % Interpolate each coordinate
    waypoints = zeros(numPoints, 3);
    for dim = 1:3
        waypoints(:,dim) = interp1(cumDist, centerline(:,dim), newDist, 'pchip');
    end
end

function tangents = calculateTangents(waypoints)
    % Calculate tangent vectors along the path
    numPoints = size(waypoints, 1);
    tangents = zeros(numPoints, 3);
    
    % Calculate finite differences for interior points
    for i = 2:numPoints-1
        tangents(i,:) = waypoints(i+1,:) - waypoints(i-1,:);
        tangents(i,:) = tangents(i,:) / norm(tangents(i,:));
    end
    
    % For endpoints, use forward/backward differences
    tangents(1,:) = waypoints(2,:) - waypoints(1,:);
    tangents(1,:) = tangents(1,:) / norm(tangents(1,:));
    
    tangents(numPoints,:) = waypoints(numPoints,:) - waypoints(numPoints-1,:);
    tangents(numPoints,:) = tangents(numPoints,:) / norm(tangents(numPoints,:));
end

function collision = checkCollision(toolVertices, pipeVertices, pipeFaces)
    % Improved collision detection using point-to-triangle distance
    collision = false;
    
    % Check if any tool vertex is too close to the pipe wall
    for i = 1:size(toolVertices, 1)
        pt = toolVertices(i,:);
        
        % Find minimum distance to pipe vertices as a quick check
        distToVerts = sqrt(sum((pipeVertices - pt).^2, 2));
        minDistToVert = min(distToVerts);
        
        % If very close to a vertex, likely a collision
        if minDistToVert < 0.03
            collision = true;
            return;
        end
        
        % For points that are moderately close, check distance to faces
        if minDistToVert < 0.2
            % Check a subset of nearby faces to save computation
            [~, nearestVertIdx] = min(distToVerts);
            
            % Find faces that contain the nearest vertex
            faceIndices = find(any(pipeFaces == nearestVertIdx, 2));
            
            % Check distance to these faces
            for j = 1:length(faceIndices)
                faceIdx = faceIndices(j);
                v1 = pipeVertices(pipeFaces(faceIdx,1),:);
                v2 = pipeVertices(pipeFaces(faceIdx,2),:);
                v3 = pipeVertices(pipeFaces(faceIdx,3),:);
                
                % Calculate point-to-triangle distance
                distToFace = pointTriangleDistance(pt, v1, v2, v3);
                
                % If too close to face, we have a collision
                if distToFace < 0.05
                    collision = true;
                    return;
                end
            end
        end
    end
end

function dist = pointTriangleDistance(p, v1, v2, v3)
    % Calculate minimum distance from point p to triangle defined by vertices v1, v2, v3
    
    % Edge vectors
    e1 = v2 - v1;
    e2 = v3 - v1;
    
    % Normal vector
    n = cross(e1, e2);
    n = n / norm(n);
    
    % Distance from point to plane
    planeDistance = abs(dot(p - v1, n));
    
    % Project point onto the plane of the triangle
    pp = p - planeDistance * n;
    
    % Check if projected point is inside the triangle
    % Use barycentric coordinates
    area = 0.5 * norm(cross(e1, e2));
    
    area1 = 0.5 * norm(cross(v2 - pp, v3 - pp));
    area2 = 0.5 * norm(cross(v3 - pp, v1 - pp));
    area3 = 0.5 * norm(cross(v1 - pp, v2 - pp));
    
    sum_area = area1 + area2 + area3;
    
    if abs(sum_area - area) < 1e-10
        % Point projects inside the triangle
        dist = planeDistance;
    else
        % Point projects outside the triangle, find min distance to edges
        dist1 = pointLineSegmentDistance(p, v1, v2);
        dist2 = pointLineSegmentDistance(p, v2, v3);
        dist3 = pointLineSegmentDistance(p, v3, v1);
        dist = min([dist1, dist2, dist3]);
    end
end

function dist = pointLineSegmentDistance(p, v1, v2)
    % Calculate minimum distance from point p to line segment v1-v2
    
    % Line segment vector
    d = v2 - v1;
    lengthSquared = sum(d.^2);
    
    if lengthSquared == 0
        % v1 and v2 are the same point
        dist = norm(p - v1);
        return;
    end
    
    % Calculate projection parameter
    t = dot(p - v1, d) / lengthSquared;
    
    if t < 0
        % Closest to v1
        dist = norm(p - v1);
    elseif t > 1
        % Closest to v2
        dist = norm(p - v2);
    else
        % Closest to point on segment
        proj = v1 + t * d;
        dist = norm(p - proj);
    end
end

% Call the main function
% simulatePipeMotion('Pipe for matlab code.STL');