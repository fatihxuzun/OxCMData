function OxCM_rotateProcessedData(rotationAngle, filename)
    % STEP 2
    % OxCM_rotateProcessedData rotates the in-plane (XY) coordinates of a 2.5D profilometry 
    % point cloud about the Z-axis, re-centers the data at the origin, and updates 
    % the text file for subsequent processing and solid meshing.
    %
    % Inputs:
    %   - rotationAngle: In-plane rotation angle about the Z-axis in degrees (default: 0)
    %   - filename: Name of text file containing [x, y, z] data to rotate and update (default: 'myData.txt')

    clc;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameter Defaults & Validation
    if nargin < 1 || isempty(rotationAngle)
        rotationAngle = 0;
    end
    if nargin < 2 || isempty(filename)
        filename = 'myData.txt';
    end

    if ~isfile(filename)
        error('File "%s" not found. Please run Step 1 (OxCM_profilometryData) first or specify a valid file.', filename);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read data file
    raw = importdata(filename, '\t');
    if isstruct(raw)
        data = raw.data';
    else
        data = raw';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply 2D rotation in the XY plane
    theta = deg2rad(rotationAngle);
    R = [cos(theta), -sin(theta); 
         sin(theta),  cos(theta)];
    
    data(1:2, :) = R * data(1:2, :);

    % Re-center data in the XY plane
    data(1,:) = data(1,:) - mean(data(1,:));
    data(2,:) = data(2,:) - mean(data(2,:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot rotated data
    d = mean(abs(data(3,:))) * 5 / ((mean(abs(data(1,:))) + mean(abs(data(2,:)))) / 2);
    
    figure(7); 
    scatter3(data(1,:), data(2,:), data(3,:), 1, '.', 'r'); 
    shading interp; view(0, 90); grid on; axis equal; axis tight; daspect([1, 1, d]); 
    title(sprintf('Rotated (\\theta = %.2f^\\circ)', rotationAngle));
    xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis');
    fontsize(16, 'points'); fontname('Arial');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Write updated point cloud to file
    myData = single(data(1:3, :)');
    writematrix(myData, filename, 'Delimiter', '\t');
    fprintf('Successfully rotated (%.2f deg) and updated "%s".\n', rotationAngle, filename);
end