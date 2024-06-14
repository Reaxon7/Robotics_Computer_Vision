function X = LinearTriangulation(K, C1, R1, C2, R2, x1, x2)
    %% LinearTriangulation
    % Find 3D positions of the point correspondences using the relative
    % position of one camera from another
    % Inputs:
    %     K - size (3 x 3) intrinsic matrix
    %     C1 - size (3 x 1) translation of the first camera pose
    %     R1 - size (3 x 3) rotation of the first camera pose
    %     C2 - size (3 x 1) translation of the second camera
    %     R2 - size (3 x 3) rotation of the second camera pose
    %     x1 - size (N x 2) matrix of points in image 1
    %     x2 - size (N x 2) matrix of points in image 2, each row corresponding
    %       to x1
    % Outputs: 
    %     C - the correct camera pose translation
    %     R - the correct camera pose rotation
    %     X - size (N x 3) matrix whose rows represent the 3D triangulated points

    % Number of points
    N = size(x1, 1);
    
    % Initialize the matrix to store the 3D points for each configuration
    Xset = cell(4, 1);
    
    % Four possible configurations of camera centers and rotations
    Cset = {C1, -C1, C2, -C2};
    Rset = {R1, R1, R2, R2};
    
    % Compute the projection matrices for both cameras
    %P1 = K * [R1, -R1 * C1];
    %P2 = K * [R2, -R2 * C2];
    
    % Loop over all configurations
    for configIdx = 1:4
        P1 = K * [R1, -R1 * C1];
        P2 = K * [Rset{configIdx}, -Rset{configIdx} * Cset{configIdx}];
        X = zeros(N, 3);
        for i = 1:N
            % Extract the corresponding points
            x1_i = [x1(i, :), 1]';
            x2_i = [x2(i, :), 1]';
            
            % Formulate the linear system
            A = [
                Vec2Skew(x1_i) * P1;
                Vec2Skew(x2_i) * P2
            ];
            
            % Solve the linear system using SVD
            [~, ~, V] = svd(A);
            X_homogeneous = V(:, end);
            
            % Convert from homogeneous coordinates to 3D coordinates
            X(i, :) = X_homogeneous(1:3)' / X_homogeneous(4);
        end
        Xset{configIdx} = X;
    end
    
    % Disambiguate the camera pose
    [~, ~, X] = DisambiguateCameraPose(Cset, Rset, Xset);
end

function [C, R, X0] = DisambiguateCameraPose(Cset, Rset, Xset)
    %% DisambiguateCameraPose
    % Find the unique camera pose by checking the cheirality condition.
    % Inputs:
    %     Cset - 4 configurations of camera centers
    %     Rset - 4 configurations of camera rotations
    %     Xset - 4 sets of triangulated points
    % Outputs:
    %     C - the correct camera pose translation
    %     R - the correct camera pose rotation
    %     X0 - the 3D triangulated points from the correct camera pose

    numConfigurations = length(Cset);
    maxPositiveDepths = 0;
    bestConfigIdx = 1;

    for i = 1:numConfigurations
        C = Cset{i};
        R = Rset{i};
        X = Xset{i};

        % Calculate the number of points in front of the camera
        numPositiveDepths = sum(checkCheirality(C, R, X));

        if numPositiveDepths > maxPositiveDepths
            maxPositiveDepths = numPositiveDepths;
            bestConfigIdx = i;
        end
    end

    % Return the best configuration
    C = Cset{bestConfigIdx};
    R = Rset{bestConfigIdx};
    X0 = Xset{bestConfigIdx};
    
end

function isPositive = checkCheirality(C, R, X)
    %% checkCheirality
    % Check if the 3D points are in front of the camera
    % Inputs:
    %     C - camera center
    %     R - camera rotation matrix
    %     X - 3D points
    % Outputs:
    %     isPositive - logical array indicating if points are in front of the camera

    X_cam = R * (X' - C); % Transform points to camera coordinate system
    X_cam = X_cam'; % Transpose to match dimensions
    isPositive = X_cam(:, 3) > 0; % Check if the Z coordinate is positive
    
end