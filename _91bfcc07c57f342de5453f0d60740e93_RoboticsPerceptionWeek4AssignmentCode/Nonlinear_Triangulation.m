function X = Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
    %% Nonlinear_Triangulation
    % Refining the poses of the cameras to get a better estimate of the points
    % 3D position
    % Inputs: 
    %     K - size (3 x 3) camera calibration (intrinsics) matrix
    %     C1, R1 - the first camera pose
    %     C2, R2 - the second camera pose
    %     C3, R3 - the third camera pose
    %     x1, x2, x3 - size (N x 2) matrices of points in images 1, 2, and 3
    %     X0 - size (N x 3) initial estimate of 3D points
    % Outputs: 
    %     X - size (N x 3) matrix of refined point 3D locations

    % Number of points
    N = size(X0, 1);

    % Initialize refined 3D points with the initial estimate
    X = X0;

    % Set the number of iterations for the optimization
    %num_iterations = 30;
    delta_X = inf;

    while delta_X > 0.001
        % Initialize the Jacobian and error vector
        J = zeros(6 * N, 3 * N);
        b = zeros(6 * N, 1);
        
        for i = 1:N
            % Get the current estimate of the 3D point
            X_i = X(i, :)';

            % Compute the reprojection errors and Jacobians for each camera
            [e1, J1] = reprojectionErrorAndJacobian(K, C1, R1, x1(i, :)', X_i);
            [e2, J2] = reprojectionErrorAndJacobian(K, C2, R2, x2(i, :)', X_i);
            [e3, J3] = reprojectionErrorAndJacobian(K, C3, R3, x3(i, :)', X_i);

            % Populate the global Jacobian and error vector
            J((i-1)*6+1:i*6, (i-1)*3+1:i*3) = [J1; J2; J3];
            b((i-1)*6+1:i*6) = [e1; e2; e3];
        end

        % Compute the update step using the normal equations
        delta_X = (J' * J) \ (J' * b);

        % Update the 3D points
        X = X + reshape(delta_X, 3, N)';
    end
end

function [error, J] = reprojectionErrorAndJacobian(K, C, R, x, X)
    % Project the 3D point into the image
    X_cam = R * (X - C);
    x_proj = K * X_cam;
    w = x_proj(3);
    x_proj = x_proj(1:2) / x_proj(3);

    % Compute the reprojection error
    error = x - x_proj;

    % Intrinsic parameters
    f = K(1, 1);
    px = K(1, 3);
    py = K(2, 3);
    
    % Rotation matrix elements
    r11 = R(1, 1); r12 = R(1, 2); r13 = R(1, 3);
    r21 = R(2, 1); r22 = R(2, 2); r23 = R(2, 3);
    r31 = R(3, 1); r32 = R(3, 2); r33 = R(3, 3);
    
    % Jacobian components
    %w = X_cam(3);
    du_dX = [f*r11 + px*r31, f*r12 + px*r32, f*r13 + px*r33];
    dv_dX = [f*r21 + py*r31, f*r22 + py*r32, f*r23 + py*r33];
    dw_dX = [r31, r32, r33];
    
    % Complete Jacobian matrix
    J = [
        (w * du_dX - x_proj(1) * dw_dX) / w^2;
        (w * dv_dX - x_proj(2) * dw_dX) / w^2
    ];
end