function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose transation
%     R - size (3 x 1) pose rotation
%
% IMPORTANT NOTE: While theoretically you can use the x directly when solving
% for the P = [R t] matrix then use the K matrix to correct the error, this is
% more numeically unstable, and thus it is better to calibrate the x values
% before the computation of P then extract R and t directly

    % Number of points
    N = size(X, 1);
    
    % Convert 2D points to homogeneous coordinates
    x_homogeneous = [x, ones(N, 1)];
    
    % Normalize the 2D points using the intrinsic matrix K
    x_normalized = (K \ x_homogeneous')';

    % Construct the matrix A
    A = [];
%     for i = 1:N
%        X_i = X(i, :);
%        x_i = x_normalized(i, 1);
%        y_i = x_normalized(i, 2);
%        A(2*i-1, :) = [-X_i, -1, 0, 0, 0, 0, x_i*X_i, x_i];
%        A(2*i, :)   = [0, 0, 0, 0, -X_i, -1, y_i*X_i, y_i];
%     end
    for i=1:N 

        Xt = [X(i,:), 1]; 

        A = [A; Vec2Skew(x_normalized(i,:)) * [Xt, zeros(1,4), zeros(1,4); zeros(1,4), Xt, zeros(1,4); zeros(1,4), zeros(1,4), Xt]]; 

    end 
    
    % Solve the linear system using SVD
    [~, ~, V] = svd(A);
    P = reshape(V(:, end), 4, 3)'; % Projection matrix P is the last column of V reshaped

    % Extract R and t from P
    R = P(:, 1:3);
    t = P(:, 4);
    
    % Ensure R is a proper rotation matrix by enforcing orthogonality
    [U, D, V] = svd(R);
    if int8(det(U*V')) == 1
        R = U * V';
        t = t/D(1,1);
        
        %keyboard
    elseif int8(det(U*V')) == -1
        R = -U * V';
        t = -t/D(1,1);    
        %keyboard
    end
    

    % Translation vector C is -R' * t
    C = -R' * t;
end




