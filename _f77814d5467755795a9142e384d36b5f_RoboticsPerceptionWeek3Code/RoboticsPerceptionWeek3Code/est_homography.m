function [ H ] = est_homography(video_pts, logo_pts)
% est_homography estimates the homography to transform each of the
% video_pts into the logo_pts
% Inputs:
%     video_pts: a 4x2 matrix of corner points in the video
%     logo_pts: a 4x2 matrix of logo points that correspond to video_pts
% Outputs:
%     H: a 3x3 homography matrix such that logo_pts ~ H*video_pts
% Written for the University of Pennsylvania's Robotics:Perception course

% YOUR CODE HERE

A = [];
for n = 1:size(video_pts,1)
    A = [A;-video_pts(n,1),-video_pts(n,2),-1,0,0,0,video_pts(n,1)*logo_pts(n,1),video_pts(n,2)*logo_pts(n,1),logo_pts(n,1)];
    A = [A;0,0,0,-video_pts(n,1),-video_pts(n,2),-1,video_pts(n,1)*logo_pts(n,2),video_pts(n,2)*logo_pts(n,2),logo_pts(n,2)];
end


[~, ~, V] = svd(A);

H = reshape(V(:,end),3,3)';

end

