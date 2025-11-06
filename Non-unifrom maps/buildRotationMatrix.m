function [rotMatrix] = buildRotationMatrix(phi1, Phi, phi2)
%     % Rotation about Z-axis (phi1)
%     R_z_phi1 = [cos(phi1), -sin(phi1), 0;
%                 +sin(phi1), cos(phi1),  0;
%                 0,         0,          1];
% 
%     % Rotation about X-axis (Phi)
%     R_x_Phi = [1,  0,          0;
%                0,  cos(Phi), -sin(Phi);
%                0,  +sin(Phi),  cos(Phi)];
% 
%     % Rotation about Z-axis (phi2)
%     R_z_phi2 = [cos(phi2), -sin(phi2), 0;
%                 +sin(phi2), cos(phi2),  0;
%                 0,         0,          1];
% 
%     % Compute final rotation matrix
%     rotMatrix = R_z_phi1 * R_x_Phi * R_z_phi2;

    % Rotation about Z-axis (phi1)
    R_z_phi1 = [cos(phi1), +sin(phi1), 0;
                -sin(phi1), cos(phi1),  0;
                0,         0,          1];

    % Rotation about X-axis (Phi)
    R_x_Phi = [1,  0,          0;
               0,  cos(Phi), +sin(Phi);
               0,  -sin(Phi),  cos(Phi)];

    % Rotation about Z-axis (phi2)
    R_z_phi2 = [cos(phi2), +sin(phi2), 0;
                -sin(phi2), cos(phi2),  0;
                0,         0,          1];

    % Compute final rotation matrix
    gmat = R_z_phi2 * R_x_Phi * R_z_phi1;
    rotMatrix = inv(gmat);

end