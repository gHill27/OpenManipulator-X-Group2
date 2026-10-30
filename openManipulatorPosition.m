function [pos_vector] = openManipulatorPosition(q, params)
% pos and joint vector is a 1 x n matrix 
    % q = joint_vector;

    % params.a0 = 0.012;
    % params.a1 = 0.024;
    % params.a2 = 0.124;
    % params.d1 = 0.0595;
    % params.d2 = 0.128;
    % params.le = 0.12;

    pos_vector = [
        params.a0 + cos(q(1))*(params.a1*cos(q(2))+ params.d2*sin(q(2)) + params.a2*cos(q(2)+q(3)) + params.le*cos(q(2)+q(3)+q(4)));
        sin(q(1))*(params.a1*cos(q(2)) + params.d2*sin(q(2)) + params.a2*cos(q(2)+q(3)) + params.le*cos(q(2)+q(3)+q(4)));
        params.d1 - params.a1*sin(q(2)) + params.d2*cos(q(2)) - params.a2*(sin(q(2)+q(3))) - params.le*sin(q(2)+q(3)+q(4));
        ];
end