function q_dot = computeRedundantControl(q, x_des, x_dot_des, x_current, K, K_bar, params)
% COMPUTEREDUNDANTCONTROL Computes joint velocities for a 4-DOF redundant manipulator
% using the exact analytical Hessian tensor (Eq. 11) for manipulability gradient ascent.
%
% Inputs:
%   q         - Current joint angles (4x1 or 1x4 vector)
%   x_des     - Desired end-effector position (3x1 vector)
%   x_dot_des - Desired end-effector velocity (3x1 vector)
%   x_current - Current end-effector position (3x1 vector)
%   K         - Task-space position gain (3x3 matrix or scalar)
%   K_bar     - Null-space manipulability gain (4x4 matrix or scalar)
%   params    - Struct with fields: a1, a2, d2, le (and optional d0)
%
% Output:
%   q_dot     - Commanded joint velocities (4x1 column vector)

    % Enforce column vectors to prevent dimension/bounds errors
    q = q(:);
    x_des = x_des(:);
    x_dot_des = x_dot_des(:);
    x_current = x_current(:);
    n = length(q); % n = 4

    % 1. Compute Analytical Jacobian J (3x4) and Hessian Tensor H (3x4x4)
    [J, H] = getJacobianAndHessian(q, params);

    % 2. Compute Yoshikawa Manipulability Index: m = sqrt(det(J * J'))
    JJT = J * J';
    det_JJT = det(JJT);
    
    % Guard against numerical underflow or exact singularity
    if det_JJT <= 1e-12
        m = 0;
         nabla_q_w = zeros(n, 1);
    else
        m = sqrt(det_JJT);
        
        % Damped inverse of (J*J') for numerical stability near singularities
        m_rows = size(J, 1);
        lambda_m = 1e-4; 
        inv_JJT = (JJT + lambda_m * eye(m_rows)) \ eye(m_rows);
        
        % Vectorize inv(J*J'): vec((J*J')^-1) -> 9x1 column vector
        vec_inv_JJT = inv_JJT(:);
        
        % 3. Compute Manipulability Jacobian J_m (4x1) via Equation (11)
        nabla_q_w = zeros(n, 1);
        for i = 1:n
            Hi = H(:, :, i);          % 3x4 Hessian slice: dJ / dq_i
            J_HiT = J * Hi';          % 3x3 matrix: J * H_i^T
            vec_J_HiT = J_HiT(:);     % vec(J * H_i^T) -> 9x1 column vector
            
            % Eq (11): m * vec(J * H_i^T)^T * vec((J * J^T)^-1)
            nabla_q_w(i) = m * (vec_J_HiT' * vec_inv_JJT);
        end
    end

    % 4. Compute Moore-Penrose Pseudoinverse (with small damping for safety)
    lambda_pinv = 1e-3;
    J_dagger = J' / (JJT + (lambda_pinv^2) * eye(size(J, 1)));

    % 5. Null-Space Projector: P = I - J_dagger * J (4x4)
    P = eye(n) - (J_dagger * J);

    % 6. Primary Task Control Law: J_dagger * (x_dot_des + K * (x_des - x))
    nu = x_dot_des + K * (x_des - x_current);
    q_dot_task = J_dagger * nu;

    % 7. Secondary Task (Null-Space Manipulability Maximization)
    q_dot_null = P * (K_bar * nabla_q_w);

    % 8. Total Joint Velocity Command (4x1)
    q_dot = q_dot_task + q_dot_null;
end