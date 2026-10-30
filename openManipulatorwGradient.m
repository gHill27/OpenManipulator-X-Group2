function [gradw] = openManipulatorwGradient(q, params)
    % Ensure q is a 4x1 column vector regardless of whether a row or column was passed
    q = q(:);
    n = length(q);

    % 1. Call your existing analytical Jacobian function (3x4)
    J = openManipulatorJacobian(q, params);

    % 2. Compute Analytical Hessian Tensor H (3x4x4) where H(:,:,i) = dJ/dq_i
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    a1 = params.a1; a2 = params.a2; d2 = params.d2; le = params.le;

    s1 = sin(q1); c1 = cos(q1);
    s2 = sin(q2); c2 = cos(q2);
    s23 = sin(q2 + q3); c23 = cos(q2 + q3);
    s234 = sin(q2 + q3 + q4); c234 = cos(q2 + q3 + q4);

    % Radial plane reach r and its first derivatives
    r      =  a1*c2 + d2*s2 + a2*c23 + le*c234;
    dr_dq2 = -a1*s2 + d2*c2 - a2*s23 - le*s234;
    dr_dq3 = -a2*s23 - le*s234;
    dr_dq4 = -le*s234;

    % Second partial derivatives of r
    d2r_dq22  = -a1*c2 - d2*s2 - a2*c23 - le*c234;
    d2r_dq2q3 = -a2*c23 - le*c234;
    d2r_dq2q4 = -le*c234;
    d2r_dq33  = -a2*c23 - le*c234;
    d2r_dq3q4 = -le*c234;
    d2r_dq44  = -le*c234;

    % Second partial derivatives of z (row 3 of Jacobian)
    d2z_dq22  =  a1*s2 - d2*c2 + a2*s23 + le*s234;
    d2z_dq2q3 =  a2*s23 + le*s234;
    d2z_dq2q4 =  le*s234;
    d2z_dq33  =  a2*s23 + le*s234;
    d2z_dq3q4 =  le*s234;
    d2z_dq44  =  le*s234;

    H = zeros(3, 4, 4);

    % H1 = dJ / dq1 (Notice this is all zeros in row 3, meaning base rotation q1 
    % has zero effect on manipulability!)
    H(:, :, 1) = [ -c1*r,  -s1*dr_dq2,  -s1*dr_dq3,  -s1*dr_dq4;
                   -s1*r,   c1*dr_dq2,   c1*dr_dq3,   c1*dr_dq4;
                       0,           0,           0,           0 ];

    % H2 = dJ / dq2
    H(:, :, 2) = [ -s1*dr_dq2,   c1*d2r_dq22,   c1*d2r_dq2q3,   c1*d2r_dq2q4;
                    c1*dr_dq2,   s1*d2r_dq22,   s1*d2r_dq2q3,   s1*d2r_dq2q4;
                            0,      d2z_dq22,      d2z_dq2q3,      d2z_dq2q4 ];

    % H3 = dJ / dq3
    H(:, :, 3) = [ -s1*dr_dq3,   c1*d2r_dq2q3,   c1*d2r_dq33,   c1*d2r_dq3q4;
                    c1*dr_dq3,   s1*d2r_dq2q3,   s1*d2r_dq33,   s1*d2r_dq3q4;
                            0,      d2z_dq2q3,      d2z_dq33,      d2z_dq3q4 ];

    % H4 = dJ / dq4
    H(:, :, 4) = [ -s1*dr_dq4,   c1*d2r_dq2q4,   c1*d2r_dq3q4,   c1*d2r_dq44;
                    c1*dr_dq4,   s1*d2r_dq2q4,   s1*d2r_dq3q4,   s1*d2r_dq44;
                            0,      d2z_dq2q4,      d2z_dq3q4,      d2z_dq44 ];

    % 3. Compute Manipulability m = sqrt(det(J * J'))
    JJT = J * J';
    det_JJT = det(JJT);

    gradw = zeros(n, 1);
    if det_JJT <= 1e-12
        return; % Singularity guard
    end

    m = sqrt(det_JJT);

    % Invert (J * J') with light damping for numerical stability near singularities
    inv_JJT = (JJT + 1e-6 * eye(3)) \ eye(3);
    vec_inv_JJT = inv_JJT(:); % vec((J*J')^-1) -> 9x1 column vector

    % 4. Evaluate Equation (11) for each joint i
    for i = 1:n
        Hi = H(:, :, i);
        J_HiT = J * Hi';
        vec_J_HiT = J_HiT(:); % vec(J * H_i^T) -> 9x1 column vector
        
        gradw(i) = m * (vec_J_HiT' * vec_inv_JJT);
    end
end