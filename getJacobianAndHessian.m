function [J, H] = getJacobianAndHessian(q, params)
    % Unpack joints and parameters
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    a1 = params.a1; a2 = params.a2; d2 = params.d2; le = params.le;

    % Precompute trig terms for speed and readability
    s1 = sin(q1); c1 = cos(q1);
    s2 = sin(q2); c2 = cos(q2);
    s23 = sin(q2 + q3); c23 = cos(q2 + q3);
    s234 = sin(q2 + q3 + q4); c234 = cos(q2 + q3 + q4);

    % Radial reach in the arm plane (r) and its partial derivatives:
    % r = a1*c2 + d2*s2 + a2*c23 + le*c234
    r      =  a1*c2 + d2*s2 + a2*c23 + le*c234;
    dr_dq2 = -a1*s2 + d2*c2 - a2*s23 - le*s234;
    dr_dq3 = -a2*s23 - le*s234;
    dr_dq4 = -le*s234;

    % Second partial derivatives of r (used in Hessian slices H2, H3, H4)
    d2r_dq22  = -a1*c2 - d2*s2 - a2*c23 - le*c234; % = -r
    d2r_dq2q3 = -a2*c23 - le*c234;
    d2r_dq2q4 = -le*c234;
    d2r_dq33  = -a2*c23 - le*c234;
    d2r_dq3q4 = -le*c234;
    d2r_dq44  = -le*c234;

    % Vertical z-row terms (dz/dq_i) and their partial derivatives:
    % z = d0 - a1*s2 + d2*c2 - a2*s23 - le*s234
    dz_dq2 = -a1*c2 - d2*s2 - a2*c23 - le*c234;
    dz_dq3 = -a2*c23 - le*c234;
    dz_dq4 = -le*c234;

    d2z_dq22  =  a1*s2 - d2*c2 + a2*s23 + le*s234;
    d2z_dq2q3 =  a2*s23 + le*s234;
    d2z_dq2q4 =  le*s234;
    d2z_dq33  =  a2*s23 + le*s234;
    d2z_dq3q4 =  le*s234;
    d2z_dq44  =  le*s234;

    % --- Construct 3x4 Jacobian Matrix J ---
    J = [ -s1*r,   c1*dr_dq2,   c1*dr_dq3,   c1*dr_dq4;
           c1*r,   s1*dr_dq2,   s1*dr_dq3,   s1*dr_dq4;
              0,      dz_dq2,      dz_dq3,      dz_dq4 ];

    % --- Construct 3x4x4 Hessian Tensor H (where H(:,:,i) = dJ / dq_i) ---
    H = zeros(3, 4, 4);

    % H1 = dJ / dq1
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
end