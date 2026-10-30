clc, clear all, close all;


syms q1 q2 q3 q4;
q_sym = [q1; q2; q3; q4];

% Geometric parameters
params.a0 = 0.012;
params.d1 = 0.0595;
params.a1 = 0.024;
params.d2 = 0.128;
params.a2 = 0.124;
params.le = 0.12;

jacob_sym = openManipulatorJacobian(q_sym, params);
w = sqrt(det(jacob_sym*jacob_sym'));
gradw_sym = gradient(w, q_sym);

%gradw = double(subs(gradw_sym, q_sym, q));

matlabFunction(gradw_sym, "File", "openManipulatorwGradient", "Vars", {q_sym});