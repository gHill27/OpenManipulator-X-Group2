clc, clear all, close all;

filename = 'Data/problem1.mat';
load('Data/problem1sym.mat');

%% Robot setup:
robot = Robot();
homingtravelTime = 2;
robot.writeTime(homingtravelTime)
robot.writeMode('p') %sets it to position command.
robot.writeJoints([0,0,0,0]) %go home
pause(2.0)
robot.writeMode('v') %sets it to velocity command
robot.writeVelocities([0; 0; 0; 0]);
pause(0.15);


%% Geometric parameters
params.a0 = 0.012;
params.d1 = 0.0595;
params.a1 = 0.024;
params.d2 = 0.128;
params.a2 = 0.124;
params.le = 0.12;

%% Desired Location of the end effector
x_des = [0.20; 0.20; 0.15];    
Kp = 1.0;
q = [0; 0; 0; 0];

timestamps = zeros(1, 400);
qHist = zeros(numel(q), 400);
qdotHist = zeros(numel(q), 400);
xHist = zeros(3, 400);
errHist = zeros(3, 400);
manipHist = zeros(1, 400);
real_q_dot_hist = zeros(4, 400);

v_des = [0; 0; 0];
tic; % Start timer
travelTime = 10.0;
index = 1; %matlab starts at 1 so keep consistent

try 
    while toc < travelTime
        timestamps(index) = toc;

        joint_matrix = robot.getJointsReadings();
        real_q_dot_hist(:,index) = deg2rad(joint_matrix(2,:));
        q = deg2rad(joint_matrix(1,:));
        x = openManipulatorPosition(q, params);
        J = openManipulatorJacobian(q, params);
        
        q_dot = ControLaw(x_des,v_des, x, q, params, Kp);
        robot.writeVelocities(rad2deg(q_dot));
        
        % Save data
        qHist(:, index) = q;
        xHist(:, index) = x;
        errHist(:, index) = x_des - x;
        qdotHist(:, index) = q_dot;
        manipHist(:,index) = sqrt(det(J*J'));

        index = index + 1;
        pause(.05);
    end
    robot.writeVelocities([0; 0; 0; 0]);
catch ME
    fprintf('error stopping: %s\n', ME.message);
    robot.writeVelocities([0; 0; 0; 0]);
end

%% Plot variables
last_non_zero_index = find(timestamps, 1, 'last');



% Subplot 1: Joint Positions
figure();
plot(timestamps(1, 1:last_non_zero_index), qHist(:, 1:last_non_zero_index), 'LineWidth', 1.5);
hold on;
ax = gca; ax.ColorOrderIndex = 1; % Resets color order so sim matches real colors
plot(t, symqHist, ':', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('q (rad)');
title('Joint Positions vs. Time');
legend('Real q1', 'Real q2', 'Real q3', 'Real q4', 'Sim q1', 'Sim q2', 'Sim q3', 'Sim q4');
grid on;

% % Subplot 2: Joint Velocities
% subplot(2,2,2);
% plot(timestamps(1, 1:last_non_zero_index), qdotHist(:, 1:last_non_zero_index), 'LineWidth', 1.5);
% hold on;
% ax = gca; ax.ColorOrderIndex = 1;
% plot(t, symQdotHist, ':', 'LineWidth', 1.5);
% hold off;
% xlabel('Time (s)');
% ylabel('qdot (rad/s)');
% title('Joint Velocities vs. Time');
% legend('qdot1 ctrl', 'qdot2 ctrl', 'qdot3 ctrl', 'qdot4 ctrl', 'qdot1 sim', 'qdot2 sim', 'qdot3 sim', 'qdot4 sim');
% grid on;

% Subplot 3: End-Effector Position
figure();
plot(timestamps(1, 1:last_non_zero_index), xHist(:, 1:last_non_zero_index), 'LineWidth', 1.5);
hold on;
ax = gca; ax.ColorOrderIndex = 1;
plot(t, symxHist, ':', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('End-Effector Position (m)');
title('End-Effector Position vs. Time');
legend('Real x', 'Real y', 'Real z', 'Sim x', 'Sim y', 'Sim z');
grid on;

% Subplot 4: Position Error Norm
figure();
errHistNorm = vecnorm(errHist, 2, 1);
symErrHistNorm = vecnorm(symerrHist, 2, 1);
plot(timestamps(1, 1:last_non_zero_index), errHistNorm(:, 1:last_non_zero_index), 'LineWidth', 1.5);
hold on;
ax = gca; ax.ColorOrderIndex = 1;
plot(t, symErrHistNorm, ':', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('Norm of Position Error (m)');
title('Position Error Norm vs. Time');
legend('Real Error', 'Sim Error');
grid on;


figure();

% Joint 1
subplot(2,2,1);
plot(timestamps(1,1:last_non_zero_index), qdotHist(1, 1:last_non_zero_index), 'LineWidth', 1.5); % controller
hold on;
plot(timestamps(1,1:last_non_zero_index), real_q_dot_hist(1, 1:last_non_zero_index), 'LineWidth', 1.5); % real
plot(t, symQdotHist(1, :), ':', 'LineWidth', 1.5); % sim
hold off;
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
title('Real, Commanded, and Sim Angular Velocity of Joint 1 vs. Time');
legend('Controller', 'Real', 'Sim');
grid on;

% Joint 2
subplot(2,2,2);
plot(timestamps(1,1:last_non_zero_index), qdotHist(2, 1:last_non_zero_index), 'LineWidth', 1.5);
hold on;
plot(timestamps(1,1:last_non_zero_index), real_q_dot_hist(2, 1:last_non_zero_index), 'LineWidth', 1.5);
plot(t, symQdotHist(2, :), ':', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('Angular Velocity');
title('Real, Commanded, and Sim Angular Velocity of Joint 2 vs. Time');
legend('Controller', 'Real', 'Sim');
grid on;

% Joint 3
subplot(2,2,3);
plot(timestamps(1,1:last_non_zero_index), qdotHist(3, 1:last_non_zero_index), 'LineWidth', 1.5);
hold on;
plot(timestamps(1,1:last_non_zero_index), real_q_dot_hist(3, 1:last_non_zero_index), 'LineWidth', 1.5);
plot(t, symQdotHist(3, :), ':', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
title('Real, Commanded, and Sim Angular Velocity of Joint 3 vs. Time');
legend('Controller', 'Real', 'Sim');
grid on;

% Joint 4
subplot(2,2,4);
plot(timestamps(1,1:last_non_zero_index), qdotHist(4, 1:last_non_zero_index), 'LineWidth', 1.5);
hold on;
plot(timestamps(1,1:last_non_zero_index), real_q_dot_hist(4, 1:last_non_zero_index), 'LineWidth', 1.5);
plot(t, symQdotHist(4, :), ':', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
title('Real, Commanded, and Sim Angular Velocity of Joint 4 vs. Time');
legend('Controller', 'Real', 'Sim');
grid on;

% Manip
figure();
plot(timestamps(1,1:last_non_zero_index), manipHist(:, 1:last_non_zero_index), 'LineWidth', 1.5);
hold on;
ax = gca; ax.ColorOrderIndex = 1;
plot(t, symmanipHist, ':', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('Translational Manipulability');
title('Translational Manipulability vs. Time');
legend('Real', 'Sim');
ylim([0.0078, .0112])
grid on;


save(filename, "qHist", "qdotHist", "real_q_dot_hist", "timestamps")
function [qdot] = ControLaw(desired_pos,desired_velocity,curr_pos, q, params, Kp)
    J = openManipulatorJacobian(q, params);
    e = desired_pos(:) - curr_pos(:) ;
    qdot = pinv(J)*(desired_velocity(:) + Kp*e);

end