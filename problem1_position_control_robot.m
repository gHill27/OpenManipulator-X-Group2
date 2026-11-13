clc, clear all, close all;

filename = 'Data/problem1.mat';
load('Data/problem1sym.mat');

%% Robot setup:
robot = Robot();
homingtravelTime = 1;
robot.writeTime(homingtravelTime)
robot.writeMode('p') %sets it to position command.
robot.writeJoints([0,0,0,0]) %go home
pause(1.0)
robot.writeMode('v') %sets it to velocity command


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

timestamps = zeros(1, 1000);
qHist = zeros(numel(q), 1000);
qdotHist = zeros(numel(q), 1000);
xHist = zeros(3, 1000);
errHist = zeros(3, 1000);
manipHist = zeros(1, 1000);
real_q_dot_hist = zeros(4, 1000);

v_des = [0; 0; 0];
tic; % Start timer
travelTime = 5.0;
index = 1; %matlab starts at 1 so keep consistent
try 
    while toc < travelTime
        joint_matrix = robot.getJointsReadings();
        real_q_dot_hist(:,index) = deg2rad(joint_matrix(2,:));
        q = joint_matrix(1,:); % in degrees
        x = openManipulatorPosition(deg2rad(q'),params);
        q_dot = ControLaw(x_des,v_des, x, deg2rad(q), params, Kp);
        robot.writeVelocities(rad2deg(q_dot));
        
        % Save data
        qHist(:, index) = deg2rad(q);
        xHist(:, index) = x;
        errHist(:, index) = x_des - x;
        qdotHist(:, index) = q_dot;
        
        J = openManipulatorJacobian(deg2rad(q), params);
        manipHist(:,index) = sqrt(det(J*J'));
        timestamps(index) = toc;
        index = index + 1;
        pause(0.005);
    end
    robot.writeVelocities([0; 0; 0; 0])
catch ME
    fprintf('error stopping: %s\n', ME.message);
    robot.writeVelocities([0; 0; 0; 0])
end

%% Plot variables
last_non_zero_index = find(timestamps, 1, 'last');

% Use the cleaner logic to find the shortest array length
if last_non_zero_index > width(symQdotHist)
    shortest = width(symQdotHist);
else
    shortest = last_non_zero_index;
end

figure('Name', 'Controller Results');

% Subplot 1: Joint Positions
subplot(2,2,1);
plot(timestamps(1, 1:shortest), qHist(:, 1:shortest), 'LineWidth', 1.5);
hold on;
ax = gca; ax.ColorOrderIndex = 1; % Resets color order so sim matches real colors
plot(timestamps(1, 1:shortest), symqHist(:, 1:shortest), ':', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('q (rad)');
title('Joint Positions vs. Time');
legend('q1 real', 'q2 real', 'q3 real', 'q4 real', 'q1 sim', 'q2 sim', 'q3 sim', 'q4 sim');
grid on;

% Subplot 2: Joint Velocities
subplot(2,2,2);
plot(timestamps(1, 1:shortest), qdotHist(:, 1:shortest), 'LineWidth', 1.5);
hold on;
ax = gca; ax.ColorOrderIndex = 1;
plot(timestamps(1, 1:shortest), symQdotHist(:, 1:shortest), ':', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('qdot (rad/s)');
title('Joint Velocities vs. Time');
legend('qdot1 ctrl', 'qdot2 ctrl', 'qdot3 ctrl', 'qdot4 ctrl', 'qdot1 sim', 'qdot2 sim', 'qdot3 sim', 'qdot4 sim');
grid on;

% Subplot 3: End-Effector Position
subplot(2,2,3);
plot(timestamps(1, 1:shortest), xHist(:, 1:shortest), 'LineWidth', 1.5);
hold on;
ax = gca; ax.ColorOrderIndex = 1;
plot(timestamps(1, 1:shortest), symxHist(:, 1:shortest), ':', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('End-Effector Position (m)');
title('End-Effector Position vs. Time');
legend('x real', 'y real', 'z real', 'x sim', 'y sim', 'z sim');
grid on;

% Subplot 4: Position Error Norm
subplot(2,2,4);
errHistNorm = vecnorm(errHist, 2, 1);
symErrHistNorm = vecnorm(symerrHist, 2, 1);
plot(timestamps(1, 1:shortest), errHistNorm(:, 1:shortest), 'LineWidth', 1.5);
hold on;
ax = gca; ax.ColorOrderIndex = 1;
plot(timestamps(1, 1:shortest), symErrHistNorm(:, 1:shortest), ':', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('Norm of Position Error (m)');
title('Position Error Norm vs. Time');
legend('Real Error', 'Sim Error');
grid on;

figure();

% Joint 1
subplot(2,2,1);
plot(timestamps(1,1:shortest), qdotHist(1, 1:shortest), 'LineWidth', 1.5); % controller
hold on;
plot(timestamps(1,1:shortest), real_q_dot_hist(1, 1:shortest), 'LineWidth', 1.5); % real
plot(timestamps(1,1:shortest), symQdotHist(1, 1:shortest), ':', 'LineWidth', 1.5); % sim
hold off;
xlabel('Time (s)');
ylabel('Qdot Angular Velocity');
title('Joint 1 Real and Controller Qdot vs. Time');
legend('controller', 'real', 'sim');
grid on;

% Joint 2
subplot(2,2,2);
plot(timestamps(1,1:shortest), qdotHist(2, 1:shortest), 'LineWidth', 1.5);
hold on;
plot(timestamps(1,1:shortest), real_q_dot_hist(2, 1:shortest), 'LineWidth', 1.5);
plot(timestamps(1,1:shortest), symQdotHist(2, 1:shortest), ':', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('Qdot Angular Velocity');
title('Joint 2 Real and Controller Qdot vs. Time');
legend('controller', 'real', 'sim');
grid on;

% Joint 3
subplot(2,2,3);
plot(timestamps(1,1:shortest), qdotHist(3, 1:shortest), 'LineWidth', 1.5);
hold on;
plot(timestamps(1,1:shortest), real_q_dot_hist(3, 1:shortest), 'LineWidth', 1.5);
plot(timestamps(1,1:shortest), symQdotHist(3, 1:shortest), ':', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('Qdot Angular Velocity');
title('Joint 3 Real and Controller Qdot vs. Time');
legend('controller', 'real', 'sim');
grid on;

% Joint 4
subplot(2,2,4);
plot(timestamps(1,1:shortest), qdotHist(4, 1:shortest), 'LineWidth', 1.5);
hold on;
plot(timestamps(1,1:shortest), real_q_dot_hist(4, 1:shortest), 'LineWidth', 1.5);
plot(timestamps(1,1:shortest), symQdotHist(4, 1:shortest), ':', 'LineWidth', 1.5);
hold off;
xlabel('Time (s)');
ylabel('Qdot Angular Velocity');
title('Joint 4 Real and Controller Qdot vs. Time');
legend('controller', 'real', 'sim');
grid on;

save(filename, "qHist", "qdotHist", "real_q_dot_hist", "timestamps")
function [qdot] = ControLaw(desired_pos,desired_velocity,curr_pos, q, params, Kp)
    J = openManipulatorJacobian(q, params);
    e = desired_pos(:) - curr_pos(:) ;
    qdot = pinv(J)*(desired_velocity(:) + Kp*e);

end