clc, clear all, close all;


filename = 'Data/problem2.mat';
load('Data/problem2sym.mat');


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
Kbar = 1000;
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
        q_dot = ControLaw(x_des,v_des, x, deg2rad(q), params, Kp, Kbar);
        robot.writeVelocities(rad2deg(q_dot));
        
        % Save data
        qHist(:, index) = q;
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
% 
%% Plot variables
last_non_zero_index = find(timestamps, 1, 'last');
figure('Name', 'Controller Results');
subplot(2,2,1);
plot(timestamps(1, 1:last_non_zero_index), qHist(:, 1:last_non_zero_index));
xlabel('Time (s)');
ylabel('q (rad)');
title('Joint Positions vs. Time');
legend('Joint 1', 'Joint 2', 'Joint 3', 'Joint 4');
grid on;

subplot(2,2,2);
plot(timestamps(1, 1:last_non_zero_index), qdotHist(:, 1:last_non_zero_index));
xlabel('Time (s)');
ylabel('qdot (rad/s)');
title('Joint Velocities vs. Time');
legend('Joint 1', 'Joint 2', 'Joint 3', 'Joint 4');
grid on;

subplot(2,2,3);
plot(timestamps(1, 1:last_non_zero_index), xHist(:, 1:last_non_zero_index));
xlabel('Time (s)');
ylabel('End-Effector Position (m)');
title('End-Effector Position vs. Time');
legend('x', 'y', 'z');
grid on;

subplot(2,2,4);
errHistNorm = vecnorm(errHist, 2, 1);
plot(timestamps(1, 1:last_non_zero_index), errHistNorm(:, 1:last_non_zero_index));
xlabel('Time (s)');
ylabel('Norm of Position Error (m)');
title('Position Error Norm vs. Time');
grid on;

figure();
plot(timestamps(1,1:last_non_zero_index), manipHist(:,1:last_non_zero_index));
xlabel('Time (s)');
ylabel('Translational Manipulability');
title('Translational Manipulability vs. Time');
ylim([0.0078, .0112])
grid on;



if width(timestamps(1,1:last_non_zero_index)) > width(symQdotHist)
    shortest = width(symQdotHist);
else
    shortest = width(timestamps(1,1:last_non_zero_index));
end


allvdot = [qdotHist(:,1:last_non_zero_index); real_q_dot_hist(:,1:last_non_zero_index)];
plot1 = [allvdot(1,1:shortest); allvdot(5,1:shortest); symQdotHist(1,1:shortest)];
plot2 = [allvdot(2,1:shortest); allvdot(6,1:shortest); symQdotHist(2,1:shortest)];
plot3 = [allvdot(3,1:shortest); allvdot(7,1:shortest); symQdotHist(3,1:shortest)];
plot4 = [allvdot(4,1:shortest); allvdot(8,1:shortest); symQdotHist(4,1:shortest)];

figure();
subplot(2,2,1);
plot(timestamps(1,1:shortest), plot1);
xlabel('Time (s)');
ylabel('Qdot Angular Velocity');
title('Joint 1 Real and Controller Qdot vs. Time');
legend('controller', 'real', 'sim');
grid on;

subplot(2,2,2);
plot(timestamps(1,1:shortest), plot2);
xlabel('Time (s)');
ylabel('Qdot Angular Velocity');
title('Joint 2 Real and Controller Qdot vs. Time');
legend('controller', 'real', 'sim');
grid on;

subplot(2,2,3);
plot(timestamps(1,1:shortest), plot3);
xlabel('Time (s)');
ylabel('Qdot Angular Velocity');
title('Joint 3 Real and Controller Qdot vs. Time');
legend('controller', 'real', 'sim');
grid on;

subplot(2,2,4);
plot(timestamps(1,1:shortest), plot4);
xlabel('Time (s)');
ylabel('Qdot Angular Velocity');
title('Joint 4 Real and Controller Qdot vs. Time');
legend('controller', 'real', 'sim');
grid on;

save(filename, "qHist", "qdotHist", "timestamps")

function [qdot] = ControLaw(desired_pos,desired_velocity,curr_pos,q, params, Kp, Kbar)
    J = openManipulatorJacobian(q, params);
    e = desired_pos(:) - curr_pos(:) ;
    primary = pinv(J)*(desired_velocity(:) + Kp*e);

    gradw = openManipulatorwGradient(q, params);
    P = (eye(4) - pinv(J) * J);

    secondary = P*Kbar*gradw;

    qdot = primary + secondary;

end