%
% FILENAME: binaryStarSystemSimulation.m
%
% DESCRIPTION: Generates Simulations of a binary star system by numerically solving twelve first
%               order equations of motion for the two stars.
%
% AUTHORS:
%  - Graham Driver (https://github.com/GrahamDrive)
%  - Assistance provided by OpenAIs ChatGPT
% 
% CREATED ON: Nov 12, 2023
% 
%

% Masse of Stars
M1 = 2;
M2 = 6;

tspan = [0 10000];
% Initial Position Vector of star 1 & 2
r01 = [-1 2 -1];
r02 = [-1 3 1];

% Velocity of each star velocity two based on conservation of momentum
v01 = [1.2 -0.6 1];
v02 = - M1/M2 .* v01;

% Calculate center of mass and make origin
distanceCOM = (M1*r01 + M2*r02)/(M1+M2);
for i = 1:length(distanceCOM)
    r01(i) = r01(i) - distanceCOM(i);
    r02(i) = r02(i) - distanceCOM(i);
end

% Using ode45 Integrate to find velocity and position of stars
y0 = [r01'; r02'; v01'; v02';];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
starDifferentialsWithMass = @(t, y) starDifferentials(t, y, M1, M2);
[t, y] = ode45(starDifferentialsWithMass, tspan, y0, options);

r1 = y(:,1:3,:);
r2 = y(:,4:6,:);
v1 = y(:,7:9,:);
v2 = y(:,10:12,:);

% Calculate Velocity Magnitude and distance between stars at all points
distance = zeros(length(r1):1);
v1Mag = zeros(length(r1):1);
v2Mag = zeros(length(r1):1);
for i = 1:length(r1)
    v1Mag(i) = norm(v1(i,:));
    v2Mag(i) = norm(v2(i,:));
    distance(i) = norm(r2(i,:) - r1(i,:));
end
distance = distance';
v1Mag = v1Mag';
v2Mag = v2Mag';

% Calculate System Energy
E = (((0.5*M1).*(v1Mag.^2)) + ((0.5*M2).*(v2Mag.^2))) - ((M1*M2)./distance);

% Largest difference in energy
deltaE = max(E) - min(E);

%Calculate Center of mass
distanceCOM = (M1*r1 + M2*r2)/(M1+M2);

%Make Directory for Graphs
path = "Simulations\"+"M1=" + M1 + " M2=" + M2 + " v01 = [" + v01(1,1) + " " + v01(1,2) + " " + v01(1,3) + "]";
mkdir(path);

% Plot Some Graphs

% Adjust this parameter to adjust the range of the 2-D graphs for better
% results. Set to length(t) for all points in graph.
endtime = 10000;

% Controls the speed of the gif by dividing up the total frames. 
% Higher = faster, Lower = slower
gifSpeed = 10000;

% 3-D graphs 

% Binary Seperationg Distance
plot(t(1:endtime), distance(1:endtime,1), "b");
title("Seperation vs Time @ M1 = " + M1 + " M2 = " + M2 + " ΔE: " + deltaE);
xlabel("time")
ylabel("Seperation Distance")
saveas(gcf, path+"\Seperation M1=" + M1 + " M2=" + M2 + ".png");
clf;

% Star 1 Veolcity Graph
plot(t(1:endtime), v1Mag(1:endtime,1), "b");
title("Star 1 Velocity Magnitude vs Time @ M1 = " + M1 + " M2 = " + M2 + " ΔE: " + deltaE);
xlabel("time")
ylabel("Velocity")
saveas(gcf, path+"\Star 1 Velocity Magnitude " + M1 + " M2=" + M2 + ".png");
clf;

% Star 2 Veolcity Graph
plot(t(1:endtime), v2Mag(1:endtime,1), "r");
title("Star 2 Velocity Magnitude @ M1 = " + M1 + " M2 = " + M2 + " ΔE: " + deltaE);
xlabel("time")
ylabel("Velocity")
saveas(gcf, path+"\Star 2 Velocity Magnitude " + M1 + " M2=" + M2 + ".png");
clf;

% Vx Graph
plot(t(1:endtime), v1(1:endtime,1), "b", t(1:endtime), v2(1:endtime,1), "r");
legend(["M1" "M2"]);
title("X-Velocity vs Time @ M1 = " + M1 + " M2 = " + M2 + " ΔE: " + deltaE);
xlabel("time")
ylabel("velocity")
saveas(gcf, path+"\X-Velocity Graph for M1=" + M1 + " M2=" + M2 + ".png");
clf;

% Vy Graph
plot(t(1:endtime), v1(1:endtime,2), "b", t(1:endtime), v2(1:endtime,2), "r");
legend(["M1" "M2"]);
title("Y-Velocity vs Time @ M1 = " + M1 + " M2 = " + M2 + " ΔE: " + deltaE);
xlabel("time")
ylabel("velocity")
saveas(gcf, path+"\Y-Velocity Graph for M1=" + M1 + " M2=" + M2 + ".png");
clf;

% Vz Graph
plot(t(1:endtime), v1(1:endtime,3), "b", t(1:endtime), v2(1:endtime,3), "r");
legend(["M1" "M2"]);
title("Z-Velocity vs Time @ M1 = " + M1 + " M2 = " + M2 + " ΔE: " + deltaE);
xlabel("time")
ylabel("velocity")
saveas(gcf, path+"\Z-Velocity Graph for M1=" + M1 + " M2=" + M2 + ".png");
clf;

% Binary System Orbit Graph without start points marked
plot3(r1(:,1), r1(:,2), r1(:,3), 'b', r2(:,1), r2(:,2), r2(:,3), 'r', ...
    distanceCOM(:,1), distanceCOM(:,2), distanceCOM(:,3), 'g');
hold on
scatter3(r01(1), r01(:,2), r01(3), 100, 'filled', 'b');
scatter3(r02(1), r02(:,2), r02(3), 100, 'filled', 'r');
scatter3(distanceCOM(1,1), distanceCOM(1,2), distanceCOM(1,3), 100, 'filled', 'g');
hold off
title("3-D Plot of Binary @ M1 = " + M1 + " M2 = " + M2 + " ΔE: " + deltaE);
legend(["M1" "M2" "COM"]);
xlabel("X Coordinate");
ylabel("Y Coordinate");
zlabel("Z Coordinate");
saveas(gcf, path+"\StarOrbit Graph with start for M1=" + M1 + " M2=" + M2 + ".png");
clf;

% Binary System Orbit Graph with start points marked
plot3(r1(:,1), r1(:,2), r1(:,3), 'b', r2(:,1), r2(:,2), r2(:,3), 'r', ...
    distanceCOM(:,1), distanceCOM(:,2), distanceCOM(:,3), 'g');
hold on
scatter3(distanceCOM(1,1), distanceCOM(1,2), distanceCOM(1,3), 10, 'filled', 'g');
hold off
title("3-D Plot of Binary Orbit @ M1 = " + M1 + " M2 = " + M2 + " ΔE: " + deltaE);
legend(["M1" "M2" "COM"]);
xlabel("X Coordinate");
ylabel("Y Coordinate");
zlabel("Z Coordinate");
saveas(gcf, path+"\StarOrbit Graph for M1=" + M1 + " M2=" + M2 + ".png");
clf;

% Note this loop will run for a very long time just stop with matlab
i = 1;
while i < 10000
    start = 1;

    plot3(r1(start:i,1), r1(start:i,2), r1(start:i,3), 'b', ...
        r2(start:i,1), r2(start:i,2), r2(start:i,3), 'r', ...
        distanceCOM(start:i,1), distanceCOM(start:i,2), distanceCOM(start:i,3), 'g');
    hold on
    scatter3(r1(i,1), r1(i,2), r1(i,3), 100, 'filled', 'b');
    scatter3(r2(i,1), r2(i,2), r2(i,3), 100, 'filled', 'r');
    scatter3(distanceCOM(1,1), distanceCOM(1,2), distanceCOM(1,3), 100, 'filled', 'g');
    hold off
    title("3-D Plot of Binary Orbit @ M1 = " + M1 + " M2 = " + M2 + " ΔE: " + deltaE);
    legend(["M1" "M2" "COM"]);
    xlabel("X Coordinate");
    ylabel("Y Coordinate");
    zlabel("Z Coordinate");
    exportgraphics(gca, path+"\StarOrbit Graph for M1=" + M1 + " M2=" + M2 + ".gif","Append",true);
    i = round(i + length(r1)/gifSpeed);
end

% Differential Function of the accelration of the stars
function dydt = starDifferentials(t, y, M1, M2)
    r1 = y(1:3);
    r2 = y(4:6);
    v1 = y(7:9);
    v2 = y(10:12);

    [F1, F2] = CalculateForce(M1, M2, r1, r2);
    
    % Acceleration
    a1 = F1./M1;
    a2 = F2./M2;

    dydt = [v1; v2; a1; a2;];
end

%  Force Calculation Function
function [F1, F2] = CalculateForce(M1, M2, r1, r2)
    rMagnitude = norm(abs(r2 - r1)).^3;
    F1 = ((M1*M2)/rMagnitude)*(r2 - r1);

    F2 = -F1;
end