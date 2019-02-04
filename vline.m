%% Function to compute velocities due to a vortex line of finite length 

% Refer Katz Chapter 10.4 page 253

function A = vline(P,L1,L2)


% P = [0; 1; 0];            % Position of the point where velocity needs to be computed
% L1 = [1; 0; 0];         % Position of one end of the vortex line
% L2 = [-1; 0; 0];         % Position of the other end of the vortex line

e =1e-6;                % Setting threshold to terminate program

%% Various vector parameters

R1 = L1 - P;
R2 = L2 - P;
R0 = L2 -L1;
R = cross(R1,R2);

%% Checking for singularities 

if norm(R1)<e || norm(R2)<e || norm(R)<e
    
    A=[0; 0; 0];            % Velocity induced is set to zero

else
    
%% Computing final velocity induced  

D1 = dot(R0,R1);
D2 = dot(R0,R2);

gamma = 1;

K = (gamma/(4*pi))*(D1/norm(R1) - D2/norm(R2))/norm(R)^2 ;

A = K*R;

end

