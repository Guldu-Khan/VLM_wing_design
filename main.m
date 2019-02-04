close all; clear all; clc;
%% Preliminaries

% This is the main function based on the vortex lattice method 
% Program is based on rectangular vortex loops to calculate velocities
% Refer Katz for coordinate frames and any other general info 

%% Input geometry parameters

% Currently written for rectangular wing

global b              % Semi-span
global c              % Chord
global n              % Number of elements in the y - direction
global m              % Number of elements in the x - direction

% Used for computing length of vortex rings
global rec_x
global rec_y

b = 20;
c = 1;
n = 20;
m = 2;

rec_y = 2*b/n;
rec_x = c/m;

lift_t = [];


for alpha = -0.5:0.1:0.5


Qinf = [1; 0; alpha];
Rho = 1;

%% Constructing the matrix to compute gamma

A = zeros(n*m);

for r=1:n*m
    for s=1:n*m

        % Coordinates of the vortex element
        i = mod(s,n);
        if i==0
            i=n;
        end
        j = ceil(s/n);      
        
        % Coordinates of the collocation point
        I = mod(r,n);
        if I==0
            I=n;
        end
        J = ceil(r/n);      
        
        [f,g] = meshcor(i,j);
        [F,G] = meshcor(I,J);
        
        A(r,s) = dot(vring(F,i,j),G);
        
    end
end

%% Constructing the RHS vector 

% This contains information about the boundary

for k = 1:n*m
    
   i = mod(k,n);
   if i==0
            i=n;
   end
   j = ceil(k/n);
   
   [f,g] = meshcor(i,j);
   Rhs(k) = dot(Qinf,g);

end

%% Invert matrices to get gamma

Gamma = (A\Rhs')';

%% Compute lift and drag

% Reconstructing gamma back into the real geometry from a vector

for k = 1:n*m
    
   i = mod(k,n);
   if i==0
            i=n;
   end
   j = ceil(k/n);

   gamma(i,j) = Gamma(k);

end

% Computing lift at each element
Lift = 0;
for i = 1:n
    for j = 1:m
        
        if i>1
            
            l(i,j) = Rho*norm(Qinf)*(gamma(i,j) - gamma(i-1,j))*rec_y;
        
        else
            
            l(i,j) = Rho*norm(Qinf)*gamma(i,j)*rec_y;
            
        end
        
        % Pressure distribution
        p(i,j) = l(i,j)/(rec_x*rec_y);
        
        % Total lift
        Lift = Lift + l(i,j);
    end
end

lift_t = [lift_t Lift];


end

%%
Cl = lift_t / (0.5*Rho*norm(Qinf)^2*c*2*b);

Alpha = linspace(-0.5,0.5,11);
plot(Alpha,Cl);