function A = vring(P,i,j)

global m;
global rec_x;
global rec_y;

% Call vline to construct the loop
% First generate the coordinates of the end points of each vline  based on n and m values 
% Currently this is for flat plate.
% Care needs to be taken for those flapped cells 

L = meshcor(i,j);                       % Gives coordinates of the center

% Reconstructing the square loop from the coordinates of the center point

L1 = L + [-rec_x/2; -rec_y/2; 0]; 
L2 = L + [-rec_x/2; rec_y/2; 0];
L3 = L + [rec_x/2; rec_y/2; 0];
L4 = L + [rec_x/2; -rec_y/2; 0];

% Wake modelling

L1inf = L + [1e6; -rec_y/2; 0];
L2inf = L + [1e6; rec_y/2; 0];

% Differentiating between wake points and interior points 

if j ~= m
    
    A = vline(P,L1,L2) + vline(P,L2,L3) + vline(P,L3,L4) + vline(P,L4,L1);
    
else
    
    A = vline(P,L1,L2) + vline(P,L1inf,L1) + vline(P,L2,L2inf);
    
end



