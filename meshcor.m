function [L,normal] = meshcor(i,j)

% Function must include formula to find out normal vector at each element
% Presently it is assumed to be z - direction

global rec_x;
global rec_y;
global b;

    L = [(j*rec_x - rec_x/2); (i*rec_y - (b+rec_y/2)); 0];
    normal = [0; 0; 1];
    
   
          