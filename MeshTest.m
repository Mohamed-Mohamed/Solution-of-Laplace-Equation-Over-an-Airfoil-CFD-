% This code is designed to test your mesh before solving as your grid must intersect the airfoil at the leading edge and trailing edge
% Fill the required data  in DATA section
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
% 12 - 5 - 2016
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all ; clc;
%% DATA
    %% loading data
    % import airfoil coordinates
    % the airfoil must be a closed countor
    data=importdata('airfoil.txt');
    for i=1:length(data(:,1))
        x1(i,1)=data(i,1);
        y1(i,1)=data(i,2);
    end
    % angle of attack in degree
    Alpha=10;
    % free stream velocity
    v_inf=100;
    % airfoil chord
    c=1;
    %% additional data
        %% grid data
        % length of domain in x-direction
        Lx=2; 
        % length of domain in y-direction
        Ly=2;          
        % interval in x-direction
        dx=0.01;       
        % interval in y-direction
        dy=0.01;       
%% meshing
[ nx, ny, XG, YG, dX, dY, x_body_index, y_body_index ] = Hgrid ( x1, y1, c, Lx, Ly, dx, dy, [1 1] );