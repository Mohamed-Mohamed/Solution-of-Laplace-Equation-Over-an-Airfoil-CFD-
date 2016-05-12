% This code is designed to test your airfoil is a closed countor or not
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
%% Airfoil plotting
figure(1);
set(gcf,'Color','w');
plot(x1, y1, 'linewidth', 2);
grid on;
xlabel('X-axis','fontsize',18);
ylabel('Y-axis','fontsize',18);
legend('Airfoil');
axis equal    
    
    
    
    
    
    
    