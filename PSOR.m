function [ output_phi ] = PSOR ( dx, dy, w, phi )
% Point Successive Over Relaxation (PSOR) method 
% This method is non-updating
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
% 12 - 5 - 2016
%% inputs
% dx      : space iteration interval in x-dirrtion
% dy      : space iteration interval in y-dirrtion
% w       : relaxation factor 0 < w < 2
% phi    : [phi_(ij), phi_(i+1j), phi_(i-1j), phi_(ij+1), phi_(ij-1)] at time n
%% outputs
% output_phi   : phi_(ij) at time n+1
%% function body
r=dx./dy;
phi_telta_ij=0.5*((1+r.^2).^(-1)).*(phi{2}+phi{3}+r.^2.*(phi{4}+phi{5}));
output_phi=phi{1}+w*(phi_telta_ij-phi{1});
end

