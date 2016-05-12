% This code is designed to solve Laplace Equation "nabla^2 Psi = 0" using Point Successive Over Relaxation (PSOR) method 
% Fill the required data  in DATA section
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
% 12 - 5 - 2016
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all ; clc;
tic;
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
        Lx=4; 
        % length of domain in y-direction
        Ly=2;
        % interval in x-direction
        dx=0.1;       
        % interval in y-direction
        dy=0.1;       
        %% PSOR data
        % relaxation factor  0 < w < 2
        w=0.5;  
        %% plotting data
        % number of countors
        no_of_countors=400;
        %% error limits data
        % order of error of stream line function
        error_order_psi=1e-10; 
        % order of error of stream line of body (airfoil)
        error_order_psi_body=1e-6; 
        %% animation data
        % NOTE THAT IT WILL TAKE A HIGH MEMORY SPACE (BE CARFULL)
        % if animation=0 -> no animation needed
        % if animation=1 -> animation will work (animation of solution of stream lines, velocity and pressure coefficient at every iteration)
        animation=0;
        % frame rate
        % number of frames = length(1:frame_rate:no_of_iteration)
        frame_rate=100;
        % to save your machine from infinite loop
        no_of_iteration=10000; 
%% meshing
[ nx, ny, XG, YG, dX, dY, x_body_index, y_body_index ] = Hgrid ( x1, y1, c, Lx, Ly, dx, dy, [1 1] );
%% Stream vel components
V_x=v_inf*cosd(Alpha);
V_y=v_inf*sind(Alpha);
%% boundary condtions and initialization
% initialization
psi=zeros(ny,nx);
psi_old=zeros(ny,nx);
vel_x=zeros(ny,nx);
vel_x_old=zeros(ny,nx);
vel_y=zeros(ny,nx);
vel_y_old=zeros(ny,nx);
vel=zeros(ny,nx);
vel_old=zeros(ny,nx);
% psi BC's
% left boundary
for leb=ny-1:-1:1
    psi(leb,1)=psi(leb+1,1)+V_x*dy;
end
% upper boundary
for upb=1:nx-1
    psi(1,upb+1)=psi(1,upb)-V_y*dx;
end
% right boundary
for rib=1:ny-1
    psi(rib+1,end)=psi(rib,end)-V_x*dy;
end
% lower boundary
for lob=nx:-1:2
    psi(end,lob-1)=psi(end,lob)+V_y*dx;
end
% velocity BC's
vel_x(1:end,1)=V_x;
vel_x(1,1:end)=V_x;
vel_x(end,1:end)=V_x;
vel_x(1:end,end)=V_x;
vel_y(1:end,1)=V_y;
vel_y(1,1:end)=V_y;
vel_y(end,1:end)=V_y;
vel_y(1:end,end)=V_y;
%% helpin parameter
psi_body(1)=5;
psi_body_error(1)=0.1;
psi_error(1)=0.1;
psi(2*(x_body_index(1)-1)+1,y_body_index(end)+1)=1;
t=1;
T=1;
%% Psi solution
while (psi_body_error(t-1+T) >= error_order_psi_body || psi_error(t-1+T) >= error_order_psi) && t <= no_of_iteration
    %% Psi solution
        if t == 1 && animation==0
            for liner_upper=2:nx-1
                psi(:,liner_upper)=linspace(psi(1,liner_upper),psi(end,liner_upper),ny);
            end
        else
            [ output_psi ] = PSOR ( dX(2:end-1,2:end-1), dY(2:end-1,2:end-1), w, {psi_old(2:end-1,2:end-1), psi_old(2:end-1,3:end), psi_old(2:end-1,1:end-2), psi_old(1:end-2,2:end-1), psi_old(3:end,2:end-1)} );
            psi(2:end-1,2:end-1)=output_psi;
        end
        % psi body
%         psi_body(t+1)=(psi(2*(x_body_index(1)-1)+1,y_body_index(end)+1));
        psi_body(t+1)=(psi(2*(x_body_index(1)-1)+1,y_body_index(end)+1)+psi(2*(x_body_index(1)-1),y_body_index(end)+1))/2;
        psi(2*(x_body_index(1)-1):2*(x_body_index(1)-1)+1,y_body_index(1):y_body_index(end))=psi_body(t+1);
        psi_error(t)=max(max(abs(psi_old-psi)));
        psi_old=psi;
        % psi body error
        psi_body_error(t)=abs(psi_body(t+1)-psi_body(t));
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if animation ==1
         % psi
         PSI{t}=psi;
         % velocity
         vel_x(2:end-1,2:end-1)=psi(1:end-2,2:nx-1)-psi(3:end,2:nx-1);
         vel_y(2:end-1,2:end-1)=-(psi(2:end-1,3:end)-psi(2:end-1,1:end-2));
         vel_x(2:end-1,2:end-1)=vel_x(2:end-1,2:end-1)./dY(2:end-1,2:end-1)/2;
         vel_y(2:end-1,2:end-1)=vel_y(2:end-1,2:end-1)./dX(2:end-1,2:end-1)/2;
         vel=sqrt(vel_x.^2+vel_y.^2);
         VEL_X{t}=vel_x;
         VEL_Y{t}=vel_y;
         VEL{t}=vel;
         % Cp
         Cp=1-(vel/v_inf).^2;
         CP{t}=Cp;
     end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if animation == 0
         no_of_iteration=no_of_iteration*2*t;
     end
     T=0;
     t=t+1;
end
%% velocity solution
vel_x(2:end-1,2:end-1)=psi(1:end-2,2:nx-1)-psi(3:end,2:nx-1);
vel_y(2:end-1,2:end-1)=-(psi(2:end-1,3:end)-psi(2:end-1,1:end-2));
vel_x(2:end-1,2:end-1)=vel_x(2:end-1,2:end-1)./dY(2:end-1,2:end-1)/2;
vel_y(2:end-1,2:end-1)=vel_y(2:end-1,2:end-1)./dX(2:end-1,2:end-1)/2;
vel=sqrt(vel_x.^2+vel_y.^2);
%% Cp solution
Cp=1-(vel/v_inf).^2;
%% Time
toc
%% PLOTTING
if animation == 0
    %% NO ANIMATION
        %% psi plotting
        figure(2);
        set(gcf,'Color','w')
        surf(XG, YG, psi/c/v_inf,'EdgeColor','none')
        hold all;
        plot3(x1+XG(1,y_body_index(1)),y1+Ly/2,max(max(abs(psi/c/v_inf)))*ones(1,length(x1)),'color',[0,0,0],'linewidth',2);
        view(0,90)
        colorbar;
        xlabel('X-axis','fontsize',18);
        ylabel('Y-axis','fontsize',18);
        zlabel('\psi/V_i_n_f/c','fontsize',18);
        title('Non-dimentional Stream function','fontsize',18);
        %% vel plotting
        figure(3);
        set(gcf,'Color','w')
        surf(XG, YG, vel/v_inf,'EdgeColor','none')
        hold all;
        plot3(x1+XG(1,y_body_index(1)),y1+Ly/2,max(max(abs(vel/v_inf)))*ones(1,length(x1)),'color',[0,0,0],'linewidth',2);
        view(0,90)
        h=colorbar;
        xlabel('X-axis','fontsize',18);
        ylabel('Y-axis','fontsize',18);
        zlabel('V/V_i_n_f','fontsize',18);
        title('Non-dimentional vel','fontsize',18);
        %% Cp plotting
        figure(4);
        set(gcf,'Color','w')
        surf(XG, YG, Cp,'EdgeColor','none')
        hold all;
        plot3(x1+XG(1,y_body_index(1)),y1+Ly/2,max(max(abs(Cp)))*ones(1,length(x1)),'color',[0,0,0],'linewidth',2);
        view(0,90)
        colorbar;
        xlabel('X-axis','fontsize',18);
        ylabel('Y-axis','fontsize',18);
        zlabel('C_P','fontsize',18);
        title('C_p Ditrbution','fontsize',18);
        %% psi body ploting
        figure(5);
        set(gcf,'Color','w')
        plot((0:length(psi_body)-1),psi_body,'linewidth',2);
        xlabel('X-axis','fontsize',18);
        ylabel('\Psi_b_o_d_y','fontsize',18);
        title('\Psi_b_o_d_y','fontsize',18);
        grid on;
        %% psi body error ploting 
        figure(6);
        set(gcf,'Color','w')
        plot((0:length(psi_body_error)-1),log10(psi_body_error),'linewidth',2);
        xlabel('X-axis','fontsize',18);
        ylabel('log_1_0(Max. |error|)','fontsize',18);
        title('Max. |error| of \Psi_b_o_d_y','fontsize',18);
        grid on;
        %% psi ploting error
        figure(7);
        set(gcf,'Color','w')
        plot((0:length(psi_error)-1),log10(psi_error),'linewidth',2);
        xlabel('X-axis','fontsize',18);
        ylabel('log_1_0(Max. |error|)','fontsize',18);
        title('Max. |error| of \Psi','fontsize',18);
        grid on;
        %% velocity over airfoil at the end of solution
        % vel_upper_airfoil
        vel_upper_airfoil=vel(2*(x_body_index(1)-1)-1,y_body_index(1)-1:y_body_index(end)+1)/v_inf;
        % vel_lower_airfoil
        vel_lower_airfoil=vel(2*(x_body_index(1)-1)+2,y_body_index(1)-1:y_body_index(end)+1)/v_inf;
        figure(8);
        set(gcf,'Color','w');
        hold all;
        hx=linspace(0,c,length(vel_upper_airfoil));
        plot((0:length(vel_upper_airfoil)-1)*(hx(2)-hx(1)),vel_upper_airfoil,'linewidth',2);
        plot((0:length(vel_lower_airfoil)-1)*(hx(2)-hx(1)),vel_lower_airfoil,'linewidth',2);
        plot(x1,y1+1)
        xlabel('Chord line','fontsize',18);
        ylabel('Non-dimentional Velocity','fontsize',18);
        legend('Upper surface','Lower surface','Airfoil');
        title('Velocity Distrbution Over Airfoil','fontsize',18)
        grid on;
        %% Cp over airfoil at the end of solution
        Cp_upper_airfoil=1-vel_upper_airfoil.^2;
        Cp_lower_airfoil=1-vel_lower_airfoil.^2;
        figure(9);
        set(gcf,'Color','w');
        hold all;
        plot((0:length(Cp_upper_airfoil)-1)*(hx(2)-hx(1)),Cp_upper_airfoil,'linewidth',2);
        plot((0:length(Cp_lower_airfoil)-1)*(hx(2)-hx(1)),Cp_lower_airfoil,'linewidth',2);
        plot(x1,y1)
        xlabel('Chord line','fontsize',18);
        ylabel('Cp','fontsize',18);
        title('C_p Distrbution Over Airfoil','fontsize',18)
        legend('Upper surface','Lower surface','Airfoil');
        grid on;
        %% Stream Lines Countors
        figure(10);
        hold all;
        set(gcf,'Color','w');
        contour(XG, YG, psi,linspace(min(min(psi)),max(max(psi)),no_of_countors));
        plot(x1+XG(1,y_body_index(1)),y1+Ly/2,'color',[0,0,0],'linewidth',2);
        colorbar;
        xlabel('X-axis','fontsize',18);
        ylabel('Y-axis','fontsize',18);
        title('Stream Lines Countors','fontsize',18)
        %% pressure Lines Countors
        figure(11);
        hold all;
        set(gcf,'Color','w');
        contour(XG, YG, Cp,linspace(min(min(Cp)),max(max(Cp)),no_of_countors))
        plot(x1+XG(1,y_body_index(1)),y1+Ly/2,'color',[0,0,0],'linewidth',2);
        colorbar;
        xlabel('X-axis','fontsize',18);
        ylabel('Y-axis','fontsize',18);
        title('Cp Countors','fontsize',18)
        %% velocity vector
        figure(12);
        hold all;
        quiver(XG, YG, vel_x, vel_y)
        plot(x1+XG(1,y_body_index(1)),y1+Ly/2,'color',[0,0,0],'linewidth',2);
        xlabel('X-axis','fontsize',18);
        ylabel('Y-axis','fontsize',18);
        title('Velocity Vector','fontsize',18)
        xlim([0, Lx]);
        ylim([0, Ly]);
        set(gcf,'Color','w');
elseif animation == 1    
    %% ANIMATION
        %% PSI plotting
        figure(2);
        for P=1:frame_rate:length(PSI)
            clf;
            set(gcf,'Color','w')
            surf(XG, YG, PSI{P}/c/v_inf,'EdgeColor','none')
            hold all;
            plot3(x1+XG(1,y_body_index(1)),y1+Ly/2,max(max(abs(PSI{end}/c/v_inf)))*ones(1,length(x1)),'color',[0,0,0],'linewidth',2);
            view(0,90)
            colorbar;
            xlabel('X-axis','fontsize',18);
            ylabel('Y-axis','fontsize',18);
            zlabel('\psi/V_i_n_f/c','fontsize',18);
            title('Non-dimentional Stream function','fontsize',18);
            pause(1/100);
        end
        %% VEL plotting
        figure(3);
        for V=1:frame_rate:length(VEL)
            clf;
            set(gcf,'Color','w')
            surf(XG, YG, VEL{V}/v_inf,'EdgeColor','none')
            hold all;
            plot3(x1+XG(1,y_body_index(1)),y1+Ly/2,max(max(abs(VEL{end}/v_inf)))*ones(1,length(x1)),'color',[0,0,0],'linewidth',2);
            view(0,90)
            colorbar;
            xlabel('X-axis','fontsize',18);
            ylabel('Y-axis','fontsize',18);
            zlabel('V/V_i_n_f','fontsize',18);
            title('Non-dimentional vel','fontsize',18);
            pause(1/100);
        end
        %% CP plotting
        figure(4);
        for C=1:frame_rate:length(CP)
            clf;
            set(gcf,'Color','w')
            surf(XG, YG, CP{C},'EdgeColor','none')
            hold all;
            plot3(x1+XG(1,y_body_index(1)),y1+Ly/2,max(max(abs(CP{end})))*ones(1,length(x1)),'color',[0,0,0],'linewidth',2);
            view(0,90)
            colorbar;
            xlabel('X-axis','fontsize',18);
            ylabel('Y-axis','fontsize',18);
            zlabel('C_P','fontsize',18);
            title('C_p Ditrbution','fontsize',18);
            pause(1/100);
        end
        %% psi body ploting
        figure(5);
        set(gcf,'Color','w')
        plot((0:length(psi_body)-1),psi_body,'linewidth',2);
        xlabel('X-axis','fontsize',18);
        ylabel('\Psi_b_o_d_y','fontsize',18);
        title('\Psi_b_o_d_y','fontsize',18);
        grid on;
        %% psi body error ploting 
        figure(6);
        set(gcf,'Color','w')
        plot((0:length(psi_body_error)-1),log10(psi_body_error),'linewidth',2);
        xlabel('X-axis','fontsize',18);
        ylabel('log_1_0(Max. |error|)','fontsize',18);
        title('Max. |error| of \Psi_b_o_d_y','fontsize',18);
        grid on;
        %% psi ploting error
        figure(7);
        set(gcf,'Color','w')
        plot((0:length(psi_error)-1),log10(psi_error),'linewidth',2);
        xlabel('X-axis','fontsize',18);
        ylabel('log_1_0(Max. |error|)','fontsize',18);
        title('Max. |error| of \Psi','fontsize',18);
        grid on;
        %% velocity over airfoil at the end of solution
        % vel_upper_airfoil
        vel_upper_airfoil=vel(2*(x_body_index(1)-1)-1,y_body_index(1)-1:y_body_index(end)+1)/v_inf;
        % vel_lower_airfoil
        vel_lower_airfoil=vel(2*(x_body_index(1)-1)+2,y_body_index(1)-1:y_body_index(end)+1)/v_inf;
        figure(8);
        set(gcf,'Color','w');
        hold all;
        hx=linspace(0,c,length(vel_upper_airfoil));
        plot((0:length(vel_upper_airfoil)-1)*(hx(2)-hx(1)),vel_upper_airfoil,'linewidth',2);
        plot((0:length(vel_lower_airfoil)-1)*(hx(2)-hx(1)),vel_lower_airfoil,'linewidth',2);
        plot(x1,y1+1,'linewidth',2);
        xlabel('Chord line','fontsize',18);
        ylabel('Non-dimentional Velocity','fontsize',18);
        legend('Upper surface','Lower surface','Airfoil');
        title('Velocity Distrbution Over Airfoil','fontsize',18)
        grid on;
        %% Cp over airfoil at the end of solution
        Cp_upper_airfoil=1-vel_upper_airfoil.^2;
        Cp_lower_airfoil=1-vel_lower_airfoil.^2;
        figure(9);
        set(gcf,'Color','w');
        hold all;
        plot((0:length(Cp_upper_airfoil)-1)*(hx(2)-hx(1)),Cp_upper_airfoil,'linewidth',2);
        plot((0:length(Cp_lower_airfoil)-1)*(hx(2)-hx(1)),Cp_lower_airfoil,'linewidth',2);
        plot(x1,y1,'linewidth',2);
        xlabel('Chord line','fontsize',18);
        ylabel('Cp','fontsize',18);
        title('C_p Distrbution Over Airfoil','fontsize',18)
        legend('Upper surface','Lower surface','Airfoil');
        grid on;
        %% Stream Lines Countors
        figure(10);
        hold all;
        set(gcf,'Color','w');
        contour(XG, YG, PSI{end},linspace(min(min(PSI{end})),max(max(PSI{end})),no_of_countors))
        plot(x1+XG(1,y_body_index(1)),y1+Ly/2,'color',[0,0,0],'linewidth',2);
        colorbar;
        xlabel('X-axis','fontsize',18);
        ylabel('Y-axis','fontsize',18);
        title('Stream Lines Countors','fontsize',18)
        %% pressure Lines Countors
        figure(11);
        hold all;
        set(gcf,'Color','w');
        contour(XG, YG, CP{end},linspace(min(min(CP{end})),max(max(CP{end})),no_of_countors))
        plot(x1+XG(1,y_body_index(1)),y1+Ly/2,'color',[0,0,0],'linewidth',2);
        colorbar;
        xlabel('X-axis','fontsize',18);
        ylabel('Y-axis','fontsize',18);
        title('Cp Countors','fontsize',18)
        %% velocity vector
        figure(12);
        hold all;
        quiver(XG, YG, VEL_X{end}, VEL_Y{end})
        plot(x1+XG(1,y_body_index(1)),y1+Ly/2,'color',[0,0,0],'linewidth',2);
        xlabel('X-axis','fontsize',18);
        ylabel('Y-axis','fontsize',18);
        title('Velocity Vector','fontsize',18)
        xlim([0, Lx]);
        ylim([0, Ly]);
        set(gcf,'Color','w');
end