function [ nx, ny, XG, YG, dX, dY, x_body_index, y_body_index  ] = Hgrid ( x, y, c, Lx, Ly, dx, dy, Plot )
% this function is used to generate H grid
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
% 12 - 5 - 2016
%% inputs
% x         : x co-ordinate of airfoil (its length must be even)
% y         : y co-ordinate of airfoil (its length must be even)
% c         : is the chord of airfoil
% Lx       : length of domain in x-direction
% Ly       : length of domain in y-direction
% dx       : interval in x-direction
% dy       : interval in y-direction
% Plot    : if Plot(0)  -->  no plot
%             if Plot(1)  -->   plot  and Plot(2) is figure number
%% outputs
% nx                        : number of grid point in x-direction
% ny                        : number of grid point in y-direction
% XG                       : the x matrix of the grid
% YG                       : the y matrix of the grid
% dX                       : the interval in x-direction as a matrix
% dY                       : the interval in y-direction as a matrix
% x_body_index     : the index of the two rows in YG matrix which contain the upper and lower surfaces of the airfoil (it's not correct, but it's correct in the main code)
% y_body_index     : the index of the columns in YG which the vertical lines of grid intersect the airfoil
%% function body
nx=floor(Lx/dx);
ny=floor(Ly/dy/2);
x_min=min(x);
x_max=max(x);
y_min=min(y);
y_max=max(y);
for n=1:nx
    xg(n)=(n-1)*dx;
end
for m=1:ny
    yg(m)=(m-1)*dy;
end
[XG,YG_upper]=meshgrid(xg,yg);
YG_lower=-YG_upper;
x_b=floor(length(xg)/2);
x_body=x+xg(x_b+1)-c/2;
y_b=floor(length(yg)/2);
y_body=y+yg(y_b*2)+dy;
NN=zeros(1,length(xg));
for xx=1:length(xg)
    [x1, y1]=curveintersect(x_body(1:length(x_body)/2), y(1:length(y)/2), xg(xx)*ones(1,2*length(yg)), [-yg(end:-1:1), yg]);
    if length(x1)==0
        x1=0;
        NN(xx)=1;
    end
    if length(y1)==0
        y1=0;
    end
    x_in_upper(xx)=x1;
    y_in_upper(xx)=y1;
end
dY_upper=(YG_upper(end,:)-y_in_upper)/(ny-1);
for k=1:ny
    YG_upper(k,:)=y_in_upper+(k-1)*dY_upper;
end

for yy=1:length(xg)
    [x1, y1]=curveintersect(x_body(length(x_body)/2:end), y(length(y)/2:end), xg(yy)*ones(1,2*length(yg)), [-yg(end:-1:1), yg]);
    if length(x1)==0
        x1=0;
    end
    if length(y1)==0
        y1=0;
    end
    x_in_lower(yy)=x1;
    y_in_lower(yy)=y1;
end
dY_lower=(YG_lower(end,:)-y_in_lower)/(ny-1);
for k=1:ny
    YG_lower(k,:)=y_in_lower+(k-1)*dY_lower;
end
LL=1;
for XxX=length(xg):-1:1
    if NN(XxX) == 1
        NN(XxX)=[];
    else
        y_body_index(LL)=XxX;
        LL=LL+1;
    end
end
y_body_index=sort(y_body_index);
x_body_index=[ny/2+1;ny/2+2];
YG_upper(1:end,:);
YG = [ans(end:-1:1,:)+Ly/2; YG_lower(1:end,:)+Ly/2];
XG = [XG(1:end,:); XG(1:end,:)];
DY_upper=abs(YG(1,:)-YG(2,:));
DY_lower=abs(YG(end-1,:)-YG(end,:));
for l=1:ny
    dY(l,:)=DY_upper;
end
for l=ny+1:2*ny
    dY(l,:)=DY_lower;
end
[ny, nx]=size(YG);
dX=ones(size(XG));
%% plotting
if Plot(1)==1
    figure(Plot(2));
    view(0,90);
    set(gcf,'Color','w')
    hold all
    mesh(XG,YG,XG*0-0.01,'edgecolor','red');
    area(x_body,y_body);
%     xlim([-c/2,Lx+c/2]);
%     ylim([-c/2,Ly+c/2]);
    xlabel('X','Fontsize',18)
    ylabel('Y','Fontsize',18)
    legend('H-Grid', 'Airfoil')
end
end