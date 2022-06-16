clc
close all
clear all

basename = 'Solver_';
global file_number

for file_number=1
clearvars -except file_number basename intps convfs
rng('shuffle');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Global Parameters

global ncell                %Number of Cells
global npoints              %Number of Cell Points
global c                    %Cell Number
global Cell_Point_Number    %Point on Cell                               

global l0_out               %Initial Length Between Cell Points, um
global l0_in                %Initial Length Between Nucleus Points, um \
global l0_in2
global l0_inout             %Initial Length Between Cell Points and Nucleus Points, um
global l0_nu                %Initial Length Between Internal Nucleus Points, um
global l_bond               %Initial Length Between Cell to Cell Binding, um
global thresh_max           %Max Cell to Cell binding Distance, um
global thresh_min           %Min Cell to Cell binding Distance, um
%global dist_nr

global k_cout               %Stiffness Value of Cell Membrane, nN/um
global k_cin                %Stiffness Value of Nucleus Membrane, nN/um
global k_nu                 %Stiffness Value of Internal Nucleus, nN/um
global k_cell               %Stiffness Value of Cell to Cell Binding, nN/um
global k_in                 %Stiffness Value of Cell Membrane to Nucleus, nN/um

global F_external           %Forces on Cell Points, nN
global F_intx               %Forces on Nucleus Points X-direction, nN
global F_inty               %Forces on Nucleus Points X-direction, nN
global force_number

global v                    %Viscosity Coefficient, nN*s/um
global vin
global Friction 

global total_t              %Total Time of Simulation

global ecm_thresh_max       %Max Cell to ECM binding Distance, um
global ecm_thresh_min       %Max Cell to ECM binding Distance, um
global ecm_point            %ECM Point Number
global ecm_x                %ECM X-Position
global ecm_y                %ECM Y-Position
%global ecm_bond             %Initial Length Between Cell to ECM Binding, um
global Conforce             %Contractile Force
global ttt                  %Counter
global internalP            %internal Pressure 
global thiccness
global area
%global pressure
global time_step
global interval
global interpStruct
global isint
global intbound
global nucpull
global attached
global simulation_length
global nucspring
global isinside
global pushoutforce
global k_cell
global k_cellv
global k_coutbase
global intspring
global mindistmembrane
global outsidedist
global isconforce
global K_outerlim
global circularity
global basedist
global PhaseShift
global ConAmp
global ConPer
global randorder
global outerviscscale
global exConScale
global Friction
global OVL
global isattached
global iswhole


global forceCap
%load wknd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialization for Acini Creation

point = 0;                                                  %Starting Cell Point
point_out = 0;                                              %Starting Nucleus Point

%Acini Structure Geometry
center = [0,0];                                             %Center of Structure
ncell = 8;                                                  %Number of Nucleus on Structure
npoints = 20;                                               %Number of Points per Nucleus must be multiple of 4

cellrad=5;                                                  %Radius of Cells um
nuclrad=2;                                                  %Radius of Nucleus um
r = nuclrad;                                                
r_out = cellrad;              
basedist=(cellrad-nuclrad)*1.1;                             

points=NaN(ncell,npoints,4);
structrad=(2*cellrad)/(2*sin(pi/ncell))*0.85;               %Calculates Structure Radius. Scaled up or down to squish/stretch cells
%centers=pointgen([0,0],ncell,structrad);
axisscaler=3;                                               %How much bigger than the structure is the output video? 
axis(axisscaler*[-structrad structrad -structrad structrad]); 
hold on;
s1=NaN;
outsidedist=structrad+cellrad*2.5;                          %Distance to outside membrane 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%run specific controls   0=false,1=true
extracell=0;                                                %decides if an extra cell will be added during the run
%isconforce=rem(file_number,2);                             %spring=0, contractile=1
isconforce=2;

nucspring=1;                                                %0=no spring, 1=just spring, 2=both. no %%spring does not work well %probably for replacing mem2nuc springs
alwayssame=1;                                               %
istrap=1;                                                   %Is the Cell a trapazoidal shape? 
roundnuc=1;                                                 %is the nucleus a circle? 
parameterset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Original Plot + Defining 'c' and 'Cell_Point_Number'
x0=NaN(1,ncell*npoints*4);
if istrap~=1
    fig1=figure;
    nucrat=nuclrad/cellrad;
    structrad=(2*cellrad)/(2*sin(pi/ncell))*1.0;  %0.65
    
    nppts=circles(ncell,npoints,structrad,cellrad,0,nucrat);
    structrad=structrad+nuclrad*1;
    if roundnuc==1
        nnppts=circles(ncell,npoints,structrad,nuclrad,0,nucrat);
    else
        nnppts=traps(ncell,npoints,structrad,nucrat,1,nucrat);
    end
        


else
    nucrat=nuclrad/cellrad;
    structrad=(2*cellrad)/(2*sin(pi/ncell))*1.30;  %0.65
    
    nppts=traps(ncell,npoints,structrad,cellrad,0,nucrat);
    structrad=structrad+nuclrad*0.5;
    if roundnuc==1
        nnppts=circles(ncell,npoints,structrad,nuclrad,0,nucrat);
    else
        nnppts=traps(ncell,npoints,structrad,nucrat,1,nucrat);
    end
     
end
    x0=NaN(1,ncell*npoints*2);
    %%% assigns output values to positions in vector
    for n=1:ncell
        x0(((n-1)*npoints+1):((n)*npoints))=nppts(((n-1)*npoints+1):(n-1)*npoints+npoints,1);
        x0(ncell*npoints+(n-1)*npoints+1:ncell*npoints +(n)*npoints)=nppts(((n-1)*npoints+1):(n-1)*npoints+npoints,2);
        x0(2*ncell*npoints+npoints*(n-1)+1:2*ncell*npoints+npoints*(n))=nnppts(((n-1)*npoints+1):(n-1)*npoints+npoints,1);
        x0(3*ncell*npoints+npoints*(n-1)+1:3*ncell*npoints+npoints*(n))=nnppts(((n-1)*npoints+1):(n-1)*npoints+npoints,2);

    end       
        pm=nppts;
        pn=nnppts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initial positions of points
for n=1:ncell
    ppnts(:,n,1)=x0(((n-1)*npoints+1):((n)*npoints));
    ppnts(:,n,2)=x0(ncell*npoints+(n-1)*npoints+1:ncell*npoints +(n)*npoints);
    plot(polyshape(ppnts(:,n,1),ppnts(:,n,2)));
    hold on;
    scatter(ppnts(:,n,1),ppnts(:,n,2),'*');
    hold on;
    ppnts(:,n,1)=x0((2*ncell*npoints+((n-1)*npoints+1)):2*ncell*npoints+((n)*npoints));
    ppnts(:,n,2)=x0((3*ncell*npoints+((n-1)*npoints+1)):3*ncell*npoints+((n)*npoints));
    plot(polyshape(ppnts(:,n,1),ppnts(:,n,2)));
    scatter(ppnts(:,n,1),ppnts(:,n,2),'*');
end
if ttt<0    %numbers each point visually. same number between nucleus and membrane means connected
for s=1:ncell*npoints

   text(x0(s),x0(s+ncell*npoints), string(s));
   text(x0(s+ncell*npoints*2),x0(s+ncell*npoints*3), string(s));
   
    
end
end


Cell_Point_Number=zeros(4*ncell*npoints,1);
for n=1:ncell
   Cell_Point_Number(((n-1)*npoints+1):((n)*npoints))=1:npoints;
   c(((n-1)*npoints+1):((n)*npoints))=n;
end
isattached=zeros(4*ncell*npoints,1);

hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Global Parameters That Rely on initial conditions

l0_in2 = sqrt((pn(npoints/2,1)-pn(npoints/2+2,1))^2+(pn(npoints/2,2)-pn(npoints/2+2,2))^2);     %Initial Length Between Nucleus Points, um
l0_inout=NaN(npoints,1);
if istrap==1                                                                                    %Calculates Distances Beteween Various Points depending on initial Shape
    l0_inout=NaN(npoints,1);
    
for s=1:npoints*ncell
l0_inout(s) = sqrt((nppts(s,1)-nnppts(s,1))^2+(nppts(s,2)-nnppts(s,2))^2);                      %Initial Length Between Cell Points and Nucleus Points, um
end


for s=1:npoints-1
    l0_out=[l0_out sqrt((nppts(s,1)-nppts(s+1,1))^2+(nppts(s,2)-nppts(s+1,2))^2)];              %Initial Length Between Cell Points, um
    l0_in=[l0_in sqrt((nnppts(s,1)-nnppts(s+1,1))^2+(nnppts(s,2)-nnppts(s+1,2))^2)];
end

l0_out=[l0_out sqrt((nppts(npoints,1)-nppts(1,1))^2+(nppts(npoints,2)-nppts(1,2))^2)];
l0_in=[l0_in sqrt((nnppts(npoints,1)-nnppts(1,1))^2+(nnppts(npoints,2)-nnppts(1,2))^2)];
l0_out=repmat(l0_out,1,ncell*3);
l0_in=repmat(l0_in,1,ncell*3);

l0_in=double(min(l0_in));
else
    l0_inout=ones(npoints*ncell,1);
    l0_out=ones(npoints*ncell,1);
   % l0_in=ones(npoints*ncell,1);
    l0_in = l0_in*sqrt((pn(1,1)-pn(2,1))^2+(pn(1,2)-pn(2,2))^2);            %Initial Length Between Nucleus Points, um had scaled down 0.75
    l0_inout(:)=l0_inout*sqrt((pm(1,1)-pn(1,1))^2+(pm(1,2)-pn(1,2))^2);
    l0_out(:) = l0_out*sqrt((pm(2,1)-pm(3,1))^2+(pm(2,2)-pm(3,2))^2);       %Initial Length Between Cell Points, um
end



mindistmembrane=0.1*mean(l0_out);                                           %sets minimum distance between membrane points
l0_nu = nuclrad*2;                                                          %Initial Length Between Internal Nucleus Points, um
l_bond = 0.1;                                                               %Initial Length Between Cell to Cell Binding, um
%ecm_bond = 0.1;                                                             %Initial Length Between Cell to ECM Binding, um
outsidedist=structrad+cellrad*3;


%Initialization for Solver
attached = zeros(1,4*ncell*npoints*3);                                      %Zeros all Cell to Cell Adhesions
F_intx = zeros(ncell*npoints*3, simulation_length/time_step);               %Zeros Force Internal in X-direction, nN
F_inty = zeros(ncell*npoints*3, simulation_length/time_step);               %Zeros Force Internal in Y-direction, nN
F_external = zeros(ncell*npoints*3,1);  %, simulation_length/time_step);    %Zeros Force External, nN
p = zeros(ncell*npoints*4*3, simulation_length/time_step*4);                %Zeros Point Positions, um
new_point = zeros(ncell*npoints*4*3, simulation_length/time_step*4);        %Zeros Point Positions, um
interval=0;                                                                 %Initial Interval
force_number = 1;                                                           %Initial Force Vector
start = 0.1;                                                                %Start Time, s




for s=1:ncell*3
    ECM(s,1)=meanforce;
    ECM(s,2)=SDforce;
end
circularity=cat(3,circularity,circularitytest(x0,ncell,npoints,attached));
F_external  = zeros(ncell*npoints,1);     

initialcellnumber=ncell;

xi=x0;                                                                      %saves original positions in xi
ttt=1;                                                                      %step count for ODE runs, used for data storage

options = odeset('RelTol',1e-2);%,'OutputFcn',@odeplot);                    %Option Settings for solver
timetaken=zeros(length(first_step:time_step:simulation_length),1);          %stores time per step

%%%%%%% sets up video 
if isconforce==1
    basenamemod='Contractile';
else
    basenamemod='Spring';
end
    

filename = [basename,string(file_number)];
Fmov=VideoWriter(strcat(basename,basenamemod,num2str(file_number)));
Fmov.Quality=95;
open(Fmov);

attached =zeros(1,4*ncell*npoints);                                         %Zeros all Cell to Cell Adhesions Attachment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%The Important For Loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for  total_t = first_step:time_step:simulation_length                       %Total Simulation Time And Step Size
    tic
    ttt=1;
    interval=interval+1;                                                    %Next Interval
    OVL=zeros(ncell*npoints,1);
    %%%%%%%%%%%%% adds extra cell if needed. position can be modified by
    %%%%%%%%%%%%% changeing angle 

% % % % %     if extracell==1
% % % % %        %make a new cell
% % % % %     if interval==11
% % % % %         newcell=zeros(npoints,2);
% % % % %         xcenter=(structrad-cellrad*2)*cos((0)*pi/npoints);
% % % % %         ycenter=(structrad-cellrad*2)*sin((0)*pi/npoints);
% % % % %         
% % % % %         newcellm=pointgen([xcenter,ycenter],npoints,cellrad);  
% % % % %         newcelln=pointgen([xcenter,ycenter],npoints,nuclrad);  
% % % % %         newcell=[newcellm,newcelln];
% % % % % 
% % % % %         xnew=addnewcell(x0,newcell,ncell,npoints);
% % % % % 
% % % % %         x0=xnew(:,1);
% % % % % 
% % % % %         ncell=ncell+1;
% % % % %         Cell_Point_Number=zeros(4*ncell*npoints,1);
% % % % %             for n=1:ncell
% % % % %                Cell_Point_Number(((n-1)*npoints+1):((n)*npoints))=1:npoints;
% % % % %                c(((n-1)*npoints+1):((n)*npoints))=n;
% % % % %             end
% % % % %     end 
% % % % % end
% % % % % if alwayssame~=0            %allows change in cell stiffness mid run
% % % % %         if interval==10
% % % % %             k_in(:,1:npoints*2)=0.5.*k_in(:,1:npoints*2);
% % % % %         end
% % % % % end
 
  %%%% Changes Membrane-Nucleus Connections at given time
  if isdelayed~=0&&total_t==isdelayed
     s1=randi([1 ncell-nbroken+1]);
     for s=((s1-1)*npoints+1):((s1+nbroken)*npoints)
                   if rem(s,isbroken)==0
                      k_in(:,s)=k_in(:,s)*brokenscale; 
                   end
     end
  end
  
    internarea(interval)=intpressuretest(x0);                               %defines internal pressure
    ppnts=zeros(npoints,ncell,2);

    isinside=zeros(ncell*3,6);                                              %zeros check for if cells are inside structure
    
    
    %%%% sets up initial area and pressure
    if interval==1
       area=NaN(ncell*3,length( first_step:time_step:simulation_length  )); 
    end
    
    
    area(1:ncell,interval)=getarea(x0,ncell,npoints);
    if rem(total_t,30)==0
    circularity=cat(3,circularity,circularitytest(x0,initialcellnumber,npoints,attached));
    end
        randorder=randperm(npoints*ncell);
    %%%%%%%%%%% First Interval Setup and Run
    if interval==1
        attached =zeros(1,4*ncell*npoints);                                                             %Zeros all Cell to Cell Adhesions Attachment
        interpStruct(interval)= PLumen;                                                                 %stores lumen pressure 
    for n=0:smalltimestep:time_step

%        [t,dxdy] = ode45(@Function_InternalP,[start:0.01:total_t],x0,options);                         %ODE45 Solver
        dxdy=Function_InternalP(start,x0);                                                              %change in position
        x0=x0+dxdy.'*smalltimestep; 
    end
    if n==ncell
    disp('initialized');
    end
    end
        
        
    for n=1:ncell
    %%%%%%%%%%%%%%% All Subsequent Intervals
        if interval>1
            try
        isinside=checkinternalstructure(x0,ncell,npoints);                                              %not in use b/c not adding new cells
            end

        internalP(n,interval)=internalP(n,1).*area(n,1).*thiccness./(area(n,interval).*thiccness);      %Calculates new internal pressure
        end
    end
        if iswhole==1
        interpStruct(interval)=abs(PLumen*internarea(1).*thiccness./(internarea(interval).*thiccness));
        else
        interpStruct(interval)=0;
        end
    
    
    if iseuler==1
        
    for n=start:smalltimestep:total_t
%        [t,dxdy] =
%        ode45(@Function_InternalP,[start:0.01:total_t],x0,options);
%        %ODE45 Solver. for smaller time steps, not really needed otherwise
        dxdy=Function_InternalP(start,x0);   
        
       x0=x0+dxdy.'*smalltimestep;                                                                      %change in position
    end
        
        new_point(1:length(x0),interval) = x0; 
        start=total_t;
    
    else
            for s=1:ncell*npoints                                                                       %scales nucleus-membrane springs if needed
                if isint(s)~=0
                    k_in(s)=basalkscale*k_in(s);
                    k_cin(s)=basalkscale*nucmembrane_stiffness;

                    %k_cin(s)=cin;
                end
                if isbroken ~=0
                   if rem(s,isbroken)==0
                      k_in(s)=0; 
                   end
                end
            end
    [t,newpoints] = ode45(@Function_InternalP,[start:0.01:total_t],x0,options);                         %ODE45 Solver
    x0 = newpoints(size(t,1),:).';                                                                      %Initial Positions Of All Points
    new_point(1:length(x0),interval) = newpoints(size(t,1),:).';                                        %Stores Position Values On Workspace
    start=total_t;
        
        
        
        
        
    end
    if interval ==1
    area(ncell:ncell*3,1)=mean(area(1:ncell,1));                                                        
    end

    if rem(total_t,30)==0                                                                               %Fixes overlap every 30 seconds and order every 60 seconds
        if rem(total_t,60)==0
            fxrd=1;
        else 
            fxrd=0;
        end
        x0=fixOverlap(x0,ncell,npoints,l_bond,fxrd);
        attached=zeros(1,4*ncell*npoints);
    end
    
    if total_t>0
    F_external  = RandForceGen(ECM,point_set);                                                         %Generates ECM Force
    end

   if rem(total_t,videostep)==0||total_t==6                                                            %Stores frames every x seconds
    close all;
    figure1=figure;
 
    newcolors={'#A2142F','#FF0000','#D95319','#EDB120', '#FFFF00','#00FF00','#77AC30','#0072BD','#0000FF','#0072BD'};
    %newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
    colororder(newcolors)
    n=1;
    %%%%% like this to make sure there are no overlaps/errors.
    if rem(total_t,30)==0
        hasfixed=0;
    else
        hasfixed=999;
    end
%     x0t=x0;
while n<=ncell
    if n==1
    colororder(newcolors)
    end
    lastwarn('');
    ppnts(:,n,1)=x0(((n-1)*npoints+1):((n)*npoints));
    ppnts(:,n,2)=x0(ncell*npoints+(n-1)*npoints+1:ncell*npoints +(n)*npoints);
    plot(polyshape(ppnts(:,n,1),ppnts(:,n,2)));

    hold on;
    if ~isempty(lastwarn) && hasfixed<ncell*2
        x0=fixorder(ncell,npoints,x0,n);
            close all;
        figure1=figure;
        %newcolors={'#A2142F','#FF0000','#D95319','#EDB120', '#FFFF00','#00FF00','#77AC30','#0072BD','#0000FF','#0072BD'};
        %colororder(newcolors)
         n=0;
         hasfixed=hasfixed+1;
         close;
         figure;
    end
    %scatter(ppnts(1,n,1),ppnts(1,n,2),'*');        %can be used to show
    %each membrane p
    if n>=1
    ppnts(:,n,1)=x0((2*ncell*npoints+((n-1)*npoints+1)):2*ncell*npoints+((n)*npoints));
    ppnts(:,n,2)=x0((3*ncell*npoints+((n-1)*npoints+1)):3*ncell*npoints+((n)*npoints));
    plot(polyshape(ppnts(:,n,1),ppnts(:,n,2)));
    end
    n=n+1;
end

axis(axisscaler*[-structrad structrad -structrad structrad]); 
timemins=strjoin([string(total_t/60),'mins']);

text(-axisscaler*0.9*structrad,-axisscaler*0.9*structrad,timemins);

ttl=[strjoin(['isbroken', string(isbroken)]), ...
    strjoin(['Conforce', string(Conforce)])];

text(-axisscaler*0.9*structrad,axisscaler*0.9*structrad,ttl);

title(ttl);

frame=getframe;
writeVideo(Fmov,frame);
end
    interval
    timeelapsed=toc
hold off

    timetaken(interval)=timeelapsed;
    currentt=now;
    d = datetime(currentt,'ConvertFrom','datenum')
    if timeelapsed>1500000000000
        break
    end
    start=total_t;
end  
close(Fmov);
close all;
filename = [basename,string(file_number),string(total_t),date];
save(strjoin(filename))

if ttt<0                    %For Testing purposes only. 
    
    
figure;
hold on
for n=1:ncell
circmat(:,n)=squeeze(circularity(3,n,:));
circrad(:,n)=squeeze(circularity(4,n,:));

circstruct(:,n)=squeeze(circularity(3,n,:));
plot(circmat)
end
ylim([0,1.5])
%plot(mean(circmat(:,:)),'g')

%plot(mean(1./circrad))

figure;
hold on
for n=1:ncell
circarea(:,n)=squeeze(circularity(1,n,:));
circperim(:,n)=squeeze(circularity(2,n,:));




%circstruct(:,n)=squeeze(circularity(3,n,:));
C=circperim./(2*pi*sqrt(circarea/pi));

plot(C)
end

end
end
%save('wknd','wknd')