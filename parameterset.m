%% Set parameters
simulation_length = 600;                %Length of Simulation, s
%28800 = 8 hours
%3600 = 1 hour
time_step = 6;                              %Time Step of Simulation, s
smalltimestep = 0.1;                        %micro step 
first_step = time_step;                     %Initial Time Step, s
videostep=time_step*10;

meanforce=20;                               %average ECM force, original value is 20
SDforce=2;                                  %standard deviation of ECM force        
pushoutforce=40;                            %force pushing interstructural cells out                               
iseuler=1;                                  %euler approximation Y=1, N=0
% if rem(file_number,2)==1                  %for changing between runs, not used
                                            % for final version 
isbroken=2;                                 %0=normal, else value = 1/n connections set to 0.3x;
nbroken=4;                                  %Number of broken Cells
brokenscale=0.01;
isdelayed=3600*2;                           %0=no, else length of delay, s
% else
% isbroken=0;                               %0=normal, else value = 1/n connections set to 0.3x;
% nbroken=0;                                %Number of broken Cells
% isdelayed=0;                              %0=no, else length of delay, s    
%     
%     
%  end


forceCap=36;                                %maximum force applied at any point

%%% Lengths

l_bond = 0.05;                              %Initial Length Between Cell to Cell Binding, um 0.1

thresh_min = 0.01;                          %Min Cell to Cell binding Distance, um 0.01
thresh_max = 0.20;                          %Max Cell to Cell binding Distance, um 1.25
ecm_thresh_max = 7;                         %Max Cell to ECM binding Distance, um
ecm_thresh_min = 1;                         %Min Cell to ECM binding Distance, um
thiccness=10;                               %simulated thickness of the cell

%%% Stiffness Values nN/um
membrane_stiffness=12;                      %stiffness of cell membrane
nucmembrane_stiffness=36;                   %stiffness of nuclear membrane
basemembrane_nucleus=6;                     %membrane to nucleus connection stiffness      
intnuc_stiffness =24;                       %internal nucleus stiffness  nN/um

cell_cell =24;                              %cell-cell adhesion stiffness


%%% Other Values
Friction=80;                                %nN*s/um
outerviscscale=5;                           %friction scale for outside points
K_outerlim=24;                             %stiffness of external membrane holding structure in shape
%%% Contractile Parameters

startPressure=0.4;                          %initial intracellular pressure
Conforce=30;                                %contractile force
PLumen=0.2;                                 %initial lumen pressure

ConPeriod=0.300;                            %period of conforce, seconds
exConScale=1;                               %scales external cell contractility relative to internal
%%% Pa*10^-3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Values Based on parameters. nothing past here should be edited 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Membrane_nucleus=findingmem2nuc(basemembrane_nucleus,cellrad-nuclrad);      %Makes Membrane-nucleus conneciton very stiff if they get to close together


k_cout = membrane_stiffness * ones(1,ncell*npoints*3);                      %Stiffness Value of Cell Membrane, nN/um was 40
k_cin  = nucmembrane_stiffness  * ones(1,ncell*npoints*3);                  %Stiffness Value of Nucleus Membrane, nN/um 14
basalkscale=0.5;                                                            %can be used to scale down the spring stiffness in basal membrane to nucleus connections
%k_in=NaN(2,ncell*npoints*3);
k_in(2,:) = Membrane_nucleus * ones(1,ncell*npoints*3);                     %Stiffness Value of Cell Membrane to Nucleus, nN/um 15
k_in(1,:)= basemembrane_nucleus * ones(1,ncell*npoints*3);

%%% Calculates the stiffness values if structure is broken initially
 if isbroken ~=0 && isdelayed==0
     s1=randi([1 ncell-nbroken+1]);
%      xi(((n-1)*nptts+1):((n)*nptts))
     for s=((s1-1)*npoints+1):((s1+nbroken)*npoints)
                   if rem(s,isbroken)==0
                      k_in(:,s)=k_in(:,s)*brokenscale; 
                   end
     end
 end
% s1=1;
k_nu =  intnuc_stiffness * ones(1,ncell*npoints*3);     %Stiffness Value of Internal Nucleus, nN/um  14

k_cell=cell_cell*ones(ncell*npoints*3,1);
v = Friction * ones(1,ncell*npoints*3);                                                   %Viscosity Coefficient, nN*s/um 20
%v((ncell*npoints*2+1):ncell*npoints*4)=Friction*10;
point_set=zeros(ncell*npoints*3,1);
intspring=Membrane_nucleus;


%nucpull=Conforce*0.2;       %contractile element of membrane to nucleus connection

ConAmp=Conforce*2/3;
%ConAmp=0;
ConPer=2*pi/ConPeriod;
PhaseShift=ones(ncell*3*npoints,1);
for n=1:ncell*npoints*3
   PhaseShift(n)=rand(1)*pi*2;
   PhaseShift(n)=0;
end


scalemean=0.4;
scalesd=0.15;




initialPressure=startPressure;

internalP=zeros(ncell*3,length(first_step:time_step:simulation_length));
ConforceV=zeros(npoints,ncell*3,2);     %contraction force on each point
internalP(:,1)=startPressure;
%pressure=zeros(ncell*3,length(first_step:time_step:simulation_length));
%pressure(:,1)=startPressure;
circularity=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initialization for Solver
attached = zeros(1,4*ncell*npoints*3);                                 %Zeros all Cell to Cell Adhesions
F_intx = zeros(ncell*npoints*3, simulation_length/time_step);          %Zeros Force Internal in X-direction, nN
F_inty = zeros(ncell*npoints*3, simulation_length/time_step);          %Zeros Force Internal in Y-direction, nN
F_external = zeros(ncell*npoints*3,1);  %, simulation_length/time_step);      %Zeros Force External, nN
p = zeros(ncell*npoints*4*3, simulation_length/time_step*4);           %Zeros Point Positions, um
new_point = zeros(ncell*npoints*4*3, simulation_length/time_step*4);   %Zeros Point Positions, um
interval=0;                                                          %Initial Interval
force_number = 1;                                                    %Initial Force Vector
start = 0.1;                                                         %Start Time, s

