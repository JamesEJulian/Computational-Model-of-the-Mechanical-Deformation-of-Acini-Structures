function dydt = Acini_Function_v14_T1(t,i)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Global Parameters

global ncell                %Number of Cells
global npoints              %Number of Cell Points
global c                    %Cell Number
global Cell_Point_Number    %Point on Cell

global l0_out               %Initial Length Between Cell Points, um
global l0_in                %Initial Length Between Nucleus Points, um
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

global Friction
global v                    %Viscosity Coefficient, nN*s/um

global total_t              %Total Time of Simulation

% % % global ecm_thresh_max       %Max Cell to ECM binding Distance, um
% % % global ecm_thresh_min       %Max Cell to ECM binding Distance, um
% % % global ecm_point            %ECM Point Number
% % % global ecm_x                %ECM X-Position
% % % global ecm_y                %ECM Y-Position
% % % global ecm_bond             %Initial Length Between Cell to ECM Binding, um
% % % global Conforce
global ttt
global internalP
global thiccness
% % % global area
% % % global pressure
global time_step
global interval
global interpStruct
global isint
global intbound
global file_number
global simulation_length
global nucpull
global nucspring
% % % global somwrong
global attached
global cell2cell
global isinside
global pushoutforce
global isconforce
global outsidedist
global K_outerlim
% % % global basedist
global PhaseShift
global mindistmembrane
global ConAmp
global ConPer
global randorder
global outerviscscale
global exConScale
global OVL
global oldT
global isattached
global forceCap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialization, setting global varliables that are regularly used to local
%variables to save time. multiplication b/c it's also slow
ncl=ncell;
nptts=npoints;
expectedist=mindistmembrane;
totalind=ncl*nptts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Force Calculation on Cell Points From Neighboring Cell Points and Nucleus Points
F_neighx=zeros(totalind,1);
F_neighy=zeros(totalind,1);
F_outinx=zeros(totalind,1);
F_outiny=zeros(totalind,1);
d0=zeros(totalind,1);
F_adhx_out=zeros(totalind,1);
F_adhy_out=zeros(totalind,1);
F_extx=zeros(totalind,1);
F_exty=zeros(totalind,1);
dist_n=zeros(totalind,1);
xi=i;
ppnts=zeros(nptts,ncl,2);
%%%% Area and perimeter calculations for each cell, to be used for force
%%%% due to pressure calculations
tot_t=rem(floor(total_t/30),2);
sideforce=zeros(nptts,ncl);
dd=zeros(nptts,ncl);                    %allocate difference matrix
ConforceV=zeros(nptts,ncl,2);
sideunitV=zeros(nptts,ncl,2);
intforce=zeros(nptts,ncl);
intforceV=zeros(nptts,ncl,2);
pushout=zeros(totalind,2);
v = Friction * ones(1,ncell*npoints*3); %Viscosity Coefficient, nN*s/um 20
Fx_out=zeros(totalind,1);
Fy_out=zeros(totalind,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Stores positions in an easier to use matrix. also calculates the center
%%% of each cell
for n=1:ncell
        ppnts(:,n,1)=xi(((n-1)*nptts+1):((n)*nptts));
        ppnts(:,n,2)=xi(totalind+(n-1)*nptts+1:totalind +(n)*nptts);
        xcent(n)=mean(ppnts(:,n,1));
        ycent(n)=mean(ppnts(:,n,2));
end


%attached =zeros(1,4*totalind);        %Zeros all Cell to Cell Adhesions Attachment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First Set of Force Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for sx=1:(totalind)                %Iteration, Points 1 Through (totalind) of Cell Points
     s=randorder(sx);
     
    if isconforce ==0 || isconforce==2      %original elastic calculations
        if Cell_Point_Number(s)>1 && Cell_Point_Number(s)<nptts                                                                             %If Cell Point # 1<CP<npoints use this for loop
        d1(s)=sqrt((xi(s)-xi(s+1))^2+(xi(s+(totalind))-xi(s+(totalind)+1))^2);                                                 %Distance from Cell Point to Cell Point Left
        d2(s)=sqrt((xi(s)-xi(s-1))^2+(xi(s+(totalind))-xi(s+(totalind)-1))^2);                                                 %Distance from Cell Point to Cell Point Right
        F_neighx(s)=F_neighx(s)+k_cout(s)*(d1(s)-l0_out(s))*(xi(s+1)-xi(s))/d1(s)+k_cout(s)*(d2(s)-l0_out(s-1))*(xi(s-1)-xi(s))/d2(s);                       %Force Calculation Between Cell Points X-Direction
        F_neighy(s)=F_neighy(s)+k_cout(s)*(d1(s)-l0_out(s))*(xi(s+(totalind)+1)-xi(s+(totalind)))/d1(s)+...                                 %Force Calculation Between Cell Points Y-Direction
        k_cout(s)*(d2(s)-l0_out(s-1))*(xi(s+(totalind)-1)-xi(s+(totalind)))/d2(s);

        elseif Cell_Point_Number(s) == 1                                                                                                 %If Cell Point # = 1 use this for loop
        d1(s)=sqrt((xi(s)-xi(s+1))^2+(xi(s+(totalind))-xi(s+(totalind)+1))^2);                                                 %Distance from Cell Point to Cell Point Left
        d2(s)=sqrt((xi(s)-xi(s+(nptts-1)))^2+(xi(s+(totalind))-xi(s+(totalind)+(nptts-1)))^2);                             %Distance from Cell Point to Cell Point Right
        F_neighx(s)=F_neighx(s)+k_cout(s)*(d1(s)-l0_out(s))*(xi(s+1)-xi(s))/d1(s)+k_cout(s)*(d2(s)-l0_out(nptts))*(xi(s+(nptts-1))-xi(s))/d2(s);             %Force Calculation Between Cell Points X-Direction
        F_neighy(s)=F_neighy(s)+k_cout(s)*(d1(s)-l0_out(s))*(xi(s+(totalind)+1)-xi(s+(totalind)))/d1(s)+...                                 %Force Calculation Between Cell Points Y-Direction
            k_cout(s)*(d2(s)-l0_out(nptts))*(xi(s+(totalind)+(nptts-1))-xi(s+(totalind)))/d2(s);
        
        elseif Cell_Point_Number(s) == nptts                                                                                           %If Cell Point # = npoints use this for loop
        d1(s)=sqrt((xi(s)-xi(s-(nptts-1)))^2+(xi(s+(totalind))-xi(s+(totalind)-(nptts-1)))^2);                             %Distance from Cell Point to Cell Point Left
        d2(s)=sqrt((xi(s)-xi(s-1))^2+(xi(s+(totalind))-xi(s+(totalind)-1))^2);                                                 %Distance from Cell Point to Cell Point Right
        F_neighx(s)=F_neighx(s)+k_cout(s)*(d1(s)-l0_out(s))*(xi(s-(nptts-1))-xi(s))/d1(s)+k_cout(s)*(d2(s)-l0_out(s-1))*(xi(s-1)-xi(s))/d2(s);
        F_neighy(s)=F_neighy(s)+k_cout(s)*(d1(s)-l0_out(s))*(xi(s+(totalind)-(nptts-1))-xi(s+(totalind)))/d1(s)+...                       %Force Calculation Between Cell Points Y-Direction
            k_cout(s)*(d2(s)-l0_out(s-1))*(xi(s+(totalind)-1)-xi(s+(totalind)))/d2(s);
        end
    end
     
    d0(s)=sqrt((xi(s)-xi(s+(2*totalind)))^2+(xi(s+(totalind))-xi(s+(3*totalind)))^2);       %Distance From Cell Point to Corresponding Nucleus Point
    
    
    %%%% calculates force of membrane-nucleus connection based on type
    %%%% chosen
    if nucspring==1 || nucspring ==2
    ktemp=(k_in(1,s));
    F_outinx(s)=ktemp*(d0(s)-l0_inout(s))*(xi(s+(2*totalind))-xi(s))/d0(s);                          %Force Calculation Between Cell Point and Nucleus Point X-Direction
    F_outiny(s)=ktemp*(d0(s)-l0_inout(s))*(xi(s+(3*totalind))-xi(s+(totalind)))/d0(s);          %Force Calculation Between Cell Point and Nucleus POint Y-Direction
    
    elseif nucspring ==0 || nucspring ==2
    
    dstnc=sqrt((xi(s+(2*totalind))-xi(s))^2+(xi(s+(3*totalind))-xi(s+(totalind)))^2);
    F_outinx(s)=(xi(s+(2*totalind))-xi(s))/dstnc*nucpull+ F_outinx(s);
    F_outiny(s)=(xi(s+(3*totalind))-xi(s+(totalind)))/dstnc*nucpull+F_outiny(s);
    F_inoutx(s)=-F_outinx(s);
    F_inouty(s)=-F_outiny(s);
    end
    

%% Force Calculation on Cell to Cell Adhesion
    
    F_adhx=0;                                   %Zeros all Cell to Cell Adhesions Forces
    F_adhy=0;                                   %Zeros all Cell to Cell Adhesions Forces

        for wx=1:totalind                                             %Iteration, Points 1 through (totalind) of Cell Structure
            w=randorder(wx);
            if tot_t==0                                               %Makes sure correct neighbors are tagged
                if rem(w,nptts)==1
                    w2=w+nptts-1;
                else
                    w2=w-1;
                end
            else                                
                if rem(w,nptts)==0
                    w2=w-nptts+1;
                else
                    w2=w+1;
                end              
            end 
            if c(w)~=c(s)
            [F_adhx,F_adhy,attached,OVL,oldT,isattached]=ADH_Force(c,s,w,w2,xi,totalind,attached,F_adhx,F_adhy,thresh_max,thresh_min,k_cell,l_bond,xcent(c(s)),ycent(c(s)),ttt,OVL,v,oldT,t,isattached);
            end
             if F_adhx~=0 || F_adhy~=0
                 break                                                      %goes to next point once once connection is found
             end
        end
        
    if F_adhx==0
        isattached(s)=0;
    else
        isattached(s)=1;
         v(s)=Friction*15;
    end

    F_adhx_out(s) = F_adhx;                                                 %Set Adhesion Force
    F_adhy_out(s) = F_adhy;                                                 %Set Adhesion Force
    d_center(s)=1;

    %% ECM FORCE
    if attached(s) ~= 0                                     
        
        F_extx (s) = 0;                                                     %Zero Force In X-Direction
        F_exty (s) = 0;                                                     %Zero Force In Y-Direction
    elseif attached(s) == 0    
                nr_x = mean(xi(((2*totalind)+((c(s)-1)*nptts)+1):(2*totalind+((c(s))*nptts))));
                nr_y = mean(xi(((3*totalind)+((c(s)-1)*nptts)+1):(3*totalind+((c(s))*nptts))));
                dist_n(s) = sqrt((xi(s)-nr_x)^2+(xi(s+totalind)-nr_y)^2);
                d_center(s)=sqrt(xi(s)^2+xi(s+totalind)^2);
                if abs(d_center(s))>outsidedist*1.1                         %adds force to hold structure together based on size
                    F_externlim(s)=K_outerlim*(-d_center(s)+outsidedist);
                else
                    F_externlim(s)=0;
                end
        if F_external(s) ~= 0  && isint(s)==0 &&attached(s)==0             %External Force Generated
            d_center(s)=1;
            F_extx (s) = F_external(s)*(xi(s)-nr_x)/dist_n(s)+F_externlim(s)*xi(s)/d_center(s);
            F_exty (s) = F_external(s)*(xi(s+(totalind))-nr_y)/dist_n(s)+F_externlim(s)*xi(s+totalind)/d_center(s);
        elseif isint(s)==0 &&attached(s)==0
            
            F_extx (s) = F_externlim(s)*xi(s)/d_center(s); 
            F_exty (s) = F_externlim(s)*xi(s+totalind)/d_center(s);
            v(s)=v(s)*outerviscscale;
        else 
              F_extx (s) = 0;                                               %Zero Force In X-Direction
              F_exty (s) = 0;                                               %Zero Force In Y-Direction

        end

    else
        
        F_extx (s) = 0;                                                     %Zero Force In X-Direction
        F_exty (s) = 0;                                                     %Zero Force In Y-Direction
        
    end
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% force calculation on membrane
for n=1:ncl
        adhp(:,n,1)=attached(((n-1)*nptts+1):((n)*nptts));
        adhp(:,n,2)=c(((n-1)*nptts+1):((n)*nptts));
end
    for(n=1:ncl)
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Intracellular Pressure
        for(s=1:nptts-1)
         dd(s,n)=sqrt((-ppnts(s,n,1)+ppnts(s+1,n,1))^2 +(-ppnts(s,n,2)+ppnts(s+1,n,2))^2);     
         if adhp(s,n,1)==0||adhp(s+1,n,1)==0
         sideforce(s,n)=internalP(n,interval).*(dd(s,n)*thiccness);
         else
            tempP=internalP(n,interval)-internalP(adhp(s,n,2),interval);
            sideforce(s,n)=tempP.*(dd(s,n)*thiccness);
         end
            ConSin(s,n)=ConAmp*sin(ConPer*t+PhaseShift((n-1)*nptts+s))+ConAmp;      %ConSin is used if it is enabled, not enabled by default so no changes
            %%%%%%%%%%%%%%%%%%%%%%% Lumen Pressure
         if isint((n-1)*nptts+s)~=0 && attached((n-1)*nptts+s)==0 %&& isint((n-1)*nptts+s+1)~=0
                intforce(s,n)=-interpStruct(interval)*dd(s,n)*thiccness;
         else
                intforce(s,n)=0;
         end
            
        end
        
        
        ConSin(s+1,n)=ConAmp*sin(ConPer*t+PhaseShift((n-1)*nptts+s+1))+ConAmp*1.5;
        dd(nptts,n)=sqrt((+ppnts(1,n,1)-ppnts(nptts,n,1))^2 +(+ppnts(1,n,2)-ppnts(nptts,n,2))^2);

        if isint((n-1)*nptts+nptts)~=0 &&attached((n-1)*nptts+s)==0 %&&isint((n-1)*nptts+1)~=0
         intforce(nptts,n)=-interpStruct(interval)*dd(nptts,n)*thiccness;
        end


        for s=1:(nptts-1)
            sideunitV(s,n,1)=(-ppnts(s,n,1)+ppnts(s+1,n,1))/dd(s,n);
            sideunitV(s,n,2)=(-ppnts(s,n,2)+ppnts(s+1,n,2))/dd(s,n);
        end
%% pressure cancels out        
        if adhp(1,n,1)==0||adhp(nptts,n,1)==0
         sideforce(nptts,n)=internalP(n,interval)*(dd(nptts,n)*thiccness);
        else
            tempP=internalP(n,interval)-internalP(adhp(nptts,n,2),interval);
            sideforce(nptts,n)=tempP.*(dd(nptts,n)*thiccness);   
        end
%%          

         sideunitV(nptts,n,1)=(-ppnts(nptts,n,1)+ppnts(1,n,1))/dd(s,n);
         sideunitV(nptts,n,2)=(-ppnts(nptts,n,2)+ppnts(1,n,2))/dd(s,n);    
         tmp=sideunitV(:,n,1);
         sideunitV(:,n,1)=sideunitV(:,n,2);
         sideunitV(:,n,2)=-tmp;

        %find distance between each points. The last point in each cell is
        %different
        for s=2:nptts
         sideforceV(s,n,1)=(sideunitV(s,n,1)*sideforce(s,n))/2+(sideunitV(s-1,n,1)*sideforce(s-1,n))/2;
         sideforceV(s,n,2)=(sideunitV(s,n,2)*sideforce(s,n))/2+(sideunitV(s-1,n,2)*sideforce(s-1,n))/2;
            if isint((n-1)*nptts+s)~=0 || isint((n-1)*nptts+s-1)~=0
                intforceV(s,n,1)=(sideunitV(s,n,1)*intforce(s,n))/2+(sideunitV(s-1,n,1)*intforce(s-1,n))/2;
                intforceV(s,n,2)=(sideunitV(s,n,2)*intforce(s,n))/2+(sideunitV(s-1,n,2)*intforce(s-1,n))/2;      
            end
         end
       for (s=2:nptts-1)

           if dd(s-1,n)>expectedist
             ConforceV(s,n,1)=(-ppnts(s,n,1)+ppnts(s-1,n,1))/dd(s-1,n)* ConSin(s-1,n);
             ConforceV(s,n,2)=(-ppnts(s,n,2)+ppnts(s-1,n,2))/dd(s-1,n)* ConSin(s-1,n);
           end

           if dd(s,n)>expectedist

            ConforceV(s,n,1)= ConforceV(s,n,1)+(-ppnts(s,n,1)+ppnts(s+1,n,1))/dd(s,n)*ConSin(s,n);
            ConforceV(s,n,2)=ConforceV(s,n,2)+(-ppnts(s,n,2)+ppnts(s+1,n,2))/dd(s,n)*ConSin(s,n);

           end
       end
       
       if dd(1,n)>expectedist
         ConforceV(1,n,:)=(-ppnts(1,n,:)+ppnts(2,n,:))./dd(1,n).*ConSin(1,n);
       end

       if dd(nptts,n)>expectedist
         ConforceV(1,n,:)=ConforceV(1,n,:)+(-ppnts(1,n,:)+ppnts(nptts,n,:))./dd(nptts,n).* ConSin(nptts,n);
         ConforceV(nptts,n,:)=(-ppnts(nptts,n,:)+ppnts(1,n,:))./dd(nptts,n).*ConSin(nptts,n);

       end

       if dd(nptts-1,n)>expectedist
          ConforceV(nptts,n,:)=ConforceV(nptts,n,:)+(-ppnts(nptts,n,:)+ppnts(nptts-1,n,:))./dd(nptts-1,n).* ConSin(nptts-1,n);
       end

          sideforceV(1,n,1)=sideunitV(nptts,n,1)*sideforce(nptts,n)/2+(sideunitV(1,n,1)*sideforce(1,n))/2;
          sideforceV(1,n,2)=sideunitV(nptts,n,2)*sideforce(nptts,n)/2+(sideunitV(1,n,2)*sideforce(1,n))/2;
       if isint((n-1)*nptts+1)~=0 
          intforceV(1,n,1)=sideunitV(nptts,n,1)*intforce(nptts,n)/2+(sideunitV(1,n,1)*intforce(1,n))/2;
          intforceV(1,n,2)=sideunitV(nptts,n,2)*intforce(nptts,n)/2+(sideunitV(1,n,2)*intforce(1,n))/2;
       end
    end
            if isinside(n,1)~=0
          pushout(((n-1)*nptts+isinside(n,1)),1)=isinside(n,3)*pushoutforce;
          pushout(((n-1)*nptts+isinside(n,1)),2)=isinside(n,4)*pushoutforce;
            end
    %%%%% puts sum of calculated forces into the correct matrix size        
    s=1;
    while(s<totalind)
        for n=1:ncl
        for ss=1:nptts

             %if attached(s)~=0
                 ConforceV(ss,n,:)=ConforceV(ss,n,:)*0.2;
                 sideforceV(ss,n,:)=sideforceV(ss,n,:)*1;
             %end
             if attached(s)==0 && isint(s)
                 ConforceV(ss,n,:)=ConforceV(ss,n,:)*exConScale;
             end    
             if isconforce==1||isconforce==2
              
             F_neighxCON(s)=ConforceV(ss,n,1)+sideforceV(ss,n,1)+intforceV(ss,n,1);
             F_neighyCON(s)=ConforceV(ss,n,2)+sideforceV(ss,n,2)+intforceV(ss,n,2);
             else
             F_neighxCON(s)=intforceV(ss,n,1);
             F_neighyCON(s)=intforceV(ss,n,2);
             end
            s=s+1;
            if s>totalind
                break
            end

        end
            if s>totalind
                break
            end
        end
        if s>totalind
                break
        end
    end




%% Force Calculation On Nucleus Points From Neighboring Nucleus Points and Cell Points
% % % % d1_in=zeros(totalind,1);
% % % % dm1_in=zeros(totalind,1);
% % % % d2_in=zeros(totalind,1);
% % % % dm2_in=zeros(totalind,1);
F_inneighx=zeros(totalind,1);
F_inneighy=zeros(totalind,1);

for q=1:(totalind) %Iteration, Points 1 Through (totalind) of Nucleus
    
    
    
    
    
   if  Cell_Point_Number(q) > 2 && Cell_Point_Number(q) < nptts-1                                                               %If Nucleus Point # 1 < CP < nptts Use This For Loop
        
    p1=q+1;
    p2=q+2;
    m1=q-1;
    m2=q-2;
    IO1=l0_in;
    IO2=l0_in2;
    k_cins=k_cin(q);
    
    [F_inneighx(q),F_inneighy(q)]=NMforces(xi,totalind,q,p1,p2,m1,m2,IO1,IO2,k_cins);
        
    ttt;
    elseif Cell_Point_Number(q) == 1                                                                                        %If Nucleus Point # = 1 Use This For Loop
        
    p1=q+1;
    p2=q+2;
    m1=q+(nptts-1);
    m2=q+(nptts-2);
    IO1=l0_in;
    IO2=l0_in2;
    k_cins=k_cin(q);
    
    [F_inneighx(q),F_inneighy(q)]=NMforces(xi,totalind,q,p1,p2,m1,m2,IO1,IO2,k_cins);
      ttt;  
    elseif Cell_Point_Number(q) == 2                                                                                        %If Nucleus Point # = 2 Use This For Loop

    p1=q+1;
    p2=q+2;
    m1=q-1;
    m2=q+(nptts-2);
    IO1=l0_in;
    IO2=l0_in2;
    k_cins=k_cin(q);
    
    [F_inneighx(q),F_inneighy(q)]=NMforces(xi,totalind,q,p1,p2,m1,m2,IO1,IO2,k_cins);   

      elseif Cell_Point_Number(q) == nptts                                                                                  %If Nucleus Point # = nptts Use This For Loop
  
    p1=q-(nptts-1);
    p2=q-(nptts-2);
    m1=q-1;
    m2=q-2;
    IO1=l0_in;
    IO2=l0_in2;
    k_cins=k_cin(q);
    
    [F_inneighx(q),F_inneighy(q)]=NMforces(xi,totalind,q,p1,p2,m1,m2,IO1,IO2,k_cins);

   
    elseif Cell_Point_Number(q) == nptts-1  
    
    p1=q+1;
    p2=q-(nptts-2);
    m1=q-1;
    m2=q-2;
    IO1=l0_in;
    IO2=l0_in2;
    k_cins=k_cin(q);
    
    [F_inneighx(q),F_inneighy(q)]=NMforces(xi,totalind,q,p1,p2,m1,m2,IO1,IO2,k_cins);

    end

    
    d0(q)=sqrt((xi(q)-xi(q+(2*totalind)))^2+(xi(q+(totalind))-xi(q+(3*totalind)))^2);           %Distance From Cell Point to Corresponding Nucleus Point
%     if d0(q)<basedist    
%     ktemp=abs(k_in(2,q))*(exp(-5*d0(q)))+abs(k_in(2,q));
%     else
    ktemp=(k_in(1,q));
    %end
    
    if nucspring==1 || nucspring ==2
   %     F_inoutx(q)=ktemp*(d0(q)-l0_inout(q))*(xi(q)-xi(q+(2*totalind)))/d0(q)+Finoutx(q);                              %Force Calculation Between Cell Point and Nucleus Point X-Direction
   %     F_inouty(q)=ktemp*(d0(q)-l0_inout(q))*(xi(q+(totalind))-xi(q+(3*totalind)))/d0(q)+Finouty(q);              %Force Calculation Between Cell Point and Nucleus Point Y-Direction
        F_inoutx(q)=-F_outinx(q);
        F_inouty(q)=-F_outiny(q);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Force Calculation On Internal Nucleus
    

    if  Cell_Point_Number(q) >= 1 && Cell_Point_Number(q) <= nptts/2                                                         %If Nucleus Point Is >=1 & <=10 Use This For Loop
        
        d_nu(q)=sqrt((xi(q+(2*totalind))-xi(q+(2*totalind)+nptts/2))^2+...                                           %Distance calculation Between Nucleus Points
            (xi(q+(2*totalind)+(totalind))-xi(q+(2*totalind)+(totalind)+nptts/2))^2);
        
        F_nx(q)=k_nu(q)*(d_nu(q)-l0_nu)*(xi(q+(2*totalind)+nptts/2)-xi(q+(2*totalind)))/d_nu(q);   %Force Calculation Between Nucleus Points X-Direction
        F_ny(q)=k_nu(q)*(d_nu(q)-l0_nu)*(xi(q+(3*totalind)+nptts/2)-xi(q+(3*totalind)))/d_nu(q);   %Force Calculation Between Nucleus Points Y-Direction
        
    elseif Cell_Point_Number(q) >= nptts/2+1 && Cell_Point_Number(q) <= nptts                                                     %If Nucleus Point Is >=11 & <=20 Use This For Loop
        
        d_nu(q)=sqrt((xi(q+(2*totalind))-xi(q+(2*totalind)-nptts/2))^2+...                                           %Distance calculation Between Nucleus Points
            (xi(q+(2*totalind)+(totalind))-xi(q+(2*totalind)+(totalind)-nptts/2))^2);
        
        F_nx(q)=k_nu(q)*(d_nu(q)-l0_nu)*(xi(q+(2*totalind)-nptts/2)-xi(q+(2*totalind)))/d_nu(q);   %Force Calculation Between Nucleus Points X-Direction
        F_ny(q)=k_nu(q)*(d_nu(q)-l0_nu)*(xi(q+(3*totalind)-nptts/2)-xi(q+(3*totalind)))/d_nu(q);   %Force Calculation Between Nucleus Points Y-Direction
        
    end
    

    F_inadhx=0;       %Zero Force In X-Direction
    F_inadhy=0;       %Zero Force In Y-Direction
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Total Nucleus Force Calculation

    Fx_in(q)= F_inneighx(q)+F_inoutx(q)+F_inadhx+F_intx(q)+F_nx(q);    %Sum Of Forces In X-Direction
    Fy_in(q)= F_inneighy(q)+F_inouty(q)+F_inadhy+F_inty(q)+F_ny(q);    %Sum Of Forces In X-Direction
   
   %%% Total Cell Point Force Calculation
     Fx_out(q)= F_neighx(q)+F_neighxCON(q)+F_outinx(q)+F_adhx_out(q)+F_extx(q)+pushout(q,1);    %Sum Of Forces In X-Direction
     Fy_out(q)= F_neighy(q)+F_neighyCON(q)+F_outiny(q)+F_adhy_out(q)+F_exty(q)+pushout(q,2);    %Sum Of Forces in Y-Direction
end

     captestin=sqrt(Fx_in.^2+Fy_in.^2);
     captestout=sqrt(Fx_out.^2+Fy_out.^2);
     CTI=find(captestin>forceCap);
     CTO=find(captestout>forceCap);
     
     for q=1:length(CTI)
        Fx_in(CTI(q))=Fx_in(CTI(q))/abs(captestin(CTI(q)))*forceCap; 
        Fy_in(CTI(q))=Fy_in(CTI(q))/abs(captestin(CTI(q)))*forceCap; 
     end
     
     for q=1:length(CTO)
        Fx_out(CTO(q))=Fx_out(CTO(q))/abs(captestout(CTO(q)))*forceCap; 
        Fy_out(CTO(q))=Fy_out(CTO(q))/abs(captestout(CTO(q)))*forceCap; 
     end



%% final force calculations and output
dydt = zeros(totalind*4, 1);                     %Zero Integration Solutions

for d = 1:totalind                             %Iteration, Points 1 Through (totalind) of Cell Points
    
    dydt(d) = Fx_out(d)/v(d).';                     %Solves New Cell Point Position For Points 1 Through nptts in X-Direction
    dydt(d+totalind) = Fy_out(d)/v(d).';       %Solves New Cell Point Position For Points 1 Through nptts in Y-Direction
    dydt(d+2*totalind) = Fx_in(d)/((v(d).'*0.6));      %Solves New Nucleus Position For Points 1 Through nptts in X-Direction
    dydt(d+3*totalind) = Fy_in(d)/((v(d).'*0.6));      %Solves New Nucleus Position For Points 1 Through nptts in Y-Direction
    
end
ttt=ttt+1;      %itterates ODE step number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% troubleshooting stuff Not normally run
if ttt<0
    %%% testing membrane forces
    
    
    figure
hold on
   for n=1:ncl
    ppnts(:,n,1)=xi(((n-1)*npoints+1):((n)*npoints));
    ppnts(:,n,2)=xi(ncell*npoints+(n-1)*npoints+1:ncell*npoints +(n)*npoints);
       
   plot(polyshape(ppnts(:,n,1),ppnts(:,n,2)));
   scatter(ppnts(npoints,n,1),ppnts(npoints,n,2),'o');
   scatter(ppnts(1,n,1),ppnts(1,n,2),'*');
   
    axis equal
   end
   
   for n=1
    %plot(ppnts(:,n,1),ppnts(:,n,2)); 
    %ssd=sideunitV(1:19,n,:)+sideunitV(2:20,n,:);
    %ssd(20,n,:)=sideunitV(1,n,:)+sideunitV(20,n,:);
     %quiver(ppnts(:,n,1),ppnts(:,n,2),ssd(:,n,1),ssd(:,n,2));
     hold on
     %quiver(ppnts(:,n,1),ppnts(:,n,2),sideforceV(:,n,1)+ConforceV(:,n,1),sideforceV(:,n,2)+ConforceV(:,n,2));
     %quiver(ppnts(:,n,1),ppnts(:,n,2),ConforceV(:,n,1),ConforceV(:,n,2));
     
     quiver(ppnts(:,n,1),ppnts(:,n,2),sideforceV(:,n,1),sideforceV(:,n,2));
     quiver(ppnts(:,n,1),ppnts(:,n,2),ConforceV(:,n,1),ConforceV(:,n,2));
     %quiver(ppnts(:,n,1),ppnts(:,n,2),ConunitV(:,n,3),ConunitV(:,n,4));
     %quiver(ppnts(:,n,1),ppnts(:,n,2),intforceV(:,n,1),intforceV(:,n,2));
     %quiver(ppnts(:,n,1),ppnts(:,n,2),sideunitV(:,n,1),sideunitV(:,n,2));
     hold on
     %quiver(ppnts(2:nptts,n,1),ppnts(2:nptts,n,2),sideunitV(1:(nptts-1),n,1),sideunitV(1:(nptts-1),n,2));
     
     
   end
   
   
   
hold on;
for s=1:ncell*npoints
    if isint(s)~=0
    scatter(xi(s),xi(s+ncell*npoints),'*');
    scatter(isint(s),isint(s+ncell*npoints),'o');
    end
end
plot(polyshape(intbound(:,1),intbound(:,2)));
 



 figure;
newcolors={'#A2142F','#FF0000','#D95319','#EDB120', '#FFFF00','#00FF00','#77AC30','#0072BD','#0000FF','#0072BD'};
%newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
colororder(newcolors)
   hold on
   for n=1:ncl
        ppnts(:,n,1)=xi(((n-1)*npoints+1):((n)*npoints));
        ppnts(:,n,2)=xi(ncell*npoints+(n-1)*npoints+1:ncell*npoints +(n)*npoints);
       % plot(polyshape(ppnts(:,n,1),ppnts(:,n,2)));
 %        xfsum=sum(dydt(((n-1)*npoints+1):((n)*npoints)));
 %        yfsum=sum(dydt(ncell*npoints+(n-1)*npoints+1:ncell*npoints +(n)*npoints));
 %        xfsum=sum(sideforceV(:,n,1));
  %       yfsum=sum(sideforceV(:,n,2));
  %       xfsum=sum(ConforceV(:,n,1));
  %       yfsum=sum(ConforceV(:,n,2));
%        xfsum=sum(intforceV(:,n,1));
%        yfsum=sum(intforceV(:,n,2)); 
%        xfsum=sum(intforceV(:,n,1));
%        yfsum=sum(intforceV(:,n,2)); 
        
 %   xfsum=sum(F_inneighx((((n-1)*npoints+1)):((n)*npoints)));
 %   yfsum=sum(F_inneighy((((n-1)*npoints+1)):((n)*npoints)));
%

        xcent=mean(ppnts(:,n,1));
        ycent=mean(ppnts(:,n,2));
        quiver(xcent,ycent,xfsum,yfsum,'linewidth',5);
        %quiver(0,0,xfsum,yfsum,'linewidth',5);
   
   
   
   end
   
   figure;
   hold on;
   for n=1:ncl
   plot(polyshape(ppnts(:,n,1),ppnts(:,n,2)));
   end
   
   for s=1:totalind
        if attached(s)~=0
         plotlocx=[xi(s),xi(attached(s))];
         plotlocy=[xi(s+totalind),xi(attached(s)+totalind)];
         plot(plotlocx,plotlocy);
         scatter(xi(s),xi(s+totalind),'*');
         scatter(xi(attached(s)),xi(attached(s)+totalind),'*');   
        end
     
       
   end
 %%% testing from failure state in solver
   figure;
   hold on;
   for n=1:ncell
    ppnts(:,n,1)=x0(((n-1)*npoints+1):((n)*npoints));
    ppnts(:,n,2)=x0(ncell*npoints+(n-1)*npoints+1:ncell*npoints +(n)*npoints);
 
   plot(polyshape(ppnts(:,n,1),ppnts(:,n,2)));
   end
   
   for s=1:ncell*npoints
        if attached(s)~=0
         plotlocx=[x0(s),x0(attached(s))];
         plotlocy=[x0(s+ncell*npoints),x0(attached(s)+ncell*npoints)];
         plot(plotlocx,plotlocy);
        scatter(x0(s),x0(s+ncell*npoints),'*');
            
        end
     
       
   end
%%% testing nucleus points
  

figure;
   hold on;
  
   for n=1:ncl
    ppnts(:,n,1)=xi(((n-1)*npoints+1):((n)*npoints));
    ppnts(:,n,2)=xi(ncell*npoints+(n-1)*npoints+1:ncell*npoints +(n)*npoints);
       
   plot(polyshape(ppnts(:,n,1),ppnts(:,n,2)));

   end
   for n=1:ncl
    ppnts(:,n,1)=xi((2*ncell*npoints+((n-1)*npoints+1)):2*ncell*npoints+((n)*npoints));
    ppnts(:,n,2)=xi((3*ncell*npoints+((n-1)*npoints+1)):3*ncell*npoints+((n)*npoints));
    plot(polyshape(ppnts(:,n,1),ppnts(:,n,2)));
   end
   axis equal
%    
%     scatter(xi(1),xi(1+totalind));
%     scatter(xi(151),xi(151+totalind));
for s=1:nptts
    %quiver(xi(s),xi(s+totalind),F_outinx(s),F_outiny(s));
    %quiver(xi(s+totalind*2),xi(s+totalind*3),F_inoutx(s),F_inouty(s));
    %quiver(xi(s+totalind*2),xi(s+totalind*3),F_inneighx(s),F_inneighy(s));            %total nucleus force
    %quiver(xi(s+totalind*2),xi(s+totalind*3),Fx_in(s),Fy_in(s));            %total nucleus force
    %quiver(xi(s),xi(s+totalind),F_neighx(s),F_neighy(s));            %membrane force
    %quiver(xi(s),xi(s+totalind),F_neighxCON(s),F_neighyCON(s));            %membrane force
    %quiver(xi(s),xi(s+totalind),Fx_out(s),Fy_out(s));            %membrane force
    %%F_neighx(s)+F_outinx(s)+F_adhx_out(s)+F_extx(s)+pushout(s,1);
    quiver(xi(s),xi(s+totalind),F_outinx(s),F_outiny(s));
    %quiver(xi(s),xi(s+totalind),F_adhx_out(s),F_adhy_out(s));
    
    %quiver(xi(s),xi(s+totalind),F_extx(s),F_exty(s));            %membrane force
    %quiver(xi(s),xi(s+totalind),dydt(s),dydt(totalind+s));

    
end


   s=1;
for n=1:ncell             
for ss=1:npoints
    con(s)=sqrt(ConforceV(ss,n,1)^2+ConforceV(ss,n,2)^2);
    side(s)=sqrt(sideforceV(ss,n,1)^2+sideforceV(ss,n,2)^2);
    int(s)=sqrt(intforceV(ss,n,1)^2+intforceV(ss,n,2)^2);
    neigh(s)=sqrt(F_neighx(s)^2+F_neighy(s)^2);
    outin(s)=sqrt(F_outinx(s)^2+F_outiny(s)^2);
    s=s+1;
end
end
figure;
hold on;
plot(con,'b');
plot(side,'r');
%plot(int,'k');
plot(neigh,'g');
plot(outin,'c');



for s=1:totalind
    if   isint(s)==0 &&attached(s)==0  
    scatter (xi(s),xi(s+totalind),'*');
    end
   
    
end
cell2cell=attached;
if interval==simulation_length/time_step-1 && rem(ttt,20)==1
save(strjoin(['forces',string(file_number)]));
end

end
cell2cell=attached;
if interval==simulation_length/time_step-1 && rem(ttt,20)==1
save(strjoin(['forces',string(file_number)]));
end
end 
