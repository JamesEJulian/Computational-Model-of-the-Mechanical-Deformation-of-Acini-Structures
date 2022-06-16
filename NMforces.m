function [Fx, Fy]=NMforces(xi,tot,base,p1,p2,m1,m2,IO1,IO2,k_cin)
    %Calculates Forces Between Nucleus and Cell Membrane
    global ttt
    global l0_in 
    global l0_in2 
    xoff=tot*2;     %offset value for x positions in matrix
    yoff=tot*3;     %offset value for y positions in matrix
    scale=0.80;     %minimum distance between point scalar
    IN=[base+xoff, base+yoff];  %holds initial points' locations
    P=[p1+xoff, p1+yoff; p2+xoff, p2+yoff; m1+xoff, m1+yoff; m2+xoff, m2+yoff]; %holds all connected points' locations
    %%%% defines initial length if first time called
    if ttt==1
        l0_in = 0.8*sqrt((xi(P(1,1))-xi(IN(1)))^2+(xi(P(1,2))-xi(IN(2)))^2);
        l0_in2 = 0.8*sqrt((xi(P(2,1))-xi(IN(1)))^2+(xi(P(2,2))-xi(IN(2)))^2);
        IO1=l0_in;
        IO2=l0_in2;
    end

%%% distances betweeen points 
     d1_in=sqrt((xi(P(1,1))-xi(IN(1)))^2+(xi(P(1,2))-xi(IN(2)))^2);
     d2_in=sqrt((xi(P(2,1))-xi(IN(1)))^2+(xi(P(2,2))-xi(IN(2)))^2);
     dm1_in=sqrt((xi(P(3,1))-xi(IN(1)))^2+(xi(P(3,2))-xi(IN(2)))^2);
     dm2_in=sqrt((xi(P(4,1))-xi(IN(1)))^2+(xi(P(4,2))-xi(IN(2)))^2);
     
 %%% calculates forcs in the direction of each point. zeros if too small
 %%% distance
    if d1_in<(IO1*scale)
       d1x=0;
       d1y=0;
    else
       d1x=k_cin*(d1_in-IO1)*(-xi(IN(1))+xi(P(1,1)))/d1_in;
       d1y=k_cin*(d1_in-IO1)*(-xi(IN(2))+xi(P(1,2)))/d1_in;
    end

    if d2_in<IO2*scale
       d2x=0;
       d2y=0;
    else
       d2x=k_cin*(d2_in-IO2)*(-xi(IN(1))+xi(P(2,1)))/d2_in;
       d2y=k_cin*(d2_in-IO2)*(-xi(IN(2))+xi(P(2,2)))/d2_in;
    end
    
    if dm1_in<IO1*scale
       dm1x=0;
       dm1y=0;
    else
       dm1x=k_cin*(dm1_in-IO1)*(-xi(IN(1))+xi(P(3,1)))/dm1_in;
       dm1y=k_cin*(dm1_in-IO1)*(-xi(IN(2))+xi(P(3,2)))/dm1_in;
    end
    
    if dm2_in<IO2*scale
       dm2x=0;
       dm2y=0;
    else
       dm2x=k_cin*(dm2_in-IO2)*(-xi(IN(1))+xi(P(4,1)))/dm2_in;
       dm2y=k_cin*(dm2_in-IO2)*(-xi(IN(2))+xi(P(4,2)))/dm2_in;
    end
%%% adds all forces in X and Y
Fx=+d1x+d2x+dm1x+dm2x;
Fy=+d1y+d2y+dm1y+dm2y; 


%% troubleshooting stuff
if xoff==1

    
    scatter(xi(IN(1)),xi(IN(2)),'*');
    scatter(xi(P(1,1)),xi(P(1,2)),'O');
    scatter(xi(P(2,1)),xi(P(2,2)),'O');
    scatter(xi(P(3,1)),xi(P(3,2)),'O');
    scatter(xi(P(4,1)),xi(P(4,2)),'O');
    
   quiver(xi(IN(1)),xi(IN(2)),d1x,d1y);
   quiver(xi(IN(1)),xi(IN(2)),d2x,d2y);
   quiver(xi(IN(1)),xi(IN(2)),dm1x,dm1y);
   quiver(xi(IN(1)),xi(IN(2)),dm2x,dm2y);
     
   quiver(xi(IN(1)),xi(IN(2)),d1x+dm1x,d1y+dm1y);
   quiver(xi(IN(1)),xi(IN(2)),d2x+dm2x,d2y+dm2y);
   
   quiver(xi(IN(1)),xi(IN(2)),d1x/(k_cin*(d1_in-IO1)),d1y/(k_cin*(d1_in-IO1)));
   quiver(xi(IN(1)),xi(IN(2)),dm1x/(k_cin*(dm1_in-IO1)),dm1y/(k_cin*(dm1_in-IO1)));
   
   
   
   %quiver(xi(IN(1)),xi(IN(2)),Fx,Fy);
    
    %sqrt((ppnts(base,n,1)-ppnts(base+1,n,1))^2+(ppnts(base,n,2)-ppnts(base+1,n,2))^2)
   
   
end
end