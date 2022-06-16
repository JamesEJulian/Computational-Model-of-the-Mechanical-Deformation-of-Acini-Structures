function [F_adhx, F_adhy,attached, OVL,oldT,isattached]=ADH_Force(c,s,w,w2,xi,totalind,attached,F_adhx,F_adhy,thresh_max,thresh_min,k_cell,l_bond,xcent,ycent, ttt, OVL,v,oldT,t,isattached)
%Checks if two points should be connected, then calculates the force           
%if c(w)~=c(s)            %If Points Not From Same Cell and Are Not Attached Proceed
                
if isattached(s)==1
    TM=thresh_max*3;
else
    TM=thresh_max;
end

OVS=-0.5;     %scaler for readjusting force to minimize overlap
cap=0.2;
                ddist=(sqrt((xi(s)-xi(w))^2+(xi(s+totalind)-xi(w+totalind))^2));    %Distance From Cell Point To All Other Points of Sorrounding Cell Points
                 if ddist<TM  && ddist>thresh_min                                   %If Cell Points Within Distance Threshold Proceed
                     dd=abs((xi(w2+totalind)-xi(w+totalind))*xi(s)-(xi(w2)-xi(w))*xi(s+totalind)+xi(w2)*xi(w+totalind)-xi(w2+totalind)*xi(w))/sqrt((xi(w2+totalind)-xi(w+totalind))^2+(xi(w2)-xi(w))^2);
                                                                                    %dd= distance between s and nearest point on neighbor side 


                     
                     attached(s)=w;                                                 %Cell Point S Becomes Attached
                         KC=k_cell(s);

                    v1=[xi(w2)-xi(w),xi(w2+totalind)-xi(w+totalind)];               %Vector between two neighbor points
                    v2=[xi(s)-xi(w),xi(s+totalind)-xi(w+totalind)];                 %Vector between target and neighbor point

                    v3=v2-(dot(v1,v2)/dot(v1,v1))*v1;                               %Vector between target and neighbor segment 

                                                                                    %OVL is used to tag points as overlapping another cell
                if or(OVL(s)==1,oldT~=floor(t*100))                                 %checks if either first step or point is tagged
                    if or(ttt==1,oldT~=floor(t*100))                                %resets OVL tags every X steps
                        OVL(s)=0;                                                   %sets OVL to zero
                        oldT=floor(t);
                        F_adh=-KC*(dd-l_bond)/l_bond;                               %Calculates adhesion force based on spring

                        F_adhx=F_adh*v3(1)/norm(v3);                                %Adhesion Forces in X-direction
                        F_adhy=F_adh*v3(2)/norm(v3);                                %Adhesion Forces in Y-direction
                    else
                        F_adh=-KC*(dd-l_bond)/l_bond;
                        F_adh=abs(F_adh)*OVS;
                        if abs(F_adh)>KC*cap
                            F_adh=(F_adh/abs(F_adh))*KC/cap; 
                        end
                        F_adhx=F_adh*v3(1)/norm(v3);                                %Adhesion Forces in X-direction
                        F_adhy=F_adh*v3(2)/norm(v3);                                %Adhesion Forces in Y-direction
                    end
                else
                   F_adh=-KC*(dd-l_bond)/l_bond;
                   if abs(F_adh)>KC/cap
                      F_adh=(F_adh/abs(F_adh))*KC/cap; 
                   end
                   F_adhx=F_adh*v3(1)/norm(v3);                                     %Adhesion Forces in X-direction
                   F_adhy=F_adh*v3(2)/norm(v3);                                     %Adhesion Forces in Y-direction
                end
            else
                F_adhx=F_adhx+0;                                                    %Adhesion Force Unchanged
                F_adhy=F_adhy+0;

                 end
end