function[basepoints]=traps(ncell, npoints,structrad,rad,isnuc,nucrat)
%Generates Trapezoidal Points 
%structrad=(2*cellrad)/(2*sin(pi/ncell))*1;  %1.15
rotateang=2*pi/(ncell);
rotatemat=[cos(rotateang),-sin(rotateang);sin(rotateang),cos(rotateang)];
% outerstructrad=structrad+cellrad;
% innerstructrad=structrad-cellrad;
% outerperm=2*outerstructrad*pi;
% innerperm=2*innerstructrad*pi;

structure=[structrad+rad,structrad-rad];

if isnuc==0
    angle=[2*pi/ncell,-2*pi/ncell].*0.95;
else
    angle=[2*pi/ncell,-2*pi/ncell]*nucrat;
end
angle=angle./2+[pi/2,pi/2];
% coords=[structure(1)*cos(angle(2)),structure(1)*sin(angle(2));...
%         structure(1)*cos(angle(1)),structure(1)*sin(angle(1));...
%         structure(2)*cos(angle(1)),structure(2)*sin(angle(1));...
%         structure(2)*cos(angle(2)),structure(2)*sin(angle(2))];
angles=linspace(angle(2),angle(1),((npoints-4)/4+2));

basepoints=NaN(npoints*ncell,2);
% basepoints(1,:)=coords(1,:);
% basepoints(npoints/4+1,:)=coords(2,:);


for s=1:((npoints-4)/4+2)
    basepoints(s,:)=[structure(1)*cos(angles(s)),structure(1)*sin(angles(s))];
     basepoints(s+2*(npoints-4)/4+2,:)=[structure(2)*cos(angles((npoints-4)/4+3-s)),structure(2)*sin(angles((npoints-4)/4+3-s))];
end
dists=linspace(structure(1),structure(2),(npoints-4)/4+2);
for s=2:length(dists)-1
    basepoints(s+(npoints-4)/4+1,:)=[dists(s)*cos(angle(1)),dists(s)*sin(angle(1))];
    basepoints(npoints+2-s,:)=[dists(s)*cos(angle(2)),dists(s)*sin(angle(2))];
end


%temp_points=basepoints(1:npoints,:);
for n=2:ncell
    for s=1:npoints
        temppoints=basepoints(s,:).*[(1^(n-1)),1];
        basepoints((n-1)*npoints+s,:)=temppoints*rotatemat^(n-1);
    end
end











% for s=1:npoints
%     basepoints((ncell-1)*npoints+s,:)=basepoints(s,:)*-rotatemat;
% end

% for n=1:ncell
%     for s=1:npoints
%         bps((n-1)*npoints+s,:)=basepoints((n-1)*npoints+s,:)*rotatemat;
%     end
% end
% basepoints=bps;




if ncell==0
scatter(basepoints(:,1),basepoints(:,2))
plot(polyshape(basepoints(npoints+1:npoints*2,1),basepoints(npoints+1:npoints*2,2)))
plot(basepoints(:,1),basepoints(:,2))

for n=1:ncell
    hold on
    plot(polyshape(basepoints(((n-1)*npoints+1):(n-1)*npoints+npoints,1),basepoints(((n-1)*npoints+1):(n-1)*npoints+npoints,2)))
    
end


end

























end