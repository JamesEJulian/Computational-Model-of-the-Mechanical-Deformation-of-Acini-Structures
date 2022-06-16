function [basepoints]=circles(ncell,npoints,structrad,rad,isnuc,nucrat)
rotateang=2*pi/(ncell);
rotatemat=[cos(rotateang),-sin(rotateang);sin(rotateang),cos(rotateang)];

basepoints=NaN(npoints,2);    
    for s=1:npoints
        aa=s*2*pi/npoints+pi/6;
        basepoints(s,1)=rad*cos(aa);
        basepoints(s,2)=structrad+rad*sin(aa);        
    end

% for n=2:ncell
%     for s=1:npoints
%         basepoints((n-1)*npoints+s,:)=basepoints((n-2)*npoints+s,:)*rotatemat;
%     end
% end

for n=2:ncell
    for s=1:npoints
        temppoints=basepoints(s,:).*[(1^(n-1)),1];
        basepoints((n-1)*npoints+s,:)=temppoints*rotatemat^(n-1);
    end
end




if ncell==0
scatter(basepoints(:,1),basepoints(:,2))
%plot(polyshape(basepoints(:,1),basepoints(:,2)))
%plot(basepoints(:,1),basepoints(:,2))

for n=1:ncell
    hold on
    plot(polyshape(basepoints(((n-1)*npoints+1):(n-1)*npoints+npoints,1),basepoints(((n-1)*npoints+1):(n-1)*npoints+npoints,2)))
    
end
end

end