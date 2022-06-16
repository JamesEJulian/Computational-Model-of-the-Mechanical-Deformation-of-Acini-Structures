function cps = pointgen(center, npts, rad)
%% Generates Circle Points
cps=zeros(npts,2);    
    for s=1:npts
        aa=s*2*pi/npts+pi/2;
        cps(s,1)=center(1)+rad*cos(aa);
        cps(s,2)=center(2)+rad*sin(aa);
        
    end
end
