function xx=addnewcell(xi,newpoints,ncell,npoints)
%%% adds a new cell to the matrix NOT IN USE   
global vin
global v
celltemp=zeros(4*(ncell+1)*npoints);
    start=1;
    ending=(ncell)*npoints;
    for n=1:4
       celltemp(start:ending)=xi(((n-1)*ncell*npoints+1):((n)*ncell*npoints));
       celltemp((ending+1):(ending+npoints))=newpoints(:,n);
       start=(ncell+1)*npoints*n+1;
       ending=start+ncell*npoints-1;
    end
    xx=celltemp;
    v(1:ncell*npoints*2)=vin;
    v((ncell*npoints*2+1):ncell*npoints*4)=vin*1.5;
end