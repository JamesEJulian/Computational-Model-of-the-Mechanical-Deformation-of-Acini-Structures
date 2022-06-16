function resarea=getarea(xi,ncell,npoints)
%Finds the area of each cell excluding the nucleus
resarea=zeros(ncell);
for n=1:ncell
    pntsc(:,n,1)=xi(((n-1)*npoints+1):((n)*npoints));
    pntsc(:,n,2)=xi(ncell*npoints+(n-1)*npoints+1:ncell*npoints +(n)*npoints);
    areac(n)=polyarea(pntsc(:,n,1),pntsc(:,n,2));  %find area of cell
    pntsn(:,n,1)=xi((2*ncell*npoints+((n-1)*npoints+1)):2*ncell*npoints+((n)*npoints));
    pntsn(:,n,2)=xi((3*ncell*npoints+((n-1)*npoints+1)):3*ncell*npoints+((n)*npoints));
    arean(n)=polyarea(pntsn(:,n,1),pntsn(:,n,2));  %find area of nucleus
end
    resarea=areac-arean;
end