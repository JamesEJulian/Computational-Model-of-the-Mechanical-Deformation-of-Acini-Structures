function neworder=fixorder(ncell,npoints,x0,n)
    %Reorders Points if they switch places relative to their index 
    neworder=x0;
    ppnts(:,1)=x0(((n-1)*npoints+1):((n)*npoints));
    ppnts(:,2)=x0(ncell*npoints+(n-1)*npoints+1:ncell*npoints +(n)*npoints);
    nppnts(:,1)=x0((2*ncell*npoints+((n-1)*npoints+1)):2*ncell*npoints+((n)*npoints));
    nppnts(:,2)=x0((3*ncell*npoints+((n-1)*npoints+1)):3*ncell*npoints+((n)*npoints));
    
    
    c=mean(ppnts,1);
    d=ppnts-c;
    th=atan2(d(:,2),d(:,1));
    [th, idx] = sort(th);   % sorting the angles 
    ppntss=ppnts(idx,:);
    
    c=mean(nppnts,1);
    d=nppnts-c;
    th=atan2(d(:,2),d(:,1));
    [th, idx] = sort(th);   % sorting the angles 
    nppntss=nppnts(idx,:);
    
    
    neworder(((n-1)*npoints+1):((n)*npoints))=ppntss(:,1);
    neworder(ncell*npoints+(n-1)*npoints+1:ncell*npoints +(n)*npoints)=ppntss(:,2);
    neworder((2*ncell*npoints+((n-1)*npoints+1)):2*ncell*npoints+((n)*npoints))=nppntss(:,1);
    neworder((3*ncell*npoints+((n-1)*npoints+1)):3*ncell*npoints+((n)*npoints))=nppntss(:,2);
    
    if n<0
         ppntss(:,1)= neworder(((n-1)*npoints+1):((n)*npoints));
         ppntss(:,2)= neworder(ncell*npoints+(n-1)*npoints+1:ncell*npoints +(n)*npoints);
  figure
hold on
    plot(polyshape(ppntss(:,1),ppntss(:,2)));
    %plot(polyshape(ppnts(:,1),ppnts(:,2)));
    end
end