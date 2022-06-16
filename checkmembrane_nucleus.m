   close all;
    figure1=figure;
 
    newcolors={'#A2142F','#FF0000','#D95319','#EDB120', '#FFFF00','#00FF00','#77AC30','#0072BD','#0000FF','#0072BD'};
    %newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
    colororder(newcolors)
    n=1;
    if rem(total_t,30)==0
        hasfixed=0;
    else
        hasfixed=999;
    end
%     x0t=x0;
while n<=ncell
    if n==1
    colororder(newcolors)
    end
    lastwarn('');
    ppnts(:,n,1)=x0(((n-1)*npoints+1):((n)*npoints));
    ppnts(:,n,2)=x0(ncell*npoints+(n-1)*npoints+1:ncell*npoints +(n)*npoints);
    plot(polyshape(ppnts(:,n,1),ppnts(:,n,2)));

    hold on;
    if ~isempty(lastwarn) && hasfixed<ncell*2
        x0=fixorder(ncell,npoints,x0,n);
            close all;
        figure1=figure;
        %newcolors={'#A2142F','#FF0000','#D95319','#EDB120', '#FFFF00','#00FF00','#77AC30','#0072BD','#0000FF','#0072BD'};
        %colororder(newcolors)
%         if sum(isnan(x0t))~=0
%             x0=x0t;
         n=0;
         hasfixed=hasfixed+1;
         close;
         figure;
%         else
%             disp('hi');
%         end
    end
    %scatter(ppnts(1,n,1),ppnts(1,n,2),'*');        %can be used to show
    %each membrane p
    if n>=1
    ppnts(:,n,1)=x0((2*ncell*npoints+((n-1)*npoints+1)):2*ncell*npoints+((n)*npoints));
    ppnts(:,n,2)=x0((3*ncell*npoints+((n-1)*npoints+1)):3*ncell*npoints+((n)*npoints));
    plot(polyshape(ppnts(:,n,1),ppnts(:,n,2)));
    end
    n=n+1;
end

for s=1:6
   plot([x0(s),x0(s+ncell*npoints*2)],[x0(s+ncell*npoints),x0(s+ncell*npoints*3)]); 
   scatter(x0(s),x0(s+ncell*npoints));
   scatter(x0(s+ncell*npoints*2),x0(s+ncell*npoints*3));
end