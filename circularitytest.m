function circ=circularitytest(x00,ncell,npoints,attached)
circ=NaN(7,ncell);
dist=zeros(npoints,1);
%circ   1=cell area
%       2=cell perimeter
%       3=cell perimeter/area
%       4=average distance from center of cell
%       5=structure area
%       6=structure perimeter
%       7=structure perimeter/area


for n=1:ncell
    ppnts(:,n,1)=x00(((n-1)*npoints+1):((n)*npoints));
    ppnts(:,n,2)=x00(ncell*npoints+(n-1)*npoints+1:ncell*npoints +(n)*npoints);
    circ(1,n)=polyarea(ppnts(:,n,1),ppnts(:,n,2));
    circ(2,n)=perimeter(polyshape(ppnts(:,n,1),ppnts(:,n,2)));
    circ(3,n)=(4*pi*circ(1,n))/(circ(2,n))^2;
    cent=[mean(ppnts(:,n,1)),mean(ppnts(:,n,2))];
    for s=1:npoints
        dist(s)=sqrt((ppnts(s,n,1)-cent(1))^2+(ppnts(s,n,2)-cent(2))^2);
    end
    circ(4,n)=mean(dist);
end

center(1)=mean(x00(1:ncell*npoints));
center(2)=mean((x00((ncell*npoints+1):(ncell*npoints*2))));
center;

%scatter(center(1),center(2));
exbound=[0,0,0];
intit=1;

ii=1;
for n=1:ncell
   dd=sqrt((mean(ppnts(:,n,1))-center(1))^2+(mean(ppnts(:,n,2))-center(2))^2);
   for s=npoints:-1:1
       d=sqrt((ppnts(s,n,1)-center(1))^2+(ppnts(s,n,2)-center(2))^2);
       %attached(ii)
       if d>=dd && (attached(((n-1)*npoints+s))==0)
           if(intit==1)
            intit=intit+1;
            exbound(1)=ppnts(s,n,1);
            exbound(2)=ppnts(s,n,2);
            exbound(3)=(n-1)*npoints+s;
           else
               exbound=[exbound;ppnts(s,n,1),ppnts(s,n,2),(n-1)*npoints+s];
           end
       end
       ii=ii+1;
       
   end
   c=mean(exbound,1);
    d=exbound-c;
    th=atan2(d(:,2),d(:,1));
    [th, idx] = sort(th);   % sorting the angles 
    exbound=exbound(idx,:);
end

master=polyshape(ppnts(:,1,1),ppnts(:,1,2));
for n=2:ncell
   master=union(master,polyshape(ppnts(:,n,1),ppnts(:,n,2)));
end
%figure;
%plot(master);
hold on;
isbad=0;
for n=1:length(exbound(:,1))

    %lin=[center(1),intbound(n,1);center(2),intbound(n,2)];
    lineseg=[center(1),center(2);exbound(n,1),exbound(n,2)];
    %plot(lineseg(:,1),linseg(:,2));
    [in,out] = intersect(master,lineseg);
    hold on;
   % plot(in(:,1),in(:,2),'b',out(:,1),out(:,2),'r')
    %figure;
    if isempty(in)~=0
    isbad=isbad+1;
    exbound(n,:)=[0,0,0];
    end
end
%figure
n=1;
while n<=length(exbound(:,1))
   if exbound(n,1)==0&&exbound(n,2)==0&&exbound(n,3)==0
       exbound(n,:)=[];
       %disp('hi');
   else
       n=n+1;
   end
end


circ(5,1)=polyarea(exbound(:,1),exbound(:,2));
circ(6,1)=perimeter(polyshape(exbound(:,1),exbound(:,2)));
circ(7,1)=(4*pi*circ(5,1))/(circ(6,1))^2;





% % hold on
% % plot(master);
% % for n=1:length(exbound(:,1))
% %     lineseg=[center(1),center(2);exbound(n,1),exbound(n,2)];
% %     plot(lineseg(:,1),lineseg(:,2));
% %     
% % end
% % 
% % isint=zeros(ncell*npoints*4,1);
% % for n=1:length(exbound(:,1))
% %     isint(exbound(n,3))=exbound(n,1);
% %     isint(exbound(n,3)+ncell*npoints)=exbound(n,2);
% % end
% % 
% % figure;
% % hold on;
% % plot(master);
% % for s=1:ncell*npoints
% %     if isint(s)~=0
% %     scatter(x00(s),x00(s+ncell*npoints),'*');
% %     scatter(isint(s),isint(s+ncell*npoints),'o');
% %     end
% % end
%plot(polyshape(intbound(:,1),intbound(:,2)));

end
