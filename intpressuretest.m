function out=intpressuretest(x00)
%% Finds lumen Pressure 
global ncell
global npoints
global intbound
global isint
global attached
global interval
global c
global iswhole
%global center
%center=NaN;
for n=1:ncell
    ppnts(:,n,1)=x00(((n-1)*npoints+1):((n)*npoints));
    ppnts(:,n,2)=x00(ncell*npoints+(n-1)*npoints+1:ncell*npoints +(n)*npoints);
    %plot(polyshape(ppnts(:,n,1),ppnts(:,n,2)));
    hold on;
end

center(1)=mean(x00(1:ncell*npoints));
center(2)=mean((x00((ncell*npoints+1):(ncell*npoints*2))));
center;

iswhole=1;
if interval<3
for n=1:ncell
    attachedpoints=attached(((n-1)*npoints+1):((n)*npoints));
    checks=[];
    for s=1:length(attachedpoints)
        if attachedpoints(s)~=0
            if length(checks)==0
                checks=[checks c(s)];
            elseif checks(length(checks))~=c(s)
                checks=[checks c(s)];
            end
        end
    end
    if length(checks)<2
        iswhole=0;
       break;
    end

end
end
%scatter(center(1),center(2));
intbound=[0,0,0];
intit=1;





ii=1;
for n=1:ncell
   dd=sqrt((mean(ppnts(:,n,1))-center(1))^2+(mean(ppnts(:,n,2))-center(2))^2)*1.25;
   for s=npoints:-1:1
       d=sqrt((ppnts(s,n,1)-center(1))^2+(ppnts(s,n,2)-center(2))^2);
       %attached(ii)
       if d<=dd && (attached(((n-1)*npoints+s))==0)
           if(intit==1)
            intit=intit+1;
            intbound(1)=ppnts(s,n,1);
            intbound(2)=ppnts(s,n,2);
            intbound(3)=(n-1)*npoints+s;
           else
               intbound=[intbound;ppnts(s,n,1),ppnts(s,n,2),(n-1)*npoints+s];
           end
       end
       ii=ii+1;
       
   end
   cd=mean(intbound,1);
    d=intbound-cd;
    th=atan2(d(:,2),d(:,1));
    [th, idx] = sort(th);   % sorting the angles 
    intbound=intbound(idx,:);
end

master=polyshape(ppnts(:,1,1),ppnts(:,1,2));
for n=2:ncell
   master=union(master,polyshape(ppnts(:,n,1),ppnts(:,n,2)));
end
%figure;
%plot(master);
hold on;
isbad=0;
for n=1:length(intbound(:,1))

    %lin=[center(1),intbound(n,1);center(2),intbound(n,2)];
    lineseg=[center(1),center(2);intbound(n,1),intbound(n,2)];
    %plot(lineseg(:,1),linseg(:,2));
    [in,out] = intersect(master,lineseg);
    hold on;
   % plot(in(:,1),in(:,2),'b',out(:,1),out(:,2),'r')
    %figure;
    if isempty(in)==0
    isbad=isbad+1;
    intbound(n,:)=[0,0,0];
    end
end
%figure
n=1;
while n<=length(intbound(:,1))
   if intbound(n,1)==0&&intbound(n,2)==0&&intbound(n,3)==0
       intbound(n,:)=[];
       %disp('hi');
   else
       n=n+1;
   end
end
%hold on
%plot(master);
% for n=1:length(intbound(:,1))
%     lineseg=[center(1),center(2);intbound(n,1),intbound(n,2)];
%     %plot(lineseg(:,1),lineseg(:,2));
%     
% end

isint=zeros(ncell*npoints*4,1);
for n=1:length(intbound(:,1))
    isint(intbound(n,3))=intbound(n,1);
    isint(intbound(n,3)+ncell*npoints)=intbound(n,2);
end

% figure;
% hold on;
% plot(master);
% for s=1:ncell*npoints
%     if isint(s)~=0
%     scatter(x00(s),x00(s+ncell*npoints),'*');
%     scatter(isint(s),isint(s+ncell*npoints),'o');
%     end
% end
% plot(polyshape(intbound(:,1),intbound(:,2)));
out=polyarea(intbound(:,1),intbound(:,2));
end
