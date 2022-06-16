function x0=fixOverlap(x0,ncell,nppts,l_bond, order)
%nppts=npoints;
%Forces Cells out of eachother if they start to overlap

totalind=ncell*nppts;
if order==0
for n=1:ncell
        ind(:,n,1)=(((n-1)*nppts+1):((n)*nppts));
        ind(:,n,2)=(totalind+(n-1)*nppts+1:totalind +(n)*nppts);
end
else
    for n=ncell:-1:1
        ind(:,n,1)=(((n-1)*nppts+1):((n)*nppts));
        ind(:,n,2)=(totalind+(n-1)*nppts+1:totalind +(n)*nppts);
    end
end


for n=1:ncell
%figure;
   temporder=1:ncell;
   temporder(n)=NaN;
   temporder=rmmissing(temporder);
   sumofshape=polyshape(x0(ind(:,temporder(1),1)),x0(ind(:,temporder(1),2)));
%hold on
%scatter(x0(ind(:,temporder(1),1)),x0(ind(:,temporder(1),2)),'*')
   for nx=2:length(temporder)
       sumofshape=union(sumofshape,polyshape(x0(ind(:,temporder(nx),1)),x0(ind(:,temporder(nx),2))));
                  
   end
%plot(sumofshape);
%title(num2str(n));
   for s=1:nppts
      lineseg=[x0(ind(s,n,1)),x0(ind(s,n,2));x0(ind(s,n,1)+totalind*2),x0(ind(s,n,2)+totalind*2)];
      [in,out]=intersect(sumofshape,lineseg);
%plot(in(:,1),in(:,2),'b',out(:,1),out(:,2),'r')
%ncell;
       if isempty(in)==0&&isempty(out)==0
           try
       newdist=sqrt((out(2,1)-out(1,1))^2+(out(2,2)-out(1,2))^2);
           catch
              break 
           end
       newvec=-out(2,:)+out(1,:);
       newdir=newvec/newdist;
       newpoint=newdir*(newdist-l_bond);
       newpoint=out(2,:)+newpoint;
%scatter(newpoint(1),newpoint(2),'*')
       x0(ind(s,n,:))=newpoint(:);
       end
%ncell;
   end
%ncell;
end
end