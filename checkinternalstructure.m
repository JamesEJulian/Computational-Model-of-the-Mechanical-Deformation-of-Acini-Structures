function isinside=checkinternalstructure(x00,ncell,npoints)
%%% is commented out b/c not currently used.
% clear isinside
% clear int
% global F_external
% global k_cell
% global k_cellv
% global k_cout
% global k_coutbase
% global intspring
% global k_in
% global v
% global vin
% global 
% isinside=ones(ncell,4); % builds output matrix [xposition, yposition, xdireciton, ydirection]
% isinside(:,3:4)=0;      %zeros direction
% cellcenters=zeros(ncell,2);
% midpoint=zeros(ncell,2);
% %%% finds centers and makes easier to use matrix for points
% for n=1:ncell
%     ppnts(:,n,1)=x00(((n-1)*npoints+1):((n)*npoints));
%     ppnts(:,n,2)=x00(ncell*npoints+(n-1)*npoints+1:ncell*npoints +(n)*npoints);
%     cellcenters(n,:)=[mean(ppnts(:,n,1)),mean(ppnts(:,n,2))];
% end
% %%% defines points as corners of a box
% centerz(:,1)=[1,1,0,-1,-1,-1,0,1];
% centerz(:,2)=[0,1,1,1,0,-1,-1,-1];
% centerz=centerz*40; %scales box size up 
% 
% master=polyshape(ppnts(:,1,1),ppnts(:,1,2)); %builds ployshape of the overall structure
% for n=2:ncell
%    master=union(master,polyshape(ppnts(:,n,1),ppnts(:,n,2)));
% end
% 
% %%% checks if each cell point can see at least one of the corners
% for n=1:ncell
%     isbad=0;
%     int=[];
%     celltest=0;
%     for s=1:npoints
%         for ss=1:8
%         lineseg=[centerz(ss,1),centerz(ss,2);ppnts(s,n,1),ppnts(s,n,2)];
%         [in,out] = intersect(master,lineseg);
% 
%             if isempty(in)==1 && celltest>=8
%                 isinside(n,:)=0;
%             elseif isempty(in)==1 && celltest<88
%                 celltest=celltest+1;
%             end
%         end
%         if isinside(n)==0
%             break
%         end
%     end
%         if isinside(n)==0
%             continue
%         end
% end
% 
%      %%% could be used to scale down spring stiffnesses of internal cell
%      %%% and its neighbors
% % k_cell(:)=k_cellv;
% % k_cout(:)=k_coutbase;
% % k_in(:)=intspring;
% % v(:)=vin;
% %%% finds the closest cells to the inside cells and makes the force vector
% %%% that pulls the inside cell out of the structure
% for n=1:ncell
%    if isinside(n)==1
%      %% could be used to scale down spring stiffnesses of internal cell
%      %% and its neighbors
%       F_external(((n-1)*npoints+1):(n)*npoints)=0;
%       k_cout(((n-1)*npoints+1):(n)*npoints)=k_coutbase*0.8;
%       k_in(((n-1)*npoints+1):(n)*npoints)=intspring*0.8;
%      k_cell(((n-1)*npoints+1):(n)*npoints)=k_cellv/2;
%      v(((n-1)*npoints+1):(n)*npoints)=vin*0.5;
%       lengths=zeros(npoints,8);
%       shortest=zeros(npoints,4);
%             for ss=1:ncell
%                 if ss~=n
%                 lengths(n,ss)=sqrt((cellcenters(ss,1)-cellcenters(n,1))^2+((cellcenters(ss,2)-cellcenters(n,2)))^2);
%                 else 
%                 lengths(n,ss)=1000000;
%                 end
%             end
%         [a b]=sort(lengths(n,:));
%         shortest(n,:)=[a(1),b(1),a(2),b(2)];
% 
%         isinside(n,1)=b(1);
%         isinside(n,2)=b(2);
%         
%      %%% finds the midpoint between the two closest cells to the internal
%      %%% cell
%      midpoint(n,:)=[(cellcenters(isinside(n,1),1)-cellcenters(isinside(n,2),1))/2,(cellcenters(isinside(n,1),2)-cellcenters(isinside(n,2),2))/2];
%      isinside(n,3:4)=-1*[midpoint(n,2),midpoint(n,1)];
%      tempdist=sqrt(midpoint(n,2)^2+midpoint(n,1)^2);
%      isinside(n,3:4)=isinside(n,3:4)/tempdist;
%     
%      midpoint(n,:)=midpoint(n,:)+[cellcenters(isinside(n,2),1),cellcenters(isinside(n,2),2) ];
% 
%      point2mid=ones(2,1)*10000;
%      %%% sets the direction to the midpoint from the internal cell
%      for s=1:npoints
%          ddist=sqrt((midpoint(n,1)-ppnts(s,n,1))^2+(midpoint(n,2)-ppnts(s,n,2))^2);
%          if ddist<point2mid(1)
%             point2mid(:)=[ddist,s]; 
%          end
%      end
%      isinside(n,1:2)=[point2mid(2),0];
% 
%      %%% could be used to scale down spring stiffnesses of internal cell
%      %%% and its neighbors
%       k_cell(((b(1)-1)*npoints+1):(b(1))*npoints)=k_cellv/10;
%       k_cell(((b(2)-1)*npoints+1):(b(2))*npoints)=k_cellv/10;
%       k_cout(((b(1)-1)*npoints+1):(b(1))*npoints)=k_coutbase/6;
%       k_cout(((b(2)-1)*npoints+1):(b(2))*npoints)=k_coutbase/6;
%       k_in(((b(1)-1)*npoints+1):(b(1))*npoints)=intspring/10;
%       k_in(((b(2)-1)*npoints+1):(b(2))*npoints)=intspring/10;
%    end
    isinside=zeros(ncell,4); % builds output matrix [xposition, yposition, xdireciton, ydirection]
end
