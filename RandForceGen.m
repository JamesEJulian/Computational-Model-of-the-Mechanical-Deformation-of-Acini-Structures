function force_matrix = RandForceGen(ECM, point_set)
%%% Calculates Random ECM Forces
global ncell
global npoints
force_matrix = zeros(ncell*npoints*3,1);
  for p=1:ncell*3                 
    pp=randi(npoints);    %Random Place Number Generated
    force_matrix(((p-1)*npoints)+pp) = normrnd(ECM(p,1),ECM(p,2));    %Force Generation Value
  end

  %force_matrix = zeros(ncell*npoints,1);       %if needed zeros forces for
  %troubleshooting
end