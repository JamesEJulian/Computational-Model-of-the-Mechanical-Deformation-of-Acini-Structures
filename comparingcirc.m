clear all;
files=3;
ff=[string('Solver_ 1 86400 05-Oct-2021.mat') string('Solver_ 2 86400 05-Oct-2021.mat') string('Solver_ 3 86400 06-Oct-2021.mat')];
for n=1:files
load(ff(n),'circularity');
load(ff(n),'s1');
title=strjoin([string('circ'),string(n)]);
title2=strjoin([string('brokenst'),string(n)]);
varname=genvarname(title);
varname2=genvarname(title2);
assignin('base',varname,circularity);
assignin('base',varname2,s1);
end

hold on;
%%% Comparing structure circularity (perimeter/area)
for n=1:files
    t1=squeeze(circ1(7,1,:));
    plot(t1,'b');
    t2=squeeze(circ2(7,1,:));
    plot(t2,'r');
    t3=squeeze(circ3(7,1,:));
    plot(t3,'g');
end
%%% comparing broken cell to its nearest nonbroken neighbor, may not always
%%% work
figure; 
hold on;
c1=squeeze(circ1(3,brokenst1,:));
plot(c1,'b');
c2=squeeze(circ1(3,12,:));
plot(c2,'r');



solocomp=circ3;

brokenmat=NaN(length(squeeze(solocomp(3,1,:))));
normat=NaN(length(squeeze(solocomp(3,1,:))));
for s=1:length(squeeze(solocomp(3,1,:)))
    brok=0;
    nor=0;
    for n=1:length(squeeze(solocomp(3,:,1)))
        if n>=brokenst1 && n<brokenst1+5
            brok=solocomp(3,n,s)+brok;
        else
             nor=solocomp(3,n,s)+nor;
        end
    end
    brokenmat(s)=brok/5;
    normat(s)=nor/7;
end
figure;
hold on;
plot(brokenmat,'b');
plot(normat,'r');