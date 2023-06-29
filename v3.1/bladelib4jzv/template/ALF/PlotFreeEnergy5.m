
Emid=(1/800):(1/400):1;
Emid2=(1/40):(1/20):1;

nsubs=load('../nsubs');
nblocks=load('../nblocks');

iG=1;
iF=1;

iblock=0;
for isite=1:length(nsubs)
jblock=iblock;
for jsite=isite:length(nsubs)
if isite==jsite

LEG={};
figure(iF)
iF=iF+1;
hold off
for i=1:nsubs(isite)
  LEG{i}=num2str(i);
  G1{i+iblock}=load(['multisite/G',num2str(iG),'.dat']);
  iG=iG+1;
  subplot(1,1,1)
  plot(Emid,G1{i+iblock})
  hold on
end
hold off
legend('location','best',LEG)

figure(iF)
iF=iF+1;
for i=1:nsubs(isite)
  for j=(i+1):nsubs(isite)
    G12{i+iblock,j+iblock}=load(['multisite/G',num2str(iG),'.dat']);
    iG=iG+1;
    subplot(nsubs(isite)-1,nsubs(isite)-1,(i-1)*(nsubs(isite)-1)+(j-1))
    hold off
    plot(Emid,G12{i+iblock,j+iblock})
  end
end

if nsubs(isite)>2
  figure(iF)
  iF=iF+1;
  for i=1:nsubs(isite)
    for j=(i+1):nsubs(isite)
      G2{i+iblock,j+iblock}=reshape(load(['multisite/G',num2str(iG),'.dat']),[20,20]);
      iG=iG+1;
      subplot(nsubs(isite)-1,nsubs(isite)-1,(i-1)*(nsubs(isite)-1)+(j-1))
      hold off
      surf(Emid2,Emid2,G2{i+iblock,j+iblock})
    end
  end
end

else

figure(iF)
iF=iF+1;
for i=1:nsubs(isite)
  for j=1:nsubs(jsite)
    G11{i+iblock,j+jblock}=reshape(load(['multisite/G',num2str(iG),'.dat']),[20,20]);
    iG=iG+1;
    subplot(nsubs(isite),nsubs(jsite),(i-1)*(nsubs(jsite))+(j))
    hold off
    surf(Emid2,Emid2,G11{i+iblock,j+jblock})
  end
end

end
jblock=jblock+nsubs(jsite);
end
iblock=iblock+nsubs(isite);
end
