% DIRs={'94a','94b','94c','94d','94e'};
DIRs={'94a','94b','94c'};
NREPS=load('../nreps');

kT=0.001987*298;
dl=0.01;
Eedge=reshape(0:dl:1,[],1);
Emid=Eedge(1:(end-1))+dl/2;

nblocks=load('../nblocks');
nsubs=load('../nsubs');

for i=1:length(nsubs)
    G1{i}=reshape(load(['G_imp/G1_',num2str(nsubs(i)),'.dat']),[],1);
    G12{i}=reshape(load(['G_imp/G12_',num2str(nsubs(i)),'.dat']),[],1);
end

for i=1:length(DIRs)
  DIR=DIRs{i}
  for j=1:NREPS
    L{i,j}=load([DIR,'/Lambda.bp',num2str(j-1),'.dat']);

ibuff=0;
for ii=1:length(nsubs)
    block=(1:nsubs(ii))+ibuff;
    for jj=block
        % fnm=['G_',num2str(jj),'.dat'];
        % tmp=load(fnm);
        % Result1{jj}=tmp-kT*G1{ii};

        tmp=histc(L{i,j}(:,jj+1),Eedge);
        tmp(end-1)=tmp(end-1)+tmp(end);
        tmp(end)=[];
        tmp=reshape(tmp,[],1);
        tmp=-kT*log(tmp);
        tmp=tmp-kT*G1{ii};

        % mask=(isfinite(tmp) & isfinite(Result1{j}));
        % meanR=mean(Result1{j}(mask));
        % meanT=mean(tmp(mask));
        % tmp=tmp-(meanT-meanR);

        Test1{i,j,jj}=tmp;
        dG(i,j,jj)=tmp(end)-tmp(1);
    end
    ibuff=ibuff+nsubs(ii);
end



  end
end


