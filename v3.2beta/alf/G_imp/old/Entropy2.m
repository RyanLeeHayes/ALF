Nmc=5e6;
Nl=20;
c=5.5;
cutlsum=0.8;
dl=1/Nl;

for Ndim=2:10

    Ndim

    histogram2=zeros([Nl,Nl]);
    offdiag=ones([Ndim,Ndim])-eye(Ndim);

    for i=1:Nmc
        ti=rand([1,Ndim]);
        unli=exp(c*sin(pi*ti-pi/2));
        li=unli/sum(unli);
        indall=floor(li*Nl)+1;
        for j=1:Ndim
            for k=(j+1):Ndim
                histogram2(indall(j),indall(k))=histogram2(indall(j),indall(k))+1;
            end
        end
    end

    S2_k=log(histogram2);
    G2_k=-S2_k+mean(S2_k(isfinite(S2_k)));

    % plot(dl/2:dl:1,histogram1)
    % plot((dl/2):dl:1,G12_k')

    save(['G2_',num2str(Ndim),'.dat'],'G2_k','-ascii')

end
