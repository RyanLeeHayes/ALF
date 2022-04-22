Nmc=5e6;
Nl=400;
c=5.5;
dl=1/Nl;

for Ndim=2:10

    Ndim

    histogram1=zeros([Nl,1]);

    for i=1:Nmc
        ti=rand([1,Ndim]);
        unli=exp(c*sin(pi*ti-pi/2));
        li=unli/sum(unli);
        indi=floor(li/dl)+1;
        for j=1:Ndim
            histogram1(indi(j))=histogram1(indi(j))+1;
        end
    end

    S1_k=log(histogram1);
    G1_k=-S1_k+S1_k(floor(Nl/2));

    % plot(dl/2:dl:1,histogram1)
    % plot(dl/2:dl:1,G1_k')

    save(['G1_',num2str(Ndim),'.dat'],'G1_k','-ascii')

end
