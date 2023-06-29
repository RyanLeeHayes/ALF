Nmc=5e6;
Nl=400;
c=5.5;
cutlsum=0.8;
dl=1/Nl;

for Ndim=2:10

    Ndim

    histogram12=zeros([Nl,1]);
    offdiag=ones([Ndim,Ndim])-eye(Ndim);

    for i=1:Nmc
        ti=rand([1,Ndim]);
        unli=exp(c*sin(pi*ti-pi/2));
        li=unli/sum(unli);
        lia=repmat(li,[Ndim,1]);
        lisum=lia+lia';
        l12a=find(offdiag & (lisum>cutlsum));
        for j=1:numel(l12a)
            l12=l12a(j);
            ind12=floor((lia(l12)/lisum(l12))/dl)+1;
            histogram12(ind12)=histogram12(ind12)+1;
        end
    end

    S12_k=log(histogram12);
    G12_k=-S12_k+S12_k(floor(Nl/2));

    % plot(dl/2:dl:1,histogram1)
    % plot((dl/2):dl:1,G12_k')

    save(['G12_',num2str(Ndim),'.dat'],'G12_k','-ascii')

end
