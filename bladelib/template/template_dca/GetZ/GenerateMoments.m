
h=load('h.dat');
J=load('J.dat');

Seq=zeros([1,10]);
nsubs=10*ones([1,10]);

Z=0;
m1=zeros([10,max(nsubs)]);
m2=zeros([10,max(nsubs),10,max(nsubs)]);
while Seq(end)<nsubs(end)
  SeqInd=Seq+1;
  for i=1:10
    SeqInd((i+1):end)=SeqInd((i+1):end)+nsubs(i);
  end
  w=exp(sum(h(SeqInd))+sum(sum(J(SeqInd,SeqInd))));
  Z=Z+w;
  for i=1:10
    m1(i,Seq(i)+1)=m1(i,Seq(i)+1)+w;
    for j=1:10
      m2(i,Seq(i)+1,j,Seq(j)+1)=m2(i,Seq(i)+1,j,Seq(j)+1)+w;
    end
  end
  Seq(1)=Seq(1)+1;
  for i=1:9
    if Seq(i)==nsubs(i)
      Seq(i)=0;
      Seq(i+1)=Seq(i+1)+1;
    end
  end
end


