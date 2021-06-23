function  [sWr,sWd]= FinderDR(W)
sWr=zeros(size(W,2));
sWd=zeros(size(W,1));
for i=1:size(W,2)
    for k=1:size(W,2)
        s=0;
        for j=1:size(W,1)
            if W(j,i)==1
                if W(j,k)==W(j,i)
                 s=s+1;
                end
            end
       end
        sWr(i,k)=s;
        sWr(k,i)=s;
    end
end
for i=1:size(W,1)
    for k=1:size(W,1)
        n=0;
        for j=1:size(W,2)
            if W(i,j)==1
                if W(k,j)==W(i,j)
                 n=n+1;
                end
            end
       end
        sWd(i,k)=n;
        sWd(k,i)=n;
    end
end  
    
end

