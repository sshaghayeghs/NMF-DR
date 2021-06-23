clc 
clear
DiDrAMat=importdata('DiDrA.txt');
k=0;
%%
while (k~=sum(sum(DiDrAMat,1),2))
    for i=1:size(DiDrAMat,1)
        for j=1:size(DiDrAMat,2)
            if (DiDrAMat(i,j)==1)
                    k=k+1;
                 KA(k,1)=DiDrAMat(i,j);
                 KA(k,2)=i;
                 KA(k,3)=j;
            end
        end
    end
end
%%
c=cvpartition(2532,'KFold',10);
%%
Testcase1=KA(test(c,1),:);Testcase2=KA(test(c,2),:);Testcase3=KA(test(c,3),:);Testcase4=KA(test(c,4),:);Testcase5=KA(test(c,5),:);
Testcase6=KA(test(c,6),:);Testcase7=KA(test(c,7),:);Testcase8=KA(test(c,8),:);Testcase9=KA(test(c,9),:);Testcase10=KA(test(c,10),:);
%%
DiDrAMatcase1=DiDrAMat; DiDrAMatcase2=DiDrAMat;DiDrAMatcase3=DiDrAMat;DiDrAMatcase4=DiDrAMat;DiDrAMatcase5=DiDrAMat;
DiDrAMatcase6=DiDrAMat;DiDrAMatcase7=DiDrAMat;DiDrAMatcase8=DiDrAMat;DiDrAMatcase9=DiDrAMat;DiDrAMatcase10=DiDrAMat;
%%
for i=1:size(Testcase1,1)
   DiDrAMatcase1(Testcase1(i,2),Testcase1(i,3))=0;
end
for i=1:size(Testcase2,1)
   DiDrAMatcase2(Testcase2(i,2),Testcase2(i,3))=0;
end
for i=1:size(Testcase3,1)
   DiDrAMatcase3(Testcase3(i,2),Testcase3(i,3))=0;
end
for i=1:size(Testcase4,1)
   DiDrAMatcase4(Testcase4(i,2),Testcase4(i,3))=0;
end
for i=1:size(Testcase5,1)
   DiDrAMatcase5(Testcase5(i,2),Testcase5(i,3))=0;
end
for i=1:size(Testcase6,1)
   DiDrAMatcase6(Testcase6(i,2),Testcase6(i,3))=0;
end
for i=1:size(Testcase7,1)
   DiDrAMatcase7(Testcase7(i,2),Testcase7(i,3))=0;
end
for i=1:size(Testcase8,1)
   DiDrAMatcase8(Testcase8(i,2),Testcase8(i,3))=0;
end
for i=1:size(Testcase9,1)
   DiDrAMatcase9(Testcase9(i,2),Testcase9(i,3))=0;
end
for i=1:size(Testcase10,1)
   DiDrAMatcase10(Testcase10(i,2),Testcase10(i,3))=0;
end
