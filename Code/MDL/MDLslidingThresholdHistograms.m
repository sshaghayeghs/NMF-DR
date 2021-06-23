function [Ldist,Lvals]=MDLslidingThresholdHistograms(V,precision,Wstore,Hstore,rValues,options)
nRuns=length(Wstore);
V=V/max(max(V));%%%%
[m,n]=size(V);
Lvals=zeros(nRuns,6);Ldist=zeros(nRuns,6);
plotNum=1;
for i=1:nRuns
    W=Wstore{i,1};H=Hstore{i,1};
    zerosNum=10;zeroBounds=linspace(0,precision,zerosNum);
    r=rValues(i);
    LwAll=zeros(zerosNum,zerosNum);LhAll=zeros(zerosNum,zerosNum);
    Lw0All=zeros(zerosNum,zerosNum);Lh0All=zeros(zerosNum,zerosNum);
    LtotAll=zeros(zerosNum,zerosNum);LeAll=zeros(zerosNum,zerosNum);
    descLengths=cell(3,2);
    minLwh=inf;
    
    %%%%%%%%%%%%%%%%%%% FOR THE W VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for kw=1:zerosNum
        W2=reshape(W,[],1);
        %%%%%%%%%%%%%%%%%% FIRST SORT OUT THE W0 TERMS %%%%%%%%%%%%%%%%%%%%
        zeroVal=zeroBounds(kw);
        W0sLoc=find(W2(:,1)<zeroVal);numW0s=length(W0sLoc);
        numW1s=length(W2)-numW0s;numWTotal=length(W2);
        W3=W2;W2(W0sLoc,:)=[];
        if(numW0s~=0)
            Lw0=-numW0s*log2(numW0s/numWTotal)-numW1s*log2(numW1s/numWTotal);
        else
            Lw0=0;
        end
        %%% NOW THE W+ TERMS
        [Wvals,Wlocs]=histcounts(W2,'binwidth',precision);
        WvalsProb=Wvals/(sum(Wvals));
        [WvalsCumu,WlocsCumu]=histcounts(W2,'binwidth',precision,'Normalization','cumcount');%%%% WvalsCumu CONTAINS comulative results
        WvalsCumu2=WvalsCumu/max(WvalsCumu);
        Wlocs2=[(Wlocs(1)-precision),Wlocs(1:end)+precision];
        Wbinned=interp1(Wlocs2,Wlocs2,W3,'nearest');
        Wbinned(W0sLoc)=0;Wbinned=reshape(Wbinned,m,r);
        [parW,~] = gamfit(Wlocs2(2:end-1),[],[],Wvals);
        yW=gampdf(Wlocs2(2:end-1),parW(1),parW(2));nBinsW=length(Wvals);
        WdistCumu=gamcdf(WlocsCumu(2:end),parW(1),parW(2));
        WdistProb=precision*yW; % The probability of being in each bin from the gamma distribution is calculated
        WdistProb3=zeros(length(Wvals),1);WvalsProb3=zeros(length(Wvals),1); % *dist* are for distributions, *vals* are for terms drawn directly from histograms
        WdistProb4=zeros(length(Wvals),1);WvalsProb4=zeros(length(Wvals),1); % *Prob3 contains the description length for each bin, *Prob4 the cumulative description length
        for j=1:nBinsW
            if(WdistProb(j)==0)
                WdistProb(j)=min(WdistProb(WdistProb>0)); % This sets zero probability terms to the smallest non-zero probability, this is an approximation.
            end
            if(WvalsProb(j)==0)
                WvalsProb(j)=min(WvalsProb(WvalsProb>0)); % For the histograms this is just to prevent taking the log of a zero.
            end
            WdistProb2=-log2(WdistProb(j));WvalsProb2=-log2(WvalsProb(j));
            WdistProb3(j,1)=Wvals(j)*WdistProb2;WvalsProb3(j,1)=Wvals(j)*WvalsProb2;
            if(j>1)
                WdistProb4(j,1)=WdistProb3(j,1)+WdistProb4(j-1,1);
                WvalsProb4(j,1)=WvalsProb3(j,1)+WvalsProb4(j-1,1);
            else
                WdistProb4(j,1)=WdistProb3(j,1);
                WvalsProb4(j,1)=WvalsProb3(j,1);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% END OF THE W VALUES CODE %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% FOR THE H VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for kh=1:zerosNum
            H2=reshape(H,[],1);
            %%%%%%%%%%%%%%%%%% FIRST SORT OUT THE H0 TERMS %%%%%%%%%%%%%%%%%%%%
            zeroValH=zeroBounds(kh);
            H0sLoc=find(H2(:,1)<zeroValH);numH0s=length(H0sLoc);
            numH1s=length(H2)-numH0s;numHTotal=length(H2);
            H3=H2;
            H2(H0sLoc,:)=[];
            if(numH0s~=0)
                Lh0=-numH0s*log2(numH0s/numHTotal)-numH1s*log2(numH1s/numHTotal);
            else
                Lh0=0;
            end
            %%% NOW THE H+ TERMS
            [Hvals,Hlocs]=histcounts(H2,'binwidth',precision);
            HvalsProb=Hvals/(sum(Hvals));
            [HvalsCumu,HlocsCumu]=histcounts(H2,'binwidth',precision,'Normalization','cumcount');%%%% WvalsCumu CONTAINS comulative results
            HvalsCumu2=HvalsCumu/max(HvalsCumu);
            Hlocs2=[(Hlocs(1)-precision),Hlocs(1:end)+precision];
            Hbinned=interp1(Hlocs2,Hlocs2,H3,'nearest');
            Hbinned(H0sLoc)=0;Hbinned=reshape(Hbinned,r,n);
            [parH,~] = gamfit(Hlocs2(2:end-1),[],[],Hvals);
            yH=gampdf(Hlocs2(2:end-1),parH(1),parH(2));nBinsH=length(Hvals);
            HdistCumu=gamcdf(HlocsCumu(2:end),parH(1),parH(2));
            HdistProb=precision*yH;
            HdistProb3=zeros(length(Hvals),1);HvalsProb3=zeros(length(Hvals),1);
            HdistProb4=zeros(length(Hvals),1);HvalsProb4=zeros(length(Hvals),1);
            for j=1:nBinsH
                if(HdistProb(j)==0)
                    HdistProb(j)=min(HdistProb(HdistProb>0));
                end
                if(HvalsProb(j)==0)
                    HvalsProb(j)=min(HvalsProb(HvalsProb>0));
                end
                HdistProb2=-log2(HdistProb(j));HvalsProb2=-log2(HvalsProb(j));
                HdistProb3(j,1)=Hvals(j)*HdistProb2;HvalsProb3(j,1)=Hvals(j)*HvalsProb2;
                if(j>1)
                    HdistProb4(j,1)=HdistProb3(j,1)+HdistProb4(j-1,1);
                    HvalsProb4(j,1)=HvalsProb3(j,1)+HvalsProb4(j-1,1);
                else
                    HdistProb4(j,1)=HdistProb3(j,1);
                    HvalsProb4(j,1)=HvalsProb3(j,1);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%% END OF THE H VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%% FOR THE E VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Vpred=Wbinned*Hbinned;
            [~,Vlocs]=histcounts(V,'binwidth',precision);
            if(Vlocs(1)~=0);Vlocs2=[0 Vlocs];else Vlocs2=Vlocs;end
            if(Vlocs(end)~=1);Vlocs2=[Vlocs2 1];end
            Vbinned=interp1(Vlocs2,Vlocs2,V,'nearest');
            E=Vbinned-Vpred;Ecol=reshape(E,1,[]);
            [Evals,Elocs]=histcounts(Ecol,'binwidth',precision);%%%% Evals CONTAINS THE NUMBER OF COUNTS IN EACH BIN, Elocs THE EDGES OF THE BINS
            Elocs2=Elocs(1:end-1)+precision;
            EvalsProb=Evals/(length(Ecol));
            [EvalsCumu,ElocsCumu]=histcounts(Ecol,'binwidth',precision,'Normalization','cumcount');%%%% EvalsCumu CONTAINS cumulative results
            EvalsCumu2=EvalsCumu/max(EvalsCumu);
            %%%% Elocs2 CONTAINS THE CENTRES OF THE BINS
            [muE,sigE,~,~] = normfit(Elocs2,[],[],Evals);% muE IS THE MEAN
            pEcumu = normcdf(Elocs(2:end),muE,sigE);
            yE = normpdf(Elocs2,muE,sigE);nBinsE=length(Evals);
            EdistProb=precision*yE;        % IF PROBABILITY IS NON-ZERO USE IT
            yEcumu=zeros(1,length(EdistProb));yEcumu(1)=EdistProb(1);
            for ii=2:length(yE)
                yEcumu(ii)=yEcumu(ii-1)+EdistProb(ii);
            end
            EdistProb3=zeros(nBinsE,1);EvalsProb3=zeros(nBinsE,1);
            EdistProb4=zeros(nBinsE,1);EvalsProb4=zeros(nBinsE,1);
            for j=1:nBinsE
                if(EdistProb(j)==0)
                    EdistProb(j)=min(EdistProb(EdistProb>0)); % IF PROB IS ZERO USE THE SMALLEST NON-ZERO VALUE TO PREVENT NAN PROBLEMS
                end
                if(EvalsProb(j)==0)
                    EvalsProb(j)=min(EvalsProb(EvalsProb>0)); % IF PROB IS ZERO USE THE SMALLEST NON-ZERO VALUE TO PREVENT NAN PROBLEMS
                end
                EdistProb2=-log2(EdistProb(j));EvalsProb2=-log2(EvalsProb(j));
                EdistProb3(j,1)=Evals(j)*EdistProb2;EvalsProb3(j,1)=Evals(j)*EvalsProb2;
                if(j>1)
                    EdistProb4(j,1)=EdistProb3(j,1)+EdistProb4(j-1,1);
                    EvalsProb4(j,1)=EvalsProb3(j,1)+EvalsProb4(j-1,1);
                else
                    EdistProb4(j,1)=EdistProb3(j,1);
                    EvalsProb4(j,1)=EvalsProb3(j,1);
                end
            end
            %%% For this am using the description lengths of the histograms
            LwAll(kh,kw)=WvalsProb4(end);LhAll(kh,kw)=HvalsProb4(end);
            Lw0All(kh,kw)=Lw0;Lh0All(kh,kw)=Lh0;
            LeAll(kh,kw)=EvalsProb4(end);
            LtotAll(kh,kw)=LwAll(kh,kw)+LhAll(kh,kw)+Lw0All(kh,kw)+...
                Lh0All(kh,kw)+LeAll(kh,kw);
            if(LtotAll(kh,kw)<minLwh)
                %%% FOR THE CUMULATIVE PROBABILITY FUNCTIONS
                ElocsCumuStore=ElocsCumu;pEcumuStore=pEcumu;EvalsCumu2Store=EvalsCumu2;
                WlocsCumuStore=WlocsCumu;WdistCumuStore=WdistCumu;WvalsCumu2store=WvalsCumu2;
                HlocsCumuStore=HlocsCumu;HdistCumuStore=HdistCumu;HvalsCumu2store=HvalsCumu2;
                %%% FOR THE PROBABILITY FUNCTIONS
                ElocsStore=Elocs;pEStore=EdistProb;EvalsStore=EvalsProb;
                WlocsStore=Wlocs;pWStore=WdistProb;WvalsStore=WvalsProb;
                HlocsStore=Hlocs;pHStore=HdistProb;HvalsStore=HvalsProb;
                %%% FOR THE DESCRIPTION LENGTH PLOTS FOR DISTRIBUTIONS AND
                %%% VALUES
                descLengths{1,1}=WdistProb4;descLengths{2,1}=WvalsProb4;
                descLengths{1,2}=HdistProb4;descLengths{2,2}=HvalsProb4;
                descLengths{1,3}=EdistProb4;descLengths{2,3}=EvalsProb4;
                %%% FOR TOTAL DESCRIPTION LENGTH WITH CHANGES IN r
                Lvals(i,2)=EvalsProb4(end);Lvals(i,3)=WvalsProb4(end);
                Lvals(i,5)=HvalsProb4(end);Lvals(i,4)=Lw0;
                Lvals(i,6)=Lh0;Lvals(i,1)=sum(Lvals(i,2:6));
                %%% FOR THE DESCRIPTION LENGTHS
                Ldist(i,2)=EdistProb4(end);Ldist(i,3)=WdistProb4(end);
                Ldist(i,5)=HdistProb4(end);Ldist(i,4)=Lw0;
                Ldist(i,6)=Lh0;Ldist(i,1)=sum(Ldist(i,2:6));
                % Reset minLwh
                minLwh=Lvals(i,1);
            end
            %%%%%%%%%%%%%%%%%%%%%%% END OF THE E VALUES CODE %%%%%%%%%%%%%%%%%%%%%%%%%%
        end %%% END OF THE kh LOOP
    end %%% END OF THE kw LOOP
    if(options.threshold)
        thresholdPlot(LtotAll,zeroBounds,plotNum)
        plotNum=plotNum+1;
    end
    if(options.eGaussianCumulative)
        eGaussianCumulativePlot(ElocsCumuStore,pEcumuStore,EvalsCumu2Store,plotNum)
        plotNum=plotNum+1;
    end
    if(options.wGammaCumulative)
        wGammaCumulativePlot(WlocsCumuStore,WdistCumuStore,WvalsCumu2store,plotNum)
        plotNum=plotNum+1;
    end
    if(options.hGammaCumulative)
        hGammaCumulativePlot(HlocsCumuStore,HdistCumuStore,HvalsCumu2store,plotNum)
        plotNum=plotNum+1;
    end
    if(options.eGaussian)
        eGaussianPlot(ElocsStore,pEStore,EvalsStore,plotNum)
        plotNum=plotNum+1;
    end
    if(options.wGamma)
        wGammaPlot(WlocsStore,pWStore,WvalsStore,plotNum)
        plotNum=plotNum+1;
    end
    if(options.hGamma)
        hGammaPlot(HlocsStore,pHStore,HvalsStore,plotNum)
        plotNum=plotNum+1;
    end
    if(options.descriptionLengths1)
        descriptionLengths1Plot(descLengths,plotNum)
        plotNum=plotNum+1;
    end
    
    %     [row,col]=find(LtotAll==min(min(LtotAll)));
    %     Ltot(i,1)=LtotAll(row(1),col(1));
    %     Le(i,1)=LeAll(row(1),col(1));
    %     Lw(i,1)=LwAll(row(1),col(1));Lh(i,1)=LhAll(row(1),col(1));
    %     Lw0(i,1)=Lw0All(row(1),col(1));Lh0(i,1)=Lh0All(row(1),col(1));
    %     display(int2str(i))
end
if(options.lValues)
    lValuesPlot(rValues,Ldist,plotNum)
    plotNum=plotNum+1;
end
if(options.lValuesNorm)
    lValuesNormPlot(rValues,Ldist,m,plotNum)
    plotNum=plotNum+1;
end

if(options.descriptionLengths2)
    descriptionLengths2Plot(rValues,Ldist,Lvals,plotNum)
end

end

function thresholdPlot(Lwh,zeroBounds,plotNum)
x=zeroBounds;y=zeroBounds;
figure(plotNum)
set(gca,'FontSize',16)
surf(x,y,Lwh)
% legend('L_{W0}+L_{W+}','L_{W0}','L_{W+}')
% legend boxoff
xlabel('W threshold value')
ylabel('H threshold value')
zlabel('Description length')
title('Zero threshold plot')
end

function eGaussianCumulativePlot(ElocsStore,pEstore,EvalsCumu2Store,plotNum)
x=ElocsStore(2:end);
yDist=pEstore; yVals=EvalsCumu2Store;
figure(plotNum)
plot(x,yVals,'kx')
hold on
plot(x,yDist,'k--')
legend('Cumulative error values','Gaussian distribution values','location','southeast')
legend boxoff
xlabel('Cumulative value')
ylabel('E terms')
title('Cumulative Probability for the Errors')
set(gca,'FontSize',16)
end

function wGammaCumulativePlot(WlocsCumuStore,WdistCumuStore,WvalsCumu2store,plotNum)
x=WlocsCumuStore(2:end);
yDist=WdistCumuStore; yVals=WvalsCumu2store;
figure(plotNum)
plot(x,yVals,'kx')
hold on
plot(x,yDist,'k--')
legend('Cumulative error values','Gamma distribution values','location','southeast')
legend boxoff
xlabel('Cumulative value')
ylabel('W terms')
title('Cumulative Probability for W')
set(gca,'FontSize',16)
end

function hGammaCumulativePlot(HlocsCumuStore,HdistCumuStore,HvalsCumu2store,plotNum)
x=HlocsCumuStore(2:end);
yDist=HdistCumuStore; yVals=HvalsCumu2store;
figure(plotNum)
plot(x,yVals,'kx')
hold on
plot(x,yDist,'k--')
legend('Cumulative error values','Gamma distribution values','location','southeast')
legend boxoff
xlabel('Cumulative value')
ylabel('H terms')
title('Cumulative Probability for H')
set(gca,'FontSize',16)
end

function eGaussianPlot(ElocsStore,pEStore,EvalsStore,plotNum)
x=ElocsStore(2:end);
yDist=pEStore; yVals=EvalsStore;
figure(plotNum)
plot(x,yVals,'kx')
hold on
plot(x,yDist,'k--')
legend('Error values','Gaussian distribution values','location','northeast')
legend boxoff
xlabel('Probability')
ylabel('E terms')
title('Probability for the Errors')
set(gca,'FontSize',16)
end

function wGammaPlot(WlocsStore,pWStore,WvalsStore,plotNum)
x=WlocsStore(2:end);
yDist=pWStore; yVals=WvalsStore;
figure(plotNum)
plot(x,yVals,'kx')
hold on
plot(x,yDist,'k--')
legend('Error values','Gamma distribution values','location','northeast')
legend boxoff
xlabel('Probability')
ylabel('W terms')
title('Probability for the W terms')
set(gca,'FontSize',16)
end

function hGammaPlot(HlocsStore,pHStore,HvalsStore,plotNum)
x=HlocsStore(2:end);
yDist=pHStore; yVals=HvalsStore;
figure(plotNum)
plot(x,yVals,'kx')
hold on
plot(x,yDist,'k--')
legend('Error values','Gamma distribution values','location','northeast')
legend boxoff
xlabel('Probability')
ylabel('H terms')
title('Probability for the H terms')
set(gca,'FontSize',16)
end

function lValuesPlot(rValues,Lvals,plotNum)
x=rValues;
y1=Lvals(:,1);y2=Lvals(:,2);y3=Lvals(:,3);y4=Lvals(:,4);y5=Lvals(:,5);
y6=Lvals(:,6);
figure(plotNum)
plot(x,y1,'k-')
hold on
plot(x,y2,'k--')
plot(x,y3,'k:')
plot(x,y4,'k-.')
plot(x,y5,'r:')
plot(x,y6,'b-.')
legend('L_{tot}','L_{E}','L_{W}','L_{W0}','L_{H}','L_{H0}','location','northeast')
legend boxoff
xlabel('r')
ylabel('Description length')
title('Description lengths')
set(gca,'FontSize',16)
end

function lValuesNormPlot(rValues,Lvals,m,plotNum)
x=rValues/m;
y1=Lvals(:,1)/max(Lvals(:,1));y2=Lvals(:,2)/max(Lvals(:,2));y3=Lvals(:,3)/max(Lvals(:,3));
y4=Lvals(:,4)/max(Lvals(:,4));y5=Lvals(:,5)/max(Lvals(:,5));y6=Lvals(:,6)/max(Lvals(:,6));
figure(plotNum)
plot(x,y1,'k-')
hold on
plot(x,y2,'k--')
plot(x,y3,'k:')
plot(x,y4,'k-.')
plot(x,y5,'r:')
plot(x,y6,'b-.')
legend('L_{tot}','L_{E}','L_{W}','L_{H}','location','northeast')
legend boxoff
xlabel('r/m')
ylabel('Normalised Description length')
title('Normalised Description lengths')
set(gca,'FontSize',16)
end

function descriptionLengths1Plot(descLengths,plotNum)
tempProbW4=descLengths{1,1};tempProbW44=descLengths{2,1};
tempProbH4=descLengths{1,2};tempProbH44=descLengths{2,2};
tempProb33=descLengths{1,3};tempProb44=descLengths{2,3};
x=linspace(0,1);x1=linspace(0,1,length(tempProb33));
y1=interp1(x1,tempProb33,x);
x2=linspace(0,1,length(tempProbH4));
y2=interp1(x2,tempProbH4,x);
x3=linspace(0,1,length(tempProbW4));
y3=interp1(x3,tempProbW4,x);
y4=y1+y2+y3;
x11=linspace(0,1,length(tempProb44));
y11=interp1(x11,tempProb44,x);
x22=linspace(0,1,length(tempProbH44));
y22=interp1(x22,tempProbH44,x);
x33=linspace(0,1,length(tempProbW44));
y33=interp1(x33,tempProbW44,x);
y44=y11+y22+y33;
figure(plotNum)
p1=plot(x,y4,'k-');
hold on
plot(x,y1,'b-')
plot(x,y3,'g-')
plot(x,y2,'r-')
p2=plot(x,y44,'k--');
plot(x,y11,'b--')
plot(x,y33,'g--')
plot(x,y22,'r--')
legend([p1 p2],'Distributions','Actual numbers','location','northwest')
legend boxoff
xlabel('Fraction of bins')
ylabel('Description Length')
title('Description lengths for approximations and values')
set(gca,'FontSize',16)
end

function descriptionLengths2Plot(rValues,Lall,Lvals,plotNum)
Ltot1=Lall(:,1);Le1=Lall(:,2);Lw1=Lall(:,3);Lw01=Lall(:,4);Lh1=Lall(:,5);Lh01=Lall(:,6);
Ltot2=Lvals(:,1);Le2=Lvals(:,2);Lw2=Lvals(:,3);Lh2=Lvals(:,5);
Lh1=Lh1+Lh01;Lw1=Lw1+Lw01;
Lw2=Lw2+Lw01;Lh3=Lh2+Lh01;
figure(plotNum)
p1=plot(rValues,Ltot1,'k-');
hold on
plot(rValues,Le1,'b-');
plot(rValues,Lw1,'g-');
plot(rValues,Lh1,'r-');
p2=plot(rValues,Ltot2,'k--');
plot(rValues,Le2,'b--');
plot(rValues,Lw2,'g--');
plot(rValues,Lh3,'r--');
legend([p1 p2],'Distributions','Actual numbers','location','northwest')
legend boxoff
xlabel('r')
ylabel('Description Length')
title('Description lengths for approximations and values for different r')
set(gca,'FontSize',16)
end




