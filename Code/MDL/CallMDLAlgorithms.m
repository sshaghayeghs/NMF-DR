clear
close all
%%% Inputs: original datamatrix V; precision of the data; list
%%% of rValues, cell of Ws for each r, cell of Hs for each r, options is a
%%% structure describing which plots to output
%%% For setThreshold also require the chosen thresholds tW,tH,
V=    ;precision = ; Wstore=  ; Hstore= ;rValues=  ; 
tW=  ; tH= ; %%% Only necessary for set thresholds code
options=struct('threshold',false,'wGamma',false,'eGaussian',false,...
    'hGamma',false,'wGammaCumulative',false,'hGammaCumulative',false,...
    'eGaussianCumulative',false,'lValues',false,'lValuesNorm',false,...
    'descriptionLengths1',false,'descriptionLengths2',false);
[~,LvalsSliding]=MDLslidingThresholdHistograms(V,precision,Wstore,Hstore,rValues,options);
[LdistSliding,~]=MDLslidingThresholdDistributions(V,precision,Wstore,Hstore,rValues,options);
[LdistSet,LvalsSet]=MDLsetThreshold(V,precision,Wstore,Hstore,rValues,tW,tH,options);

