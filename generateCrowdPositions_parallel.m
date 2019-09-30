function [modelIrrev, matchRev, rev2irrev, irrev2rev, solution, EnzConSamples] = generateCrowdPositions_parallel(model,C,biomass,num)
% INPUT
%  model             COBRA model structure with 3 additional vectors of same size as 'rxns':
%                    kcat_f, kcat_b, molwt (if any of the value unknown,
%                    provide '0')
%                    kcat units should be '1/s' and mol wt in 'Dalton'
%  frac              fraction of enzymatic mass in overall dry cell weight   
%  biomass           name of biomass reaction (to be excluded from enzyme
%                    capacity flux constraint)
% 
% OUTPUT
%  modelIrrev    Model in irreversible format with enzyme constraint added
%                as pseudo reaction
%  solution
%    f         Objective value
%    x         Primal
%    y         Dual
%    w         Reduced costs
%    s         Slacks
%    stat      Solver status in standardized form
%               1   Optimal solution
%               2   Unbounded solution
%               0   Infeasible
%              -1  No solution reported (timelimit, numerical problem etc)
% Lokanand Koduru            10/03/18
% Meiyappan Lakshmanan       10/04/18 Generalized the code

%% convert to irreversible format
[modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(model);
% sol = optimizeCbModel(modelIrrev);
% modelIrrev = changeRxnBounds(modelIrrev,'biomass_reaction',0.1*sol.f,'b');
modelIrrev.description='model_with_EnzCon';
revFlag='false';
%% identify reactions with kcat values
kcat_rxns = model.rxns((model.kcat)~=0);
%% filter out the mol wt of reactions with kcat values
mw = model.mw(find(ismember(model.rxns,kcat_rxns)));
%% identify reactions in Irreversible model with kcat values
kcat_rxns_irrev = modelIrrev.rxns(ismember(irrev2rev,find(ismember(model.rxns,kcat_rxns))));
%% assign kcat values to reactions in Irreversible model with kcat values
kcat_irrev = zeros(length(kcat_rxns_irrev),1);
for i=1:1:length(kcat_rxns_irrev)
    if model.kcat(irrev2rev(find(ismember(modelIrrev.rxns,kcat_rxns_irrev(i))))) ~= 0
        kcat_irrev(i) = model.kcat(irrev2rev(find(ismember(modelIrrev.rxns,kcat_rxns_irrev(i)))));
    end
end
kcat_mw = zeros(length(kcat_rxns_irrev),1);
%% calculate kcat/mw values
for i=1:1:length(kcat_rxns_irrev)
    kcat_mw(i) = (((mw(find(ismember(kcat_rxns,model.rxns(irrev2rev(find(ismember(modelIrrev.rxns,kcat_rxns_irrev(i))))))))/1000)*0.73)/(kcat_irrev(i)*3600))*C;
end

kcat_mw = kcat_mw(~isnan(kcat_mw));
kcat_rxns_irrev = kcat_rxns_irrev(~isnan(kcat_mw));
%% add Dummy reaction to model repreenting the enzyme constraint
lowerBound=0;
upperBound=1;
modelIrrev=addReaction(modelIrrev,'EnzCon_C_Rxn',{'EnzCon_Met[c]'},-1,revFlag,lowerBound,upperBound,0,'','','','');
EnzCon_MetInd=find(ismember(modelIrrev.mets,'EnzCon_Met[c]'));
Rxn_WithEnzConInd=find(ismember(modelIrrev.rxns,kcat_rxns_irrev));
selExc=findExcRxns(modelIrrev);
ExRxnInd=find(selExc);
BiomassRxnInd=find(ismember(modelIrrev.rxns,biomass));
Solutions_EnzConSamples=zeros(length(modelIrrev.rxns),num);
% Solutions_EnzConSamples=zeros(length(modelIrrev.rxns)+1,num);
EnzConSamples=zeros(length(modelIrrev.rxns),num);
RandomEnzConCoeff=zeros(length(modelIrrev.rxns),num);
for i=1:1:num
    RandomEnzConCoeff(:,i)=randsample(kcat_mw,length(modelIrrev.rxns),true);
end

parfor j=1:1:num
    modelIrrev1 = modelIrrev;
    kcat_mw1 = kcat_mw;
    modelIrrev1.S(EnzCon_MetInd,:)=transpose(RandomEnzConCoeff(:,j));
    %% Removing EnzCon coeff for biomass equation
    modelIrrev1.S(EnzCon_MetInd,BiomassRxnInd)=0;

    %% Replacing the randomly assigned EnzCon coeff with the original coefficients for those reactions whose EnzCon coeffs are available (calculated list) 
    for k=1:length(Rxn_WithEnzConInd)
        if(isempty(modelIrrev1.S(EnzCon_MetInd,Rxn_WithEnzConInd(k))))~=0
         modelIrrev1.S(EnzCon_MetInd,Rxn_WithEnzConInd(k))=kcat_mw1(find(ismember(kcat_rxns_irrev,modelIrrev1.rxns(Rxn_WithEnzConInd(k)))));
        end
    end
    %% Removing EnzCon coeff for exchange reactions
    modelIrrev1.S(EnzCon_MetInd,ExRxnInd)=0;
    modelIrrev1.S(EnzCon_MetInd,find(ismember(modelIrrev1.rxns,'EnzCon_C_Rxn')))=-1;
    changeCobraSolver('gurobi7');
    solution=optimizeCbModel(modelIrrev1,'max');
    if(isempty(solution.x))==0
        Solutions_EnzConSamples(:,j)=solution.x;
        EnzConSamples(:,j) = modelIrrev1.S(EnzCon_MetInd,:)';
    end
    
%     if(isempty(solution.x))==0
%         modelIrrev1.lb(modelIrrev1.c==1) = solution.f;
%         modelIrrev1.S(end+1,1:end-1) = ones(size(modelIrrev1.S(1,1:end-1)));
%         modelIrrev1.b(end+1) = 0;
%         modelIrrev1.mets{end+1} = 'fluxMeasure';
%         modelIrrev1 = addReaction(modelIrrev1,'netFlux',{'fluxMeasure'},[-1],false,0,1000,0,'','');
%         modelIrrev1.c = zeros(length(modelIrrev1.rxns),1);
%         modelIrrev1 = changeObjective(modelIrrev1, 'netFlux');
%         MinimizedFlux = optimizeCbModel(modelIrrev1,'min');
%         if(isempty(MinimizedFlux.x))==0
%             Solutions_EnzConSamples(:,j)=MinimizedFlux.x;
%         end
%     end
end
solution = Solutions_EnzConSamples;
end