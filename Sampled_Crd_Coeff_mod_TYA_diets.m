function [targetProds, Solutions_GR, Solutions_ProdRate] = Sampled_Crd_Coeff_mod_TYA_diets(model,crd_val,ex_values, ex_fluxes, diet_name)
temp_model1 = addReaction(model,'EX_TempProd',{'A'},[-1],false);
% [selExc,~] = findExcRxns(temp_model);
% temp_model.rev(selExc) = 1;
%replace metabolite names
new_mets = {'LTA[c]','LTA[c]','LTA[c]','LTA[c]','LTA[c]','LTA[c]','peptido[c]','peptido[c]','peptido[c]','peptido[c]','peptido[c]','CRB[c]','CRB[c]','CRB[c]','CRB[c]','CRB[c]','CRB[c]'};
original_mets = {'LTA_LCA[c]','LTA_LPL[c]','LTA_LME[c]','LTA_LSA[c]','LTA_LFE[c]','LTAAlaGal_LLA[c]','peptido_LCA[c]','peptido_LPL[c]','peptido_LME[c]','peptido_LSA[c]','peptido_LFE[c]','CRB_LCA[c]','CPS_LPL2[c]','CRB_LME[c]','CRB_LSA[c]','CRB_LFE[c]','CPS_LLA[c]'};
for i=1:length(new_mets)
    met = original_mets(1,i);
    A = find(ismember(temp_model1.mets,met));
    if ~isempty(A)
        temp_model1.mets(A,1)= new_mets(1,i);
    end
end
for i=1:1:length(diet_name)
    temp_model = temp_model1;
    temp_model = constraindiets(temp_model,ex_fluxes,ex_values(:,i));
    temp_model = changeRxnBounds(temp_model,'EX_o2(e)',0,'l');
    temp_model = changeRxnBounds(temp_model,'ATPH',0,'l');
    [modelIrrev, matchRev, rev2irrev, irrev2rev, solution, EnzConSamples] = generateCrowdPositions_parallel(temp_model,1/crd_val,'biomass',5000);
    NonZero_BiomassCrdInds=find((solution(find(ismember(irrev2rev,find(ismember(temp_model.rxns,'biomass')))),:)));
    L_LactateExchInd=find(ismember(irrev2rev,find(ismember(temp_model.rxns,'EX_lac-L(e)'))));
    D_LactateExchInd=find(ismember(irrev2rev,find(ismember(temp_model.rxns,'EX_lac-D(e)'))));
    NonZero_Total_LactateCrdInds=unique([find((solution(L_LactateExchInd,:)));find((solution(D_LactateExchInd,:)))])';
    NonZero_Biomass_Lactate_CrdInds = intersect(NonZero_BiomassCrdInds,NonZero_Total_LactateCrdInds);
    [Nbins,Edges]=histcounts(log(solution(1,NonZero_Biomass_Lactate_CrdInds)),80);
    minGR=Edges((find(Nbins==max(Nbins))-2));
    maxGR=Edges((find(Nbins==max(Nbins))));
    MinMaxGRIndsOf_NonZero_Biomass_Lactate_CrdInds=[];
    for h=1:length(NonZero_Biomass_Lactate_CrdInds)
        if log(solution(1,NonZero_Biomass_Lactate_CrdInds(h)))>=minGR
            if log(solution(1,NonZero_Biomass_Lactate_CrdInds(h)))<=maxGR
                MinMaxGRIndsOf_NonZero_Biomass_Lactate_CrdInds=[MinMaxGRIndsOf_NonZero_Biomass_Lactate_CrdInds,NonZero_Biomass_Lactate_CrdInds(h)];
            end
        end
    end
    Best_Crd_positions = MinMaxGRIndsOf_NonZero_Biomass_Lactate_CrdInds;
    [targetProds(i).targetProds, Solutions_GR(i).Solutions_GR, Solutions_ProdRate(i).Solutions_ProdRate] = postbiotics_eval(temp_model, EnzConSamples, Best_Crd_positions);
end
