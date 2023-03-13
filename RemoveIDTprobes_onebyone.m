function [OverallMCC,ProbesRemoved,Opt_Thresh_probes,Good_probes_final,R_Cond_final,ProbesRemoved_names] = RemoveIDTprobes_onebyone(Good_probes,R_Cond_final,noise,thresh,names)
%Inputs:    Good_probes- List of good probes from RemoveProbes_onebyone
%           R_Cond_final- Condensed reference matrix for Good_probes
%           noise- noise level 
%           thresh- vector of thresholds for binary classification
%           names- list of names for fluorescent proteins 
%Outputs:   OverallMCC- Vector containing the overall MCC value from each loop
%           ProbesRemoved- List of probes removed in each loop
%           Opt_Thresh_probes- Vector containing the optimum threshold for binary classification for each probe
%           Good_probes_final- list of final good probes (Overall MCC = 1)
%           R_Cond_Barcodes_final- Condensed reference matrix for Good_probes_final
%           ProbesRemoved_names- list of probe names for Good_probes_final

ProbesRemoved = zeros(size(Good_probes));
OverallMCC = zeros(length(Good_probes),1);
    for i=1:length(Good_probes)
        [fhatmat,fmat]=f_and_fhat_5probes_hybrid(Good_probes,R_Cond_final,noise,3);
        fhatbin = BinClass(fhatmat,thresh);
        conf_mats = calc_confusion(fhatbin,fmat);
        MCC_Mat = GetMCC(conf_mats);
        [Opt_Thresh_probes,MCCMax] = Get_OptThresh(MCC_Mat,thresh);
        [MCC_Overall]=CalcOverallClassificationMetrics(fmat,fhatmat,Opt_Thresh_probes);
        [~,MinInd] = min(MCCMax);
        ProbesRemoved(i,:) = Good_probes(MinInd,:);
        Good_probes(MinInd,:)=[];
        OverallMCC(i)=(MCC_Overall);
        R_Cond_final(:,MinInd)=[];
        if MCC_Overall == 1
            break
        end 
    end
    ProbesRemoved_names=GetIndividualProbeNames(names,ProbesRemoved);
    [fhatmat,fmat]=f_and_fhat_twobarcodes_hybrid(Good_probes,R_Cond_final,noise,3);
    fhatbin = BinClass(fhatmat,thresh);
    conf_mats = calc_confusion(fhatbin,fmat);
    MCC_Mat = GetMCC(conf_mats);
    [Opt_Thresh_probes] = Get_OptThresh(MCC_Mat,thresh);
    R_Cond_final = R_Cond_final;
    Good_probes_final=Good_probes;