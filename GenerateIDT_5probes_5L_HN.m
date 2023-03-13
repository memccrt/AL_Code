%IDT PROBES
%5L HN
load DataForSims_IDT.mat
[OverallMCC_5L_HN,ProbesRemoved_5L_HN,Opt_Thresh_probes_5L_HN,Good_probes_final_5L_HN,R_Cond_final_5L_HN,ProbesRemoved_names_5L_HN] = RemoveTwoBarcodes_onebyone(probes,R_5L,high_noise,thresh,names)
save('IDT_5Probes_5L_HN_6-15-22(3)')