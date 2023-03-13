%IDT PROBES
%3L HN
load DataForSims_IDT.mat
[OverallMCC_3L_HN,ProbesRemoved_3L_HN,Opt_Thresh_probes_3L_HN,Good_probes_final_3L_HN,R_Cond_final_3L_HN,ProbesRemoved_names_3L_HN] = RemoveIDTprobes_onebyone(probes,R_3L,high_noise,thresh,names)
save('IDT_5Probes_3L_HN_6-15-22(3)')