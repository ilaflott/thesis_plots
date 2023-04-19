//// MISC RELEVANT STRINGS
//ETABIN TAGS, RELEVANT TO SYSTEMATICS STUFF
const std::string ETABIN_TAG_array[]={
  "00eta05",
  "05eta10",
  "10eta15",
  "15eta20"
};
const int ETABIN_TAG_Nstrs=sizeof(ETABIN_TAG_array)/sizeof(std::string);




//// PY8 UNFOLDING RESULTS
//DATA, NOTE: JEC/KITER SYSTEMATICS IN THESE FILES
const std::string PY8_UNFDIR_DATA   ="PY8unf/data/";
//const std::string PY8_UNF_DATA_file_base="ak4PFJets_wjtID_anabins_Bayes_PY8_FullRECO_11.26.19_SMPbins_withJECsysv2_SMPbins_ptLo56_withLumiCorr_semifinalv2_LOMC_";
const std::string PY8_UNF_DATA_file_base="ak4PFJets_wjtID_anabins_Bayes_PY8_FullRECO_11.26.19_SMPbins_withJECsysv2_SMPbins_ptLo56_withLumiCorr_semifinalv2_LowHLT40Thresh_LOMC_";
const std::string PY8_UNFDIR_DATA_file_array[]={
  PY8_UNF_DATA_file_base+"00eta05",
  PY8_UNF_DATA_file_base+"05eta10",
  PY8_UNF_DATA_file_base+"10eta15",
  PY8_UNF_DATA_file_base+"15eta20"
};
const std::string PY8_UNFDIR_DATA_shortdir_array[]={
  "YBIN0",
  "YBIN1",
  "YBIN2",
  "YBIN3"
};
const int PY8_UNFDIR_DATA_Nfiles=sizeof( PY8_UNFDIR_DATA_file_array )/sizeof(std::string);
//CLOSURE
const std::string PY8_UNFDIR_CLOSURE="PY8unf/closure/";
const std::string PY8_UNF_CLOSURE_file_base="ak4PFJets_wjtID_anabins_Bayes_Closure_PY8_FullRECO_06.25.19_SMPbins_LOMC_";
const std::string PY8_UNFDIR_CLOSURE_file_array[]={
  PY8_UNF_CLOSURE_file_base+"00eta05",
  PY8_UNF_CLOSURE_file_base+"05eta10",
  PY8_UNF_CLOSURE_file_base+"10eta15",
  PY8_UNF_CLOSURE_file_base+"15eta20"
};
const std::string PY8_UNFDIR_CLOSURE_shortdir_array[]={
  "YBIN0",
  "YBIN1",
  "YBIN2",
  "YBIN3"
};
const int PY8_UNFDIR_CLOSURE_Nfiles=sizeof( PY8_UNFDIR_CLOSURE_file_array )/sizeof(std::string);





//// SMEARED NLO UNFOLDING RESULTS
//DATA, NOTE: FILES w/ JEC/KITER SYSTEMATICS ARE USED HERE
const std::string NLO_UNFDIR_DATA   ="NLOunf/data/";
//const std::string NLO_UNF_DATA_file_base="ak4PFJets_wjtID_anabins_11.26.19_Bayes_NNPDF_NLO_sigmu_noJERscales_withLumiCorr_SMPbins_withJECsysv2_ptLo56_semifinalv2_JECsysv2_NLOMC_wNP_";
const std::string NLO_UNF_DATA_file_base="ak4PFJets_wjtID_anabins_11.26.19_Bayes_NNPDF_NLO_sigmu_noJERscales_withLumiCorr_SMPbins_withJECsysv2_ptLo56_semifinalv2_LowHLT40Thresh_JECsysv2_NLOMC_wNP_";
const std::string NLO_UNFDIR_DATA_file_array[]={ 
  NLO_UNF_DATA_file_base+"00eta05",
  NLO_UNF_DATA_file_base+"05eta10",
  NLO_UNF_DATA_file_base+"10eta15",
  NLO_UNF_DATA_file_base+"15eta20"
};
const std::string NLO_UNFDIR_DATA_shortdir_array[]={ 
  "YBIN0",
  "YBIN1",
  "YBIN2",
  "YBIN3"
};
const int   NLO_UNFDIR_DATA_Nfiles=sizeof(   NLO_UNFDIR_DATA_file_array )/sizeof(std::string);
//CLOSURE
const std::string NLO_UNFDIR_CLOSURE="NLOunf/closure/";
const std::string NLO_UNF_CLOSURE_file_base="ak4PFJets_wjtID_anabins_06.25.19_Bayes_Closure_NNPDF_NLO_sigmu_noJERscales_semifinal_NLOMC_wNP_";
const std::string NLO_UNFDIR_CLOSURE_file_array[]={ 
  NLO_UNF_CLOSURE_file_base+"00eta05",
  NLO_UNF_CLOSURE_file_base+"05eta10",
  NLO_UNF_CLOSURE_file_base+"10eta15",
  NLO_UNF_CLOSURE_file_base+"15eta20"
};
const std::string NLO_UNFDIR_CLOSURE_shortdir_array[]={ 
  "YBIN0",
  "YBIN1",
  "YBIN2",
  "YBIN3"
};
const int   NLO_UNFDIR_CLOSURE_Nfiles=sizeof(   NLO_UNFDIR_CLOSURE_file_array )/sizeof(std::string);

//DATA SYST, NOTE: FILES HAVE ETA STRING AT END REMOVED SO ETABIN MUST BE SPECIFIED IN SCRIPTS
//const std::string NLO_UNF_DATA_SYST_file_base="ak4PFJets_wjtID_anabins_11.26.19_Bayes_NNPDF_NLO_sigmu_noJERscales_withLumiCorr_SMPbins_withJECsysv2_ptLo56_semifinalv2_";
const std::string NLO_UNF_DATA_SYST_file_base="ak4PFJets_wjtID_anabins_11.26.19_Bayes_NNPDF_NLO_sigmu_noJERscales_withLumiCorr_SMPbins_withJECsysv2_ptLo56_semifinalv2_LowHLT40Thresh_";
const std::string NLO_UNFDIR_DATA_SYST_file_array[]={ 
  NLO_UNF_DATA_SYST_file_base+"JECsys_NLOMC_wNP_",
  NLO_UNF_DATA_SYST_file_base+"JERsys_NLOMC_wNP_",
  NLO_UNF_DATA_SYST_file_base+"PDFsysupdown_NLOMC_wNP_",
  NLO_UNF_DATA_SYST_file_base+"PDF12_NLOMC_wNP_",
  NLO_UNF_DATA_SYST_file_base+"NPsysupdown_NLOMC_wNP_",
  NLO_UNF_DATA_SYST_file_base+"NPsys12_NLOMC_wNP_"
};
const int   NLO_UNFDIR_DATA_SYST_Nfiles=(sizeof(   NLO_UNFDIR_DATA_SYST_file_array )/sizeof(std::string))*(ETABIN_TAG_Nstrs);//MULTIPLY BY # OF ETABINS


//OLD, THIS SHOULD BE REDONE
const std::string NLO_UNFDIR_DATA_JOHNS_BINS   ="NLOunf/MY_RESULTS_JOHNS_BINS/";
const std::string NLO_UNF_DATA_JOHNS_BINS_file_base="ak4PFJets_wjtID_anabins_08.06.19_Bayes_NNPDF_NLO_sigmu_noJERscales_withLumiCorr_JOHNbins_semifinal_JECsys_NLOMC_wNP_";
const std::string NLO_UNFDIR_DATA_JOHNS_BINS_file_array[]={ 
  NLO_UNF_DATA_JOHNS_BINS_file_base+"00eta05",
  NLO_UNF_DATA_JOHNS_BINS_file_base+"05eta10",
  NLO_UNF_DATA_JOHNS_BINS_file_base+"10eta15",
  NLO_UNF_DATA_JOHNS_BINS_file_base+"15eta20"
};
const int   NLO_UNFDIR_DATA_JOHNS_BINS_Nfiles=sizeof(   NLO_UNFDIR_DATA_JOHNS_BINS_file_array )/sizeof(std::string);

//COMPARISON W JOHNS RESULTS
//const std::string NLO_UNFDIR_JOHNS_DATA="smearNLO_unfolding/comparison/JOHNS_DATA_OLDJEC/";//tail discrepency
//const std::string MEAS_JOHNS_DATA_SPECTRA="mergeHistos.root";
const std::string NLO_UNF_JOHNS_DATA_SPECTRA="allUnfoldedPlots.root";
const std::string NLO_UNF_JOHNS_DATA_RATIOS="allUnfoldedRatios.root";
//const std::string NLO_UNFDIR_JOHNS_DATA="smearNLO_unfolding/comparison/JOHNS_DATA_1/";//spike in me/john ratio at low pt
//const std::string MEAS_JOHNS_DATA_SPECTRA="mergeHistos_newJEC.root";
const std::string NLO_UNFDIR_JOHNS_DATA="NLOunf/comparison/JOHNS_DATA_2/";
const std::string MEAS_JOHNS_DATA_SPECTRA="makeSpectra.root";//latest, lets see

//// MC JERS
const std::string JERDIR_MC  ="JER/";
const std::string JERDIR_MC_file_array[]={ 
  "ak4PF_PY8JER_00eta05_06.25.19_sigmu_semifinal",
  "ak4PF_PY8JER_05eta10_06.25.19_sigmu_semifinal",
  "ak4PF_PY8JER_10eta15_06.25.19_sigmu_semifinal",
  "ak4PF_PY8JER_15eta20_06.25.19_sigmu_semifinal"
};
const std::string JERDIR_MC_shortdir_array[]={ 
  "YBIN0",
  "YBIN1",
  "YBIN2",
  "YBIN3"
};
const int   JERDIR_MC_Nfiles=sizeof(   JERDIR_MC_file_array )/sizeof(std::string);

/// /MC VTX WGTS
const std::string VTXWDIR_MC  ="PY8vzwgts/";
const std::string VTXWDIR_MC_file_array[]={ 
  "Py8_vzEvtWeights_baryshift_neg"
};
const std::string VTXWDIR_MC_shortdir_array[]={ 
  ""
};
const int  VTXWDIR_MC_Nfiles=sizeof(  VTXWDIR_MC_file_array )/sizeof(std::string);



//// JET PLOTS QA
const std::string JETQA  ="jetPlots/";
const std::string JETQA_file_array[]={ 
  //"ak4PF_HPtJetTrig_semiOffPy8_00eta20_11.26.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_jetPlots"
  "ak4PF_HPtJetTrig_semiOffPy8_00eta20_11.26.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh_jetPlots"
};
const std::string JETQA_shortdir_array[]={ 
  //"ak4PF_HPtJetTrig_semiOffPy8_00eta20_11.26.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_jetPlots"
  ""
};
const int       JETQA_Nfiles=sizeof(       JETQA_file_array )/sizeof(std::string);

//// JET ID QA
const std::string JETIDQA  ="JetIDHIN/";
const std::string JETIDQA_file_array[]={ 
  //"ak4PF_HPtJetTrig_00eta20_11.26.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_jetIDPlots",
  "ak4PF_HPtJetTrig_00eta20_11.26.19_SMPbins_wJECsysv2_withLumiCorr_semifinalv2_LowHLT40Thresh_jetIDPlots",
  "ak4PF_Py8_CUETP8M1_00eta20_06.25.19_SMPbins_semifinal_jetIDPlots"
};
const std::string JETIDQA_shortdir_array[]={ 
  "DATA",
  "PY8"
};
const int     JETIDQA_Nfiles=sizeof(     JETIDQA_file_array )/sizeof(std::string);

//// JET TRIG QA
const std::string JETTRIGQA  ="jetTrigEff/";
const std::string JETTRIGQA_file_array[]={ 
  //"ak4PF_HPtJetTrig_00eta51_11.26.19_HLT60Eff_useHLT40Ref_semifinalv2_jetTrigEff",
  "ak4PF_HPtJetTrig_00eta51_11.26.19_HLT60Eff_useHLT40Ref_semifinalv2_LowHLT40Thresh_jetTrigEff",
  //"ak4PF_HPtJetTrig_00eta51_11.26.19_HLT80Eff_useHLT60Ref_semifinalv2_jetTrigEff"
  "ak4PF_HPtJetTrig_00eta51_11.26.19_HLT80Eff_useHLT60Ref_semifinalv2_LowHLT40Thresh_jetTrigEff"
};
const std::string JETTRIGQA_shortdir_array[]={ 
  "HLT60",
  "HLT80"
};
const int   JETTRIGQA_Nfiles=sizeof(   JETTRIGQA_file_array )/sizeof(std::string);
//const int Nfiles=sizeof(file_array)/sizeof(std::string);

//// NON PERTURBATIVE CORRECTIONS
const std::string NPC_DIR="NLOSpectra_NPCs/NPCorr5TeV/";
const std::string NPC_FILE="NLOpNP_InclusiveJets5TeV";






