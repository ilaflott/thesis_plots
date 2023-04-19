void printEvtCountTable(std::string outdir);
void printpTbins(std::string outdir);
void  makeSMPInclJetXsec_fracstaterrData_ratios(std::string outdir, TFile* fout=NULL);
void printEvtCountTable(std::string);
void printpTbins(std::string);
void  makeSMPInclJetXsec_covmatData_onePadOneEta (std::string outdir, TFile* fout=NULL);
void  makeSMPInclJetXsec_PY8vNLOunfdata_ratios (std::string outdir, TFile* fout=NULL);//SMEARED NLO v PY8 unf data, ratios
void  makeSMPInclJetXsec_PY8vNLOunfdata_modelsys_ratios (std::string outdir, TFile* fout=NULL, bool drawSysEnvelope=false, bool drawSymmEnvelope=false );//SMEARED NLO v PY8 unf data, ratios
//void  makeSMPInclJetXsec_PY8vNLOunfdata_modelsys_ratios (std::string outdir, TFile* fout=NULL);//SMEARED NLO v PY8 unf data, ratios











//--------------------------------------------------------------------------------------------------------------------------------

void printEvtCountTable(std::string outdir){
  std::cout<<"running printEvtCountTable"<<std::endl;

  bool funcDebug=true;
  if(funcDebug)std::cout<<"in printEvtCountTable"<<std::endl;
  if(funcDebug)std::cout<<"outdir is"<<outdir <<std::endl;
  
  TFile* fin =TFile::Open((EVTCOUNTS+EVTCOUNTS_file+".root").c_str(),"READ");
  TFile* fin2=TFile::Open((JETQA+JETQA_file_array[0]+".root").c_str(),"READ");
  

  TH1D* h_NEvents_Jet80_read   =(TH1D*)fin->Get("NEvents_Jet80_read"   );
  TH1D* h_NEvents_Jet80_skimCut=(TH1D*)fin->Get("NEvents_Jet80_skimCut");
  TH1D* h_NEvents_Jet80_vzCut  =(TH1D*)fin->Get("NEvents_Jet80_vzCut"  );
  TH1D* h_NEvents_Jet80_PFMETfracCut=(TH1D*)fin->Get("NEvents_Jet80_PFMETfracCut"  );
  TH1D* h_NEvents_Jet80_is80   =(TH1D*)fin->Get("NEvents_Jet80_is80"   );

  double Nevts_jet80_read   = h_NEvents_Jet80_read   ->GetBinContent(1);
  double Nevts_jet80_skimcut= h_NEvents_Jet80_skimCut->GetBinContent(1);
  double Nevts_jet80_vzcut  = h_NEvents_Jet80_vzCut  ->GetBinContent(1);
  double Nevts_jet80_PFMETfracCut  = h_NEvents_Jet80_PFMETfracCut  ->GetBinContent(1);
  double Nevts_jet80_is80   = h_NEvents_Jet80_is80   ->GetBinContent(1);

  double perc_Nevts_jet80_read   =100.*(Nevts_jet80_read/Nevts_jet80_read); 
  double perc_Nevts_jet80_skimcut=100.*(Nevts_jet80_skimcut/Nevts_jet80_read); 
  double perc_Nevts_jet80_vzcut  =100.*(Nevts_jet80_vzcut/Nevts_jet80_read); 
  double perc_Nevts_jet80_PFMETfracCut  =100.*(Nevts_jet80_PFMETfracCut/Nevts_jet80_read); 
  double perc_Nevts_jet80_is80   =100.*(Nevts_jet80_is80/Nevts_jet80_read); 
  
  std::cout<<std::endl;
  std::cout<<" ----- Jet80 PD Event Count Report ----- "<<std::endl;
  std::cout<<"# Events ............................................ in Dataset: "<<  Nevts_jet80_read   <<", Fraction of total: "<< perc_Nevts_jet80_read   <<std::endl;
  std::cout<<"# Events after HBHENoise, BeamScraping, and PV  quality  filters: "<<  Nevts_jet80_skimcut<<", Fraction of total: "<< perc_Nevts_jet80_skimcut<<std::endl;
  std::cout<<"# Events after ................................. |vz|< 24 cm cut: "<<  Nevts_jet80_vzcut  <<", Fraction of total: "<< perc_Nevts_jet80_vzcut  <<std::endl;
  std::cout<<"# Events after ............................. PFMETFrac < 0.3 cut: "<<  Nevts_jet80_PFMETfracCut  <<", Fraction of total: "<< perc_Nevts_jet80_PFMETfracCut  <<std::endl;
  std::cout<<"# Events after ............... passing trigger jet req for HLT80: "<<  Nevts_jet80_is80   <<", Fraction of total: "<< perc_Nevts_jet80_is80   <<std::endl;
  
  TH1D* h_NEvents_LowJets_read   =(TH1D*)fin->Get("NEvents_LowJets_read"   );
  TH1D* h_NEvents_LowJets_skimCut=(TH1D*)fin->Get("NEvents_LowJets_skimCut");
  TH1D* h_NEvents_LowJets_vzCut  =(TH1D*)fin->Get("NEvents_LowJets_vzCut"  );
  TH1D* h_NEvents_LowJets_PFMETfracCut=(TH1D*)fin->Get("NEvents_LowJets_PFMETfracCut"  );
  TH1D* h_NEvents_LowJets_is40   =(TH1D*)fin->Get("NEvents_LowJets_is40"   );
  TH1D* h_NEvents_LowJets_is60   =(TH1D*)fin->Get("NEvents_LowJets_is60"   );

  double Nevts_lowjets_read   = h_NEvents_LowJets_read   ->GetBinContent(1);
  double Nevts_lowjets_skimcut= h_NEvents_LowJets_skimCut->GetBinContent(1);
  double Nevts_lowjets_vzcut  = h_NEvents_LowJets_vzCut  ->GetBinContent(1);
  double Nevts_lowjets_PFMETfracCut  = h_NEvents_LowJets_PFMETfracCut  ->GetBinContent(1);
  double Nevts_lowjets_is60   = h_NEvents_LowJets_is60   ->GetBinContent(1);
  double Nevts_lowjets_is40   = h_NEvents_LowJets_is40   ->GetBinContent(1);
  double Nevts_lowjets_is60or40   = Nevts_lowjets_is60+Nevts_lowjets_is40;
 
  double perc_Nevts_lowjets_read   =100.*(Nevts_lowjets_read/Nevts_lowjets_read); 
  double perc_Nevts_lowjets_skimcut=100.*(Nevts_lowjets_skimcut/Nevts_lowjets_read); 
  double perc_Nevts_lowjets_vzcut  =100.*(Nevts_lowjets_vzcut/Nevts_lowjets_read);
  double perc_Nevts_lowjets_PFMETfracCut  =100.*(Nevts_lowjets_PFMETfracCut/Nevts_lowjets_read); 
  double perc_Nevts_lowjets_is60   =100.*(Nevts_lowjets_is60/Nevts_lowjets_read); 
  double perc_Nevts_lowjets_is40   =100.*(Nevts_lowjets_is40/Nevts_lowjets_read); 
  double perc_Nevts_lowjets_is60or40   =100.*(Nevts_lowjets_is60or40/Nevts_lowjets_read);  

  std::cout<<std::endl;
  std::cout<<" ----- LowerJets PD Event Count Report ----- "<<std::endl;
  std::cout<<"# Events ............................................ in Dataset: "<< Nevts_lowjets_read     <<", Fraction of total: "<< perc_Nevts_lowjets_read     <<std::endl;
  std::cout<<"# Events after HBHENoise, BeamScraping, and PV  quality  filters: "<< Nevts_lowjets_skimcut  <<", Fraction of total: "<< perc_Nevts_lowjets_skimcut  <<std::endl;
  std::cout<<"# Events after ..................................|vz|< 24 cm cut: "<< Nevts_lowjets_vzcut    <<", Fraction of total: "<< perc_Nevts_lowjets_vzcut    <<std::endl;
  std::cout<<"# Events after ............................. PFMETFrac < 0.3 cut: "<<  Nevts_lowjets_PFMETfracCut  <<", Fraction of total: "<< perc_Nevts_lowjets_PFMETfracCut  <<std::endl;
  std::cout<<"# Events after ...... passing trigger jet req for HLT60 or HLT40: "<< Nevts_lowjets_is60or40 <<", Fraction of total: "<< perc_Nevts_lowjets_is60or40 <<std::endl;
  std::cout<<"# Events ......................... passing trigger jet req for HLT40 only: "<<Nevts_lowjets_is40 <<", Fraction of total: "<< perc_Nevts_lowjets_is40 <<std::endl;
  std::cout<<"# Events ......................... passing trigger jet req for HLT60 only: "<<Nevts_lowjets_is60 <<", Fraction of total: "<< perc_Nevts_lowjets_is60 <<std::endl;
  
  TH1D* h_NEvents_MC_read   =(TH1D*)fin->Get("NEvents_MC_read"   );
  TH1D* h_NEvents_MC_skimCut=(TH1D*)fin->Get("NEvents_MC_skimCut");
  TH1D* h_NEvents_MC_vzCut  =(TH1D*)fin->Get("NEvents_MC_vzCut"  );
  TH1D* h_NEvents_MC_PFMETfracCut  =(TH1D*)fin2->Get("MCTH1_hpthatWeightedVz_rebin"  );
  
  double Nevts_mc_read   = h_NEvents_MC_read   ->GetBinContent(1);
  double Nevts_mc_skimcut= h_NEvents_MC_skimCut->GetBinContent(1);
  double Nevts_mc_vzcut  = h_NEvents_MC_vzCut  ->GetBinContent(1);
  double Nevts_mc_PFMETfracCut  = h_NEvents_MC_PFMETfracCut  ->GetEntries();

  double perc_Nevts_mc_read   = 100.*(Nevts_mc_read/Nevts_mc_read);
  double perc_Nevts_mc_skimcut= 100.*(Nevts_mc_skimcut/Nevts_mc_read);
  double perc_Nevts_mc_vzcut  = 100.*(Nevts_mc_vzcut/Nevts_mc_read);
  double perc_Nevts_mc_PFMETfracCut  = 100.*(Nevts_mc_PFMETfracCut/Nevts_mc_read);

  std::cout<<std::endl;
  std::cout<<" ----- PYTHIA8 PD Event Count Report ----- "<<std::endl;
  std::cout<<"# Events ............................................ in Dataset: "<< Nevts_mc_read     <<", Fraction of total: "<< perc_Nevts_mc_read     <<std::endl;
  std::cout<<"# Events after HBHENoise, BeamScraping, and PV  quality  filters: "<< Nevts_mc_skimcut  <<", Fraction of total: "<< perc_Nevts_mc_skimcut  <<std::endl;
  std::cout<<"# Events after ................................. |vz|< 24 cm cut: "<< Nevts_mc_vzcut    <<", Fraction of total: "<< perc_Nevts_mc_vzcut    <<std::endl;
  std::cout<<"# Events after ............................. PFMETFrac < 0.3 cut: "<< Nevts_mc_PFMETfracCut    <<", Fraction of total: "<< perc_Nevts_mc_PFMETfracCut    <<std::endl;


  h_NEvents_Jet80_read    ->Delete();
  h_NEvents_Jet80_skimCut ->Delete();
  h_NEvents_Jet80_vzCut	  ->Delete();
  h_NEvents_Jet80_PFMETfracCut	  ->Delete();
  h_NEvents_Jet80_is80    ->Delete();

  h_NEvents_LowJets_read   ->Delete();
  h_NEvents_LowJets_skimCut->Delete();
  h_NEvents_LowJets_vzCut  ->Delete();
  h_NEvents_LowJets_PFMETfracCut	  ->Delete();
  h_NEvents_LowJets_is40   ->Delete();
  h_NEvents_LowJets_is60   ->Delete();

  h_NEvents_MC_read      ->Delete();
  h_NEvents_MC_skimCut	 ->Delete();
  h_NEvents_MC_vzCut     ->Delete();
  h_NEvents_MC_PFMETfracCut     ->Delete();
  
  fin->Close();
  fin2->Close();

  //// ============================================================================================ //
  //// ======================================= Jet Counts ========================================= //
  //// ============================================================================================ //
  //std::cout<<std::endl;
  //std::cout<<" ----- Jet80 PD Event Count Report ----- "<<std::endl;
  //std::cout<<"# Events ............................................ in Dataset: "<<  Nevts_jet80_read   <<", Fraction of total: "<< perc_Nevts_jet80_read   <<std::endl;
  //std::cout<<"# Events after HBHENoise, BeamScraping, and PV  quality  filters: "<<  Nevts_jet80_skimcut<<", Fraction of total: "<< perc_Nevts_jet80_skimcut<<std::endl;
  //std::cout<<"# Events after ..................................|vz|< 24 cm cut: "<<  Nevts_jet80_vzcut  <<", Fraction of total: "<< perc_Nevts_jet80_vzcut  <<std::endl;
  //std::cout<<"# Events after ............... passing trigger jet req for HLT80: "<<  Nevts_jet80_is80   <<", Fraction of total: "<< perc_Nevts_jet80_is80   <<std::endl;
  //
  //std::cout<<std::endl;
  //std::cout<<" ----- LowerJets PD Event Count Report ----- "<<std::endl;
  //std::cout<<"# Events ............................................ in Dataset: "<< Nevts_lowjets_read     <<", Fraction of total: "<< perc_Nevts_lowjets_read     <<std::endl;
  //std::cout<<"# Events after HBHENoise, BeamScraping, and PV  quality  filters: "<< Nevts_lowjets_skimcut  <<", Fraction of total: "<< perc_Nevts_lowjets_skimcut  <<std::endl;
  //std::cout<<"# Events after ..................................|vz|< 24 cm cut: "<< Nevts_lowjets_vzcut    <<", Fraction of total: "<< perc_Nevts_lowjets_vzcut    <<std::endl;
  //std::cout<<"# Events after ...... passing trigger jet req for HLT60 or HLT40: "<< Nevts_lowjets_is60or40 <<", Fraction of total: "<< perc_Nevts_lowjets_is60or40 <<std::endl;
  //std::cout<<"# Events ......................... passing trigger jet req for HLT40 only: "<<Nevts_lowjets_is40 <<", Fraction of total: "<< perc_Nevts_lowjets_is40 <<std::endl;
  //std::cout<<"# Events ......................... passing trigger jet req for HLT60 only: "<<Nevts_lowjets_is60 <<", Fraction of total: "<< perc_Nevts_lowjets_is60 <<std::endl;
  //
  //
  //std::cout<<std::endl;
  //
  //std::cout<<" ----- PYTHIA8 PD Event Count Report ----- "<<std::endl;
  //std::cout<<"# Events ............................................ in Dataset: "<< Nevts_mc_read     <<", Fraction of total: "<< perc_Nevts_mc_read     <<std::endl;
  //std::cout<<"# Events after HBHENoise, BeamScraping, and PV  quality  filters: "<< Nevts_mc_skimcut  <<", Fraction of total: "<< perc_Nevts_mc_skimcut  <<std::endl;
  //std::cout<<"# Events after ..................................|vz|< 24 cm cut: "<< Nevts_mc_vzcut    <<", Fraction of total: "<< perc_Nevts_mc_vzcut    <<std::endl;

  return;
}

void printpTbins(std::string outdir){
  std::cout<<"running printpTbins"<<std::endl;
  
  bool funcDebug=true;
  if(funcDebug)std::cout<<"in printpTbins"<<std::endl;
  if(funcDebug)std::cout<<"outdir is"<<outdir <<std::endl;
  
  const int netabins=NLO_UNFDIR_DATA_Nfiles;  
  if(netabins!=n_etabin_strs){
    printetabinerrormessage();
    return;
  }

  TH1D* spectra[netabins]={};		
  
  //first get the plots, clone + divide accordingly. binning should be set for me already, essentially
  for(int i=0; i<netabins; i++){
    std::string filepath= NLO_UNFDIR_DATA + NLO_UNFDIR_DATA_file_array[i] + ".root";//also has JEC systematics
    //std::cout<<"filepath="<<filepath<<std::endl;
    TFile* file=TFile::Open(( filepath).c_str(), "READ");

    //USE BINNING IN UNF DATA HIST
    spectra[i] =  (TH1D*)(
			  (
    			      (TH1D*)file->Get("Data_unf") 
    			      )->Clone( 
    				       ("Data_unf_ybin"+std::to_string(i)).c_str() 
    					)
    			     );    
    int nbins=spectra[i]->GetNbinsX();
    cout<<"printing pT bins for "<<const_ybin_strs[i]<<endl;
    
    for(int j=1; j<=nbins; j++){
      float ptbinlo_f=spectra[i]->GetBinLowEdge(j) ;
      float ptbinhi_f=ptbinlo_f+spectra[i]->GetBinWidth(j);
      cout<< (int)ptbinlo_f <<" - "<< (int)ptbinhi_f<<endl;          }    
    cout<<endl;
  }
  
  for(int i=0; i<netabins;i++)    spectra[i]->Delete();		
  return;
}












//--------------------------------------------------------------------------------------------------------------------------------
void  makeSMPInclJetXsec_fracstaterrData_ratios (std::string outdir, TFile* fout){
  std::cout<<"running makeSMPInclJetXsec_fracstaterrData_ratios"<<std::endl;
  
  const int netabins=NLO_UNFDIR_DATA_Nfiles;  
  if(netabins!=n_etabin_strs){
    printetabinerrormessage();
    return;
  }

  TH1D* spectra[netabins]={};		
	                                
  TH1D* measspectra[netabins]={};	
  TH1D* ratios_fracstatunc[netabins]={};
  
  //first get the plots, clone + divide accordingly. binning should be set for me already, essentially
  for(int i=0; i<netabins; i++){
    std::string filepath= NLO_UNFDIR_DATA + NLO_UNFDIR_DATA_file_array[i] + ".root";//also has JEC systematics
    //std::cout<<"filepath="<<filepath<<std::endl;
    TFile* file=TFile::Open(( filepath).c_str(), "READ");




    measspectra[i] =  (TH1D*)(
			      (
			       (TH1D*)file->Get("Data_meas_1GeVbins") 
			       )->Clone( 
					("Data_meas_1GeVbins_ybin"+std::to_string(i)).c_str() 
				     )
			      );


    
    ////USE BINNING IN HEADER
    //int nbins=-1;
    //double* ptbins=setBinning_etabin(i, &nbins);
    //measspectra[i]=(TH1D*)measspectra[i]->TH1::Rebin(nbins, ((std::string)measspectra[i]->GetName()+"_rebin").c_str(), ptbins );
    ////END USE BINNING IN HEADER

    //USE BINNING IN UNF DATA HIST
    spectra[i] =  (TH1D*)(
    			     (
    			      (TH1D*)file->Get("Data_unf") 
    			      )->Clone( 
    				       ("Data_unf_ybin"+std::to_string(i)).c_str() 
    					)
    			     );    
    //gonna have to get the binning here and rebin the stuff from the fastnlo files respectively
    int nbins=spectra[i]->GetNbinsX();
    std::vector<double> ptbins;
    for(int j=1;j<=nbins;j++)
      ptbins.push_back(spectra[i]->GetBinLowEdge(j));
    ptbins.push_back(spectra[i]->GetBinLowEdge(nbins)+spectra[i]->GetBinWidth(nbins));
    measspectra[i]=(TH1D*)measspectra[i]->TH1::Rebin(nbins, ((std::string)measspectra[i]->GetName()+"_rebin").c_str(), ((double*)ptbins.data()) );
    //END USE BINNING IN UNF DATA HIST

    //DATA STAT UNC
    ratios_fracstatunc[i]=(TH1D*) measspectra[i]->Clone(("Data_FracStatUnc_ratio_ybin"+std::to_string(i)).c_str());
    ratios_fracstatunc[i]->Reset("MICES");
    for(int bin=1;bin<=nbins;bin++){
      ratios_fracstatunc[i]->SetBinContent(bin, measspectra[i]->GetBinError(bin));
      ratios_fracstatunc[i]->SetBinError(bin, 1.e-30);
      measspectra[i]->SetBinError(bin,1.e-30);
    }
    ratios_fracstatunc[i]->Divide(measspectra[i]);
    ratios_fracstatunc[i]->SetMarkerSize(0.);  ratios_fracstatunc[i]->SetMarkerColor(kBlack);   ratios_fracstatunc[i]->SetMarkerStyle(kFullCircle);
    ratios_fracstatunc[i]->SetLineColor(kBlack);    
    
    ptbins.clear();
  }
  
  
  


  //TCanvas* canv=makeSMPRatioCanvas("fracstaterrData_SMPInclJetXsec_ratios");  
  TCanvas* canv=makeSMPRatioCanvas("fracstaterrData_SMPInclJetXsec_currentunfbinning_ratios");  
  //TCanvas* canv=makeSMPRatioCanvas("fracstaterrData_SMPInclJetXsec_lt50perc_nogaps_ratios");  
  //TCanvas* canv=makeSMPRatioCanvas("fracstaterrData_covmatErrors_SMPInclJetXsec_ratios");  

  TPaveText* SMPtitle=makePrelimPaveTextTitleRatio();


  for(int i=0; i<netabins; i++){
    
    canv->cd(i+1)->SetLogy(1);
    canv->cd(i+1)->SetLogx(1);
    canv->cd(i+1);

    ratios_fracstatunc[i]->SetTitle("");    
    ratios_fracstatunc[i]->SetMaximum(6.);
    ratios_fracstatunc[i]->GetYaxis()->CenterTitle(true);
    ratios_fracstatunc[i]->GetXaxis()->SetNoExponent(true);
    ratios_fracstatunc[i]->GetXaxis()->SetMoreLogLabels(true);
    setHistLabels(ratios_fracstatunc[i],"RECO Jet p_{T} [GeV]","Fractional Stat. Uncertainty");

    ratios_fracstatunc[i]->Draw("HIST");    
    
    float xlo=ratios_fracstatunc[i]->GetBinLowEdge(1);
    float xhi=
      ratios_fracstatunc[i]->GetBinLowEdge(ratios_fracstatunc[i]->GetNbinsX()) +   
      ratios_fracstatunc[i]->GetBinWidth(  ratios_fracstatunc[i]->GetNbinsX() );
    TLine* one     =makeTLine(xlo, 0.5 , xhi, 0.5);    
    one     ->Draw();       
   
    if(i==0){
      TLegend* leg=makeLegend(0.67, 0.71, 0.99, 0.92);
      leg->SetHeader("Merged SMP Inclusive Jet Bins","");
      leg->AddEntry(ratios_fracstatunc[i],"Fractional Stat. Unc.","l");
      leg->AddEntry(one, "50% Unc.","l");
      leg->Draw();
    }
    
    TPaveText* desc=NULL;
    if(i==0)desc=makePaveText( 0.15, 0.71, 0.45, 0.92);
    else    desc=makePaveText( 0.15, 0.78, 0.45, 0.92);
    desc->AddText(etabin_strs[i].c_str());
    std::string ptrange=ptcuts_lo;
    ptrange+=std::to_string( (int)xhi)+" GeV";
    desc->AddText(ptrange.c_str());
    if(i==0)desc->AddText(jettype.c_str());
    desc->Draw();
    
    SMPtitle->Draw();  
    
  }

  //makeSMPInclJetXsec_fracstaterrData_ratios
  saveCanv(outdir, canv, fout);
  for(int i=0; i<netabins;i++){
    spectra[i]->Delete();		
    measspectra[i]->Delete();	
    ratios_fracstatunc[i]->Delete();
  }

  return;
}






//--------------------------------------------------------------------------------------------------------------------------------
void  makeSMPInclJetXsec_covmatData_onePadOneEta (std::string outdir, TFile* fout){
  std::cout<<"running makeSMPInclJetXsec_covmatData_onePadOneEta"<<std::endl;
  
  const int netabins=PY8_UNFDIR_DATA_Nfiles;
  if(netabins!=n_etabin_strs){
    printetabinerrormessage();
    return;
  }
  
  TH2D* covmatrix[netabins]={};
  float zmax=-1.;
  float zmin=99999999999999999999.;
  
  //first get the plots, clone + divide accordingly. binning should be set for me already, essentially
  for(int i=0; i<netabins; i++){
    std::string filepath= PY8_UNFDIR_DATA + PY8_UNFDIR_DATA_file_array[i] + ".root";
    TFile* file=TFile::Open(( filepath).c_str(), "READ");
    
    covmatrix[i] =  (TH2D*)(
			(
			 (TH2D*)file->Get("Data_covmat_rebin") 
			 )->Clone( 
				  ("Data_covmat_ybin"+std::to_string(i)).c_str() 
				   )
			    );
    float th2max=covmatrix[i]->GetMaximum();
    float th2min=getnonzeromin((TH2*)covmatrix[i]);
    if(th2max>zmax)zmax=th2max;
    if(th2min<zmin)zmin=th2min;
    
  }
  
  //TCanvas* canv=makeSMPRatioCanvas("covmatData_onePadOneEta");
  TCanvas* canv=makeSMPTH2Canvas("covmatData_onePadOneEta");

  TPaveText* SMPtitle=makePrelimPaveTextTitleTH2();
  bool dologx=true, dology=true, dologz=true;
  for(int i=0; i<netabins; i++){
      
    canv->cd(i+1)->SetLogx(dologx);//maybe set to 1?
    canv->cd(i+1)->SetLogy(dology);//maybe set to 1?
    canv->cd(i+1)->SetLogz(dologz);//maybe set to 1?
    canv->cd(i+1);
    
    setHistLabels((TH2D*)covmatrix[i], "RECO Jet p_{T} [GeV]","RECO Jet p_{T} [GeV]");//,"Covariance [pb^{-2}]");
    
    if(dologx)covmatrix[i]->GetXaxis()->SetNoExponent(true);
    if(dologx)covmatrix[i]->GetXaxis()->SetMoreLogLabels(true);
    if(dology)covmatrix[i]->GetYaxis()->SetNoExponent(true);
    if(dology)covmatrix[i]->GetYaxis()->SetMoreLogLabels(true);
    covmatrix[i]->SetAxisRange(zmin*.9,zmax*1.1, "Z");
    covmatrix[i]->Draw("COLZ");    

    TPaveText* desc=makePaveText( 0.27, 0.93, 0.66, 1.00);
    std::string desctext=     "ak4PF Jets       "+etabin_strs[i];
    //etabin_strs[i]+"       ak4PF Jets";


    //desc->SetTextSize();
    //desc->AddText(etabin_strs[i].c_str());
    //desc->AddText(jettype.c_str());
    desc->AddText(desctext.c_str());
    desc->Draw();
    SMPtitle->Draw();
    
  }
  
  //makeSMPInclJetXsec_covmatData_onePadOneEta
  saveCanv(outdir, canv, fout);
  for(int i=0; i<netabins;i++){
    covmatrix[i]->Delete();
  }
  
  return;
}













//--------------------------------------------------------------------------------------------------------------------------------
void  makeSMPInclJetXsec_PY8vNLOunfdata_ratios (std::string outdir, TFile* fout){
  std::cout<<"running makeSMPInclJetXsec_PY8vNLOunfdata_ratio"<<std::endl;

  const int netabins=NLO_UNFDIR_DATA_Nfiles;  
  if(netabins!=n_etabin_strs){
    printetabinerrormessage();
    return;
  }
  
  TH1D* NLOunfspectra[netabins]={};
  TH1D* PY8unfspectra[netabins]={};
  //TH1D* smrPY8unfspectra[netabins]={};

  TH1D* ratios[netabins]={};
  
  //first get the plots, clone + divide accordingly. binning should be set for me already, essentially
  for(int i=0; i<netabins; i++){
    std::string NLOunffilepath= NLO_UNFDIR_DATA + NLO_UNFDIR_DATA_file_array[i] + ".root";//also has JEC systematics
    TFile* NLOunffile=TFile::Open(( NLOunffilepath).c_str(), "READ");
    
    //NLOSPECTRA ONLY, NO SYST
    NLOunfspectra[i] =  (TH1D*)(
			     (
			      (TH1D*)NLOunffile->Get("Data_unf") 
			      )->Clone( 
				       ("Data_unf_ybin"+std::to_string(i)).c_str() 
					)
			     );
    

    std::string PY8unffilepath= PY8_UNFDIR_DATA + PY8_UNFDIR_DATA_file_array[i] + ".root";
    TFile* PY8unffile=TFile::Open(( PY8unffilepath).c_str(), "READ");
    
    PY8unfspectra[i] =  (TH1D*)(
			     (
			      (TH1D*)PY8unffile->Get("Data_unf") 
			      )->Clone( 
				       ("Data_unf_ybin"+std::to_string(i)).c_str() 
					)
			     );

    ratios[i]=(TH1D*) PY8unfspectra[i]->Clone(("Data_NLO_PY8_ratio_ybin"+std::to_string(i)).c_str());
    ratios[i]->Divide(NLOunfspectra[i]);
    ratios[i]->SetMarkerSize(1.2);  ratios[i]->SetMarkerColor(kBlack);   ratios[i]->SetMarkerStyle(kFullCircle);
    ratios[i]->SetLineColor(kBlack);    
    ratios[i]->SetMinimum(0.6);    ratios[i]->SetMaximum(1.4);
    //setRatioHistLabels((TH1D*)ratios[i], "PY8 Unf. / NLO #otimes NP Unf.");
    setRatioHistLabels((TH1D*)ratios[i], "Ratio to NLO #otimes NP Unf. Data");
    for(int j=1; j<=ratios[i]->GetNbinsX(); j++)
      ratios[i]->SetBinError(j,1.e-30);

    //std::string smrPY8unffilepath= PY8_UNFDIR_DATA + PY8_UNFDIR_DATA_file_array[i] + ".root";
    //TFile* smrPY8unffile=TFile::Open(( smrPY8unffilepath).c_str(), "READ");
    //
    //smrPY8unfspectra[i] =  (TH1D*)(
    //				   (
    //				    (TH1D*)smrPY8unffile->Get("Data_unf") 
    //				    )->Clone( 
    //					     ("Data_unf_ybin"+std::to_string(i)).c_str() 
    //					      )
    //				   );
    //
    //ratios[i]=(TH1D*) PY8unfspectra[i]->Clone(("Data_NLO_PY8_ratio_ybin"+std::to_string(i)).c_str());
    //ratios[i]->Divide(NLOunfspectra[i]);
    //ratios[i]->SetMarkerSize(1.2);  ratios[i]->SetMarkerColor(kBlack);   ratios[i]->SetMarkerStyle(kFullCircle);
    //ratios[i]->SetLineColor(kBlack);    
    //ratios[i]->SetMinimum(0.6);    ratios[i]->SetMaximum(1.4);
    ////setRatioHistLabels((TH1D*)ratios[i], "PY8 Unf. / NLO #otimes NP Unf.");
    //setRatioHistLabels((TH1D*)ratios[i], "Ratio to NLO #otimes NP Unf. Data");
    

  }
  
  
 
  TCanvas* canv=makeSMPRatioCanvas("PY8vNLOunfdata_SMPInclJetXsec_ratio");
  
  TPaveText* SMPtitle=makePrelimPaveTextTitleRatio();
  for(int i=0; i<netabins; i++){


    float xlo=ratios[i]->GetBinLowEdge(1);
    float xhi=
      ratios[i]->GetBinLowEdge(ratios[i]->GetNbinsX()) +   
      ratios[i]->GetBinWidth(  ratios[i]->GetNbinsX() );
    TLine* one     =makeTLine(xlo, 1. , xhi, 1.);    
    
    canv->cd(i+1);
    ratios[i]->Draw("HIST E ][ ");    
    one     ->Draw();       
    ratios[i]->Draw("HIST E ][ SAME");    
    

    


    if(i==0){
      TLegend* leg=makeLegend(0.59, 0.63, 0.87, 0.87);
      //leg->AddEntry(ratios[i],"Full SIM PY8 unf.","lp");
      leg->AddEntry(ratios[i],"PY8 Unf. Data","lp");
      //leg->AddEntry(ratios2[i],"Smeared PY8 unf.","lp");
      leg->Draw();}
    
    TPaveText* desc=NULL;
    if(i==0)desc=makePaveText( 0.15, 0.67, 0.45, 0.88);
    else    desc=makePaveText( 0.15, 0.74, 0.45, 0.88);
    desc->AddText(etabin_strs[i].c_str());
    std::string ptrange=ptcuts_lo;
    ptrange+=std::to_string( (int)xhi)+" GeV";
    desc->AddText(ptrange.c_str());
    if(i==0)desc->AddText(jettype.c_str());
    desc->Draw();
        
    SMPtitle->Draw();  
    
  }
  
  
  //makeSMPInclJetXsec_PY8vNLOunfdata_ratios
  saveCanv(outdir, canv, fout);
  for(int i=0; i<netabins;i++){
    NLOunfspectra[i]->Delete();
    PY8unfspectra[i]->Delete()  ;
    ratios[i]->Delete();
  }
  
  return;
}


//--------------------------------------------------------------------------------------------------------------------------------
void  makeSMPInclJetXsec_PY8vNLOunfdata_modelsys_ratios (std::string outdir, TFile* fout, bool drawSysEnvelope, bool drawSymmEnvelope){
  std::cout<<"running makeSMPInclJetXsec_PY8vNLOunfdata_modelsys_ratios"<<std::endl;
  
  const int netabins=NLO_UNFDIR_DATA_Nfiles;  
  if(netabins!=n_etabin_strs){
    printetabinerrormessage();
    return;
  }
  
  TH1D* NLOunfspectra[netabins]={};//for the main result i want to construct and unc band for. 
  TH1D* altNLOunfspectra[NLO_UNFDIR_DATA_DIFFMODELS_Nfiles][netabins]={};//alt toy nlo unfolding that i want to use to define the deviation from the primary measurement
  TH1D* PY8unfspectra[netabins]={}; //another alt unfolding for comparison.
  //TH1D* smrPY8unfspectra[netabins]={};
  
  TH1D* altratios[NLO_UNFDIR_DATA_DIFFMODELS_Nfiles][netabins]={}; //what gets cloned into this + other ratios first should be my "favorite" spectra (gauss core unfold w/ ansatz fits)
  TH1D* PY8ratios[netabins]={};                                    //what gets cloned into this + other ratios first should be my "favorite" spectra (gauss core unfold w/ ansatz fits)
  
  //std::vector<std::vector<double>> sysup={{}, {}, {}, {} };//will record max dev. from 1 in pos direction.
  std::vector<std::vector<double>> sysup;//will record max dev. from 1 in pos direction.
  TH1D* modsysup[netabins]={};
  std::vector<std::vector<double>> sysdown;//will record max dev from 1 in neg direction.
  TH1D* modsysdown[netabins]={};
  
  //first get the plots, clone + divide accordingly. binning should be set for me already, essentially
  for(int i=0; i<netabins; i++){
    std::string NLOunffilepath= NLO_UNFDIR_DATA + NLO_UNFDIR_DATA_file_array[i] + ".root";//also has JEC systematics
    TFile* NLOunffile=TFile::Open(( NLOunffilepath).c_str(), "READ");
    
    //NLOSPECTRA ONLY, NO SYST
    NLOunfspectra[i] =  (TH1D*)(
				(
				 (TH1D*)NLOunffile->Get("Data_unf") 
				 )->Clone( 
					  ("Data_unf_ybin"+std::to_string(i)).c_str() 
					   )
				);
    
    
    std::string PY8unffilepath= PY8_UNFDIR_DATA + PY8_UNFDIR_DATA_file_array[i] + ".root";
    TFile* PY8unffile=TFile::Open(( PY8unffilepath).c_str(), "READ");
    
    PY8unfspectra[i] =  (TH1D*)(
				(
				 (TH1D*)PY8unffile->Get("Data_unf") 
				 )->Clone( 
					  ("Data_PY8unf_ybin"+std::to_string(i)).c_str() 
					   )
				);

    if(drawSysEnvelope){
      modsysup[i]  =(TH1D*)PY8unfspectra[i]->Clone((  "modsysup_ybin"+std::to_string(i)).c_str());
      modsysup[i]->Reset("MICES");
      modsysdown[i]=(TH1D*)PY8unfspectra[i]->Clone(("modsysdown_ybin"+std::to_string(i)).c_str());
      modsysdown[i]->Reset("MICES");
    }
      
    
    
    
    PY8ratios[i]=(TH1D*) PY8unfspectra[i]->Clone(("Data_NLO_PY8_ratio_ybin"+std::to_string(i)).c_str());
    PY8ratios[i]->Divide(NLOunfspectra[i]);
    PY8ratios[i]->SetMarkerSize(0);  PY8ratios[i]->SetMarkerColor(kBlack);   PY8ratios[i]->SetMarkerStyle(kFullCircle);
    PY8ratios[i]->SetLineColor(kGreen);   //   PY8ratios[i]->SetLineColor(kBlack);    
    PY8ratios[i]->SetMinimum(0.85);    PY8ratios[i]->SetMaximum(1.25);
    //setRatioHistLabels((TH1D*)PY8ratios[i], "PY8 Unf. / NLO #otimes NP Unf.");
    setRatioHistLabels((TH1D*)PY8ratios[i], "Ratio to NLO #otimes NP Unf. Data");
    
    //here we will set the bin error to negligible (to draw the markers, and b.c. the unc are correlated and so don't make inherentsense)
    std::vector<double> up  ={};    std::vector<double> down={};
    for(int j=1; j<=PY8ratios[i]->GetNbinsX(); j++){
      PY8ratios[i]->SetBinError(j,1.e-30);
      if(drawSysEnvelope){
	double dev=PY8ratios[i]->GetBinContent(j)-1.;
	if(dev>0.){
	  up.push_back(dev);
	  down.push_back(0.);
	  //std::cout<<"POS dev=="<<dev<<" > 0 push_back #"<<j<<"/"<<PY8ratios[i]->GetNbinsX()<<"\n";
	  std::cout<<"POS dev=="<<up.at(j-1)<<" > 0 push_back #"<<j<<"/"<<PY8ratios[i]->GetNbinsX()<<"\n";
	}
	else{//if dev is neg or zero
	  up.push_back(0.);
	  down.push_back(dev);
	  std::cout<<"NEG dev=="<<down.at(j-1)<<" < 0 push_back #"<<j<<"/"<<PY8ratios[i]->GetNbinsX()<<"\n";
	}
	
      }
      
    }//loop over PY8 ratio bin contents/errors
    
    
    
    
    for(int j=0;j<NLO_UNFDIR_DATA_DIFFMODELS_Nfiles;j++){
      std::string altNLOunffilepath= NLO_UNFDIR_DATA + NLO_UNFDIR_DATA_DIFFMODELS_file_array[j] + ETABIN_TAG_array[i] + ".root";//also has JEC systematics
      TFile* altNLOunffile=TFile::Open(( altNLOunffilepath).c_str(), "READ");
      
      //NLOSPECTRA ONLY, NO SYST
      altNLOunfspectra[j][i] =  (TH1D*)(
					(
					 (TH1D*)altNLOunffile->Get("Data_unf") 
					 )->Clone( 
						  ("Data_unf_alt"+std::to_string(j)+"_ybin"+std::to_string(i)).c_str() 
						   )
					);
      
      
      altratios[j][i]=(TH1D*) altNLOunfspectra[j][i]->Clone(("Data_NLO_v_alt"+std::to_string(j)+"_ratio_ybin"+std::to_string(i)).c_str());
      altratios[j][i]->Divide(NLOunfspectra[i]);
      altratios[j][i]->SetMarkerSize(0);
      if(     j==0)
	altratios[j][i]->SetLineColor(kRed);
      else if(j==1)
	altratios[j][i]->SetLineColor(kBlue);
      else
	altratios[j][i]->SetLineColor(kViolet);
      
      for(int k=1; k<=altratios[j][i]->GetNbinsX(); k++){
	//std::cout<<"k=="<<k<<
	altratios[j][i]->SetBinError(k,1.e-30);
	if(drawSysEnvelope){
	  double dev= altratios[j][i]->GetBinContent(k)-1.;
	  if(dev>0.){
	    if( up.at(k-1)<dev) {
	      std::cout<<"found a bigger +dev from 1 than up["<<k-1<<"]="<<up[k-1]<<"\n";
	      up[k-1]=dev;
	      std::cout<<"Now up["<<k-1<<"]="<<up[k-1]<<"\n";	      
	    }
	  }
	  else {//if dev is neg or zero
	    if(down.at(k-1)>dev){
	      std::cout<<"found a bigger -dev from 1 than down["<<k-1<<"]="<<down[k-1]<<"\n";
	      down[k-1]=dev;
	      std::cout<<"Now down["<<k-1<<"]="<<down[k-1]<<"\n";	      	      
	    }
	  }
	}//end draw sys envelope
      }//loop over bin contents/errors
      
    }//end loop over diff toy model unfoldings
  
  
  
    //return;
  
    
    if(drawSysEnvelope)
      {
	if(drawSymmEnvelope)
	  {
	    for(unsigned int j=0;j<up.size();j++)
	      {
		if(      fabs(  up.at(j))>fabs(down.at(j)) )
		  {
		    down[j]=-1.*  up.at(j);
		  }
		else if( fabs(down.at(j))>fabs(  up.at(j)) )
		  {
		    up[j]  =-1.*down.at(j);
		  }	  
	      }
	  }
	sysup.push_back(up);
	sysdown.push_back(down);
      }
    
  }

    
 
  
  
  std::string canvname="PY8vNLOunfdata_SMPInclJetXsec_";
  if(drawSysEnvelope){
    if(drawSymmEnvelope)
      canvname+="Symmodelsys_";
    else
      canvname+="Asymmodelsys_";    
  }
  canvname+="ratio";
  
  //TCanvas* canv=makeSMPRatioCanvas("PY8vNLOunfdata_SMPInclJetXsec_modelsys_ratio");
  TCanvas* canv=makeSMPRatioCanvas(canvname);
  
  TPaveText* SMPtitle=makePrelimPaveTextTitleRatio();
  for(int i=0; i<netabins; i++){
    
    if(drawSysEnvelope){
      for(int j=1; j<=modsysup[i]->GetNbinsX();j++){
	modsysup[i]->SetBinContent(j,1.+sysup[i][j-1]);
	modsysup[i]->SetBinError(  j,1.e-30     );
	modsysdown[i]->SetBinContent(j,1.+sysdown[i][j-1]);
	modsysdown[i]->SetBinError(  j,1.e-30);
      }
    }
    
    
    float xlo=PY8ratios[i]->GetBinLowEdge(1);
    float xhi=
      PY8ratios[i]->GetBinLowEdge(PY8ratios[i]->GetNbinsX()) +   
      PY8ratios[i]->GetBinWidth(  PY8ratios[i]->GetNbinsX() );
    TLine* one     =makeTLine(xlo, 1. , xhi, 1.);    
    
    TLegend* leg=NULL;
    if(i==0){
      leg=makeLegend(0.55, 0.60, 0.99, 0.90);
      //leg->AddEntry(ratios[i],"Full SIM PY8 unf.","lp");
      leg->AddEntry(PY8ratios[i],"PY8 Unf. Data","lp");
      //leg->AddEntry(ratios2[i],"Smeared PY8 unf.","lp");
    }
    
    canv->cd(i+1);
    PY8ratios[i]->Draw("HIST ][ ");    
    one     ->Draw();       
    for(int j=0;j<NLO_UNFDIR_DATA_DIFFMODELS_Nfiles;j++){
      altratios[j][i]->Draw("HIST ][ SAME");
      if(i==0){
	if(j==0)
	  leg->AddEntry(altratios[j][i],"NLO Unf. Data (Gauss Smear w/ DSCB)","lp");
	else if(j==1)
	  leg->AddEntry(altratios[j][i],"NLO Unf. Data (Smear with DSCB TF1s)","lp");
	else
	  leg->AddEntry(altratios[j][i],"alt toy model","lp");
      }
    }
    
    if(drawSysEnvelope){
      modsysup[i]->SetMarkerSize(0);      modsysup[i]->SetMarkerColor(kBlack);      modsysup[i]->SetMarkerStyle(26);
      modsysup[i]->SetLineColor(kBlack);
      //modsysup[i]->SetLineStyle(10);

      modsysdown[i]->SetMarkerSize(0);      modsysdown[i]->SetMarkerColor(kBlack);      modsysdown[i]->SetMarkerStyle(32);
      modsysdown[i]->SetLineColor(kBlack);
      //modsysdown[i]->SetLineStyle(10);
      if(i==0){
	if(drawSymmEnvelope)
	  leg->AddEntry(modsysup[i],"Model Sys. Unc. (Symmetric)","l");
	else
	  leg->AddEntry(modsysup[i],"Model Sys. Unc. (Asymmetric)","l");	
      }
      modsysup[i]->Draw("HIST ][ SAME");
      modsysdown[i]->Draw("HIST ][ SAME");
    }
    
    //PY8ratios[i]->Draw("HIST E ][ SAME");    


    if(i==0){
      leg->Draw();
    }
    
    TPaveText* desc=NULL;
    if(i==0)desc=makePaveText( 0.15, 0.67, 0.45, 0.88);
    else    desc=makePaveText( 0.15, 0.74, 0.45, 0.88);
    desc->AddText(etabin_strs[i].c_str());
    std::string ptrange=ptcuts_lo;
    ptrange+=std::to_string( (int)xhi)+" GeV";
    desc->AddText(ptrange.c_str());
    if(i==0)desc->AddText(jettype.c_str());
    desc->Draw();
        
    SMPtitle->Draw();  
    
  }
  
  
  //makeSMPInclJetXsec_PY8vNLOunfdata_modelsys_ratios
  saveCanv(outdir, canv, fout);
  for(int i=0; i<netabins;i++){
    NLOunfspectra[i]->Delete();
    PY8unfspectra[i]->Delete()  ;
    PY8ratios[i]->Delete();
    for(int j=0;j<NLO_UNFDIR_DATA_DIFFMODELS_Nfiles;j++){
      altNLOunfspectra[j][i]->Delete();
      altratios[j][i]->Delete();
    }

  }
  
  return;
}


