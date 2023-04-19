void  makeSMPInclJetXsec_PY8JERmu_onePadOneEta (   std::string outdir, TFile* fout=NULL, bool drawFit=false, bool drawFitUnc=false,bool drawJohnDSCB=false, bool drawPatDSCB=false);
void  makeSMPInclJetXsec_PY8JERsigma_onePadOneEta (std::string outdir, TFile* fout=NULL, bool drawFit=false, bool drawFitUnc=false, bool drawJohnDSCB=false, bool drawPatDSCB=false);
void  makeSMPInclJetXsec_PY8JERkLR_onePadOneEta (std::string outdir, TFile* fout=NULL, bool drawJohnDSCB=false    );
void  makeSMPInclJetXsec_PY8JERalphaLR_onePadOneEta (std::string outdir, TFile* fout=NULL, bool drawJohnDSCB=false);

void  makeSMPInclJetXsec_NLOSmearingFits_onePadAllEta (std::string outdir, TFile* fout=NULL);
void  makeSMPInclJetXsec_NLOSmearingFitswSFs_onePadAllEta (std::string outdir, TFile* fout=NULL);
void  makeSMPInclJetXsec_NLOSmearingFits_onePadOneEta (std::string outdir, TFile* fout=NULL);
void  makeSMPInclJetXsec_NLOSmearingFitswSFs_onePadOneEta (std::string outdir, TFile* fout=NULL, bool drawJohnDSCB=false);







//--------------------------------------------------------------------------------------------------------------------------------
void  makeSMPInclJetXsec_PY8JERmu_onePadOneEta (std::string outdir, TFile* fout, bool drawFit, bool drawFitUnc, bool drawJohnDSCB, bool drawPatDSCB){
  std::cout<<"running makeSMPInclJetXsec_PY8JERmu_onePadOneEta"<<std::endl;  

  const int netabins=JERDIR_MC_Nfiles;  
  if(netabins!=n_etabin_strs){
    printetabinerrormessage();
    return;
  }
  
  TH1D* mu_datapoints[netabins]={};
  TF1* mu_fit[netabins]={};
  TF1* mu_fit_sysup[netabins]={};
  TF1* mu_fit_sysdown[netabins]={};
  TF1* mu_fit_10xsysup[netabins]={};
  TF1* mu_fit_10xsysdown[netabins]={};

  TH1D* mu_datapoints_johndscb[netabins]={};
  TF1* mu_fit_dscb[netabins]={};

  TH1D* mu_datapoints_patdscb[netabins]={};
  TF1* mu_fit_patdscb[netabins]={};
  
  float lineptmax[netabins]={0.};
  float ptmax[netabins]={0.};
  float lineptmin[netabins]={0.};
  float ptmin[netabins]={0.};
  
  //first get the plots, clone + divide accordingly. binning should be set for me already, essentially
  for(int i=0; i<netabins; i++){
    std::string filepath= JERDIR_MC + JERDIR_MC_file_array[i] + ".root";
    TFile* file=TFile::Open(( filepath).c_str(), "READ");

    std::string filepath2=JERDIR_MC + JERDIR2_MC_file_array[i] + ".root";
    TFile* patfile=NULL;
    if(drawPatDSCB||drawJohnDSCB) 
      patfile=TFile::Open( (filepath2).c_str(),"READ");    
    
    if(drawFit){
      mu_fit[i]=(TF1*)(
		       (
			(TF1*)file->Get("MuFit_f")
			)->Clone(
				 ("MuFit_f_ybin"+std::to_string(i)).c_str()
				 )
		       );

      mu_fit[i]->SetLineColorAlpha( mu_fit[i]->GetLineColor(), 0.6);
      mu_fit[i]->SetLineWidth(1);

      if(drawJohnDSCB&&false){
	mu_fit_dscb[i]=(TF1*)(
			      (
			       (TF1*)patfile->Get("")
			       )->Clone(
					("MuFit_dscb_ybin"+std::to_string(i)).c_str()
					)
			      );
      }

      if(drawPatDSCB&&false){
	mu_fit_patdscb[i]=(TF1*)(
				 (
				  (TF1*)patfile->Get("")
				  )->Clone(
					   ("MuFit_patdscb_ybin"+std::to_string(i)).c_str()
					   )
				 );
	
      }
      
      if(drawFitUnc){
	mu_fit_sysup[i]=(TF1*)(
			       (
				(TF1*)file->Get("MuFit_f_sysup")
				)->Clone(
					 ("MuFit_f_sysup_ybin"+std::to_string(i)).c_str()
					 )
			       );			
	mu_fit_sysdown[i]=(TF1*)(
			       (
				(TF1*)file->Get("MuFit_f_sysdown")
				)->Clone(
					 ("MuFit_f_sysdown_ybin"+std::to_string(i)).c_str()
					 )
			       );			

	mu_fit_10xsysup[i]=(TF1*)(
			       (
				(TF1*)file->Get("MuFit_f_10xsysup")
				)->Clone(
					 ("MuFit_f_10xsysup_ybin"+std::to_string(i)).c_str()
					 )
			       );			
	mu_fit_10xsysdown[i]=(TF1*)(
			       (
				(TF1*)file->Get("MuFit_f_10xsysdown")
				)->Clone(
					 ("MuFit_f_10xsysdown_ybin"+std::to_string(i)).c_str()
					 )
				  );			
      }//end drawFitUnc
      
    }//end drawFit
    
    mu_datapoints[i] =  (TH1D*)(
				(
				 (TH1D*)file->Get("hMean_fit") 
				 )->Clone( 
					  ("hMean_datapoints_ybin"+std::to_string(i)).c_str() 
					   )
				);        
    mu_datapoints[i]->SetMarkerStyle(kOpenCircle);
    mu_datapoints[i]->SetMarkerColor(mu_datapoints[i]->GetLineColor());
    mu_datapoints[i]->SetMarkerSize(1);
    
    if(drawJohnDSCB){
      mu_datapoints_johndscb[i] =  (TH1D*)(
				       (
					(TH1D*)patfile->Get("hJohnDSCBMean") 
					)->Clone( 
						 ("hJohnDSCBMean_clone_ybin"+std::to_string(i)).c_str() 
						  )
				       );        
      mu_datapoints_johndscb[i]->SetMarkerStyle(kOpenSquare);
      mu_datapoints_johndscb[i]->SetMarkerSize(1);
      for(int p=1;p<=mu_datapoints_johndscb[i]->GetNbinsX();p++)
	mu_datapoints_johndscb[i]->SetBinError(p, 1.e-20);//so the marker actually draws.
    }
    
    if(drawPatDSCB){
      mu_datapoints_patdscb[i] =  (TH1D*)(
				       (
					(TH1D*)patfile->Get("hDSCBMean") 
					)->Clone( 
						 ("hJohnDSCBMean_clone_ybin"+std::to_string(i)).c_str() 
						  )
				       );        
      mu_datapoints_patdscb[i]->SetMarkerStyle(kOpenCross);
    }

    //std::cout<<"hello bin #"<<i<<std::endl;
    mu_datapoints[i]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
    //mu_datapoints[i]->GetYaxis()->SetTitle("Gaussian Core Fit #mu");
    mu_datapoints[i]->GetYaxis()->SetTitle("#mu Extracted From Fit");
    mu_datapoints[i]->GetYaxis()->SetTitleOffset(0.95);
    //    std::cout<<"mu_datapoints["<<i<<"]->GetYaxis()->GetTitleOffset()="<<mu_datapoints[i]->GetYaxis()->GetTitleOffset()<<std::endl;   
    //    std::cout<<"mu_datapoints["<<i<<"]->GetYaxis()->GetTitleSize()  ="<<mu_datapoints[i]->GetYaxis()->GetTitleSize()<<std::endl;
    //std::cout<<"mu_datapoints["<<i<<"]->GetYaxis()->GetLabelOffset()="<<mu_datapoints[i]->GetYaxis()->GetLabelOffset()<<std::endl;
    //std::cout<<"mu_datapoints["<<i<<"]->GetYaxis()->GetLabelSize()  ="<<mu_datapoints[i]->GetYaxis()->GetLabelSize()<<std::endl;    
    //std::cout<<"mu_datapoints["<<i<<"]->GetXaxis()->GetTitleOffset()="<<mu_datapoints[i]->GetXaxis()->GetTitleOffset()<<std::endl;
    //std::cout<<"mu_datapoints["<<i<<"]->GetXaxis()->GetTitleSize()  ="<<mu_datapoints[i]->GetXaxis()->GetTitleSize()<<std::endl;
    //std::cout<<"mu_datapoints["<<i<<"]->GetXaxis()->GetLabelOffset()="<<mu_datapoints[i]->GetXaxis()->GetLabelOffset()<<std::endl;
    //std::cout<<"mu_datapoints["<<i<<"]->GetXaxis()->GetLabelSize()  ="<<mu_datapoints[i]->GetXaxis()->GetLabelSize()<<std::endl;
    mu_datapoints[i]->GetYaxis()->CenterTitle(true);
    bool minimumset=false;
    for(int j=1; j<=mu_datapoints[i]->GetNbinsX();j++){
      //std::cout<<"hello pt bin #"<<i<<std::endl;
      float lowedge  =mu_datapoints[i]->GetBinLowEdge(j);          //std::cout<<"lowedge ="<<lowedge <<std::endl;
      float highedge =mu_datapoints[i]->GetBinWidth(j)+lowedge;    //std::cout<<"highedge="<<highedge<<std::endl;
      float mu       =mu_datapoints[i]->GetBinContent(j);	      //std::cout<<"mu      ="<<mu      <<std::endl;
      float muerr    =mu_datapoints[i]->GetBinError(j);      	      //std::cout<<"muerr   ="<<muerr   <<std::endl;
      if(!(mu>0.))continue;
      //if(j<8)continue;//cause i expanded the pT range and now this macro plots everything. i don't want that!
      if(  j<8)	continue;//cause i expanded the pT range and now this macro plots everything. i don't want that!
      if(!minimumset){
	ptmin[i]=lowedge;
	//lineptmin[i]=lowedge-(highedge-lowedge);//...why do i do this again? guess we will see...
	lineptmin[i]=lowedge;
	std::cout<<"MINIMUM SET: ptmin[i]="<<ptmin[i]<<std::endl;
	minimumset=true;
      }
      if(ptmax[i]<highedge){
	ptmax[i]=highedge;
	lineptmax[i]=highedge+mu_datapoints[i]->GetBinWidth(j+1);//...why do i do this again? guess we will see...	  
      }
    }
    
    if(drawPatDSCB)
      for(int j=1; j<=mu_datapoints_patdscb[i]->GetNbinsX();j++)
	if((mu_datapoints_patdscb[i]->GetBinCenter(j)>ptmax[i]) ||
	   (mu_datapoints_patdscb[i]->GetBinCenter(j)<ptmin[i]) 	 )//FOR COSMETICS (points drawing off plot... why?!)
	  mu_datapoints_patdscb[i]->SetBinContent(j,1.e+20);
  }
  
  
  std::string canvname="PY8JERmu_onePadOneEta";
  if(!drawFit)canvname+="_datapointsOnly";
  TCanvas* canv=makeSMPRatioCanvas(canvname);

  TPaveText* SMPtitle=makeSimPaveTextTitleRatio();
  bool dologx=true;
  for(int i=0; i<netabins; i++){
    
    canv->cd(i+1)->SetLogx(dologx);//maybe set to 1?
    canv->cd(i+1)->SetLogy(0);
    canv->cd(i+1)->SetTicky(1);
    canv->cd(i+1);
    
    float xlo=ptmin[i];
    float xhi=ptmax[i];    
    float linexlo=lineptmin[i];    
    float linexhi=lineptmax[i];    
    
    mu_datapoints[i]->SetTitle("");
    if(dologx)mu_datapoints[i]->GetXaxis()->SetNoExponent(true);
    if(dologx)mu_datapoints[i]->GetXaxis()->SetMoreLogLabels(true);
    mu_datapoints[i]->SetAxisRange(linexlo, xhi, "X");
    mu_datapoints[i]->SetAxisRange(0.98,1.03, "Y");
    mu_datapoints[i]->Draw("HIST E ][");    

    if(drawFit){
      mu_fit[i]->Draw("L SAME");
      if(drawFitUnc){
	mu_fit_sysup[i]->Draw("L SAME");
	mu_fit_sysdown[i]->Draw("L SAME");
	mu_fit_10xsysup[i]->Draw("L SAME");
	mu_fit_10xsysdown[i]->Draw("L SAME");	
      }
    }
    
    //TLine* onep01     =makeTLine( linexlo, 1.01 ,linexhi , 1.01);    onep01     ->Draw();       
    TLine* one        =makeTLine( linexlo, 1.00 ,linexhi , 1.00);    one     ->Draw();       
    //TLine* p99        =makeTLine( linexlo, 0.99 ,linexhi , 0.99);    p99     ->Draw();       
    
    mu_datapoints[i]->Draw("HIST E ][ SAME");    
    if(drawJohnDSCB){
      mu_datapoints_johndscb[i]->Draw("HIST E ][ SAME");
      if(drawFit&&false)
	mu_fit_dscb[i]->Draw("L SAME");      
    }
    if(drawPatDSCB){
      mu_datapoints_patdscb[i]->Draw("HIST E ][ SAME");
      if(drawFit&&false)
	mu_fit_patdscb[i]->Draw("L SAME");      
    }
    
    
    if(i==1){
      TLegend* leg=makeLegend(0.65, 0.55, 1.00, 0.90);
      leg->AddEntry((TObject*)0,"PYTHIA8","");
      leg->AddEntry(mu_datapoints[i],"Gaussian Core Fits","lpe");
      if(drawFit){
	leg->AddEntry(mu_fit[i],"Fit of Gauss Core #mu", "l");
	if(drawFitUnc){
	  leg->AddEntry(mu_fit_sysup[i],"Fit Unc.", "l");	  
	  leg->AddEntry(mu_fit_10xsysup[i],"10#times(Fit Unc.)", "l");	  
	}
      }
      if(drawJohnDSCB){
	leg->AddEntry(mu_datapoints_johndscb[i],"from DSCB+Extra Fits","lpe");
	if(drawFit&&false)
	  leg->AddEntry(mu_fit_dscb[i],"Fit of DSCB+Extra #mu","l");	
      }  
      if(drawPatDSCB){
	leg->AddEntry(mu_datapoints_patdscb[i],"from DSCB Fits","lpe");
	if(drawFit&&false)
	  leg->AddEntry(mu_fit_patdscb[i],"Fit of DSCB #mu","l");	
      }  
      leg->Draw();
    }
    
    TPaveText* desc=NULL;
    if(i==0)desc=makePaveText( 0.15, 0.67, 0.45, 0.88);
    else    desc=makePaveText( 0.15, 0.74, 0.45, 0.88);
    desc->AddText(etabin_strs[i].c_str());
    std::string ptrange=std::to_string( (int)xlo) +
                            " GeV < Jet p_{T} < " +
                        std::to_string( (int)xhi) + " GeV";
    desc->AddText(ptrange.c_str());
    if(i==0)desc->AddText(jettype.c_str());
    desc->Draw();
    
    SMPtitle->Draw();  
  }
  
  //makeSMPInclJetXsec_PY8JERmu_onePadOneEta
  saveCanv(outdir, canv, fout);
  for(int i=0; i<netabins;i++){
    mu_datapoints[i]->Delete();
    if(drawJohnDSCB)mu_datapoints_johndscb[i]->Delete();
    if(drawPatDSCB)mu_datapoints_patdscb[i]->Delete();
  }
  //canv->Delete();

  return;
}

//--------------------------------------------------------------------------------------------------------------------------------
void  makeSMPInclJetXsec_PY8JERsigma_onePadOneEta (std::string outdir, TFile* fout, bool drawFit, bool drawFitUnc, bool drawJohnDSCB, bool drawPatDSCB){
  std::cout<<"running makeSMPInclJetXsec_PY8JERsigma_onePadOneEta"<<std::endl;

  const int netabins=JERDIR_MC_Nfiles;  
  if(netabins!=n_etabin_strs){
    printetabinerrormessage();
    return;
  }
  
  TH1D* sigma_datapoints[netabins]={};
  TF1* sigma_fit[netabins]={};
  TF1* sigma_fit_sysup[netabins]={};
  TF1* sigma_fit_sysdown[netabins]={};
  TF1* sigma_fit_10xsysup[netabins]={};
  TF1* sigma_fit_10xsysdown[netabins]={};

  TH1D* sigma_datapoints_johndscb[netabins]={};
  TF1* sigma_fit_dscb[netabins]={};

  TH1D* sigma_datapoints_patdscb[netabins]={};
  TF1* sigma_fit_patdscb[netabins]={};

  float lineptmax[netabins]={0.};
  float ptmax[netabins]={0.};
  float lineptmin[netabins]={0.};
  float ptmin[netabins]={0.};
  
  //first get the plots, clone + divide accordingly. binning should be set for me already, essentially
  for(int i=0; i<netabins; i++){
    std::string filepath= JERDIR_MC + JERDIR_MC_file_array[i] + ".root";
    TFile* file=TFile::Open(( filepath).c_str(), "READ");

    std::string sigma_fit_name="SigmaFit_f";
    //if(drawJohnDSCB)sigma_fit_name="SigmaFit_dscb";
    std::string sigma_datapoints_name="hSigma_fit";
    //if(drawJohnDSCB)sigma_datapoints_name="hJohnDSCBSigma_clone";

    std::string filepath2=JERDIR_MC + JERDIR2_MC_file_array[i] + ".root";
    TFile* patfile=NULL;
    if(drawPatDSCB||drawJohnDSCB) 
      patfile=TFile::Open( (filepath2).c_str(),"READ");
    
    if(drawFit){
      sigma_fit[i]=(TF1*)(
			  (
			   (TF1*)file->Get(sigma_fit_name.c_str())
			   )->Clone(
				    (sigma_fit_name+"_ybin"+std::to_string(i)).c_str()
				    )
			  );
      sigma_fit[i]->SetLineColorAlpha( sigma_fit[i]->GetLineColor(), 0.6);
      sigma_fit[i]->SetLineWidth(2);

      if(drawFitUnc){
	//file->ls();
	//return;
	sigma_fit_sysup[i]=(TF1*)(
				  (
				   (TF1*)file->Get("SigmaFit_sysup_f")
				   )->Clone(
					    ("SigmaFit_sysup_f_ybin"+std::to_string(i)).c_str()
					    )
				  );
	
	sigma_fit_sysdown[i]=(TF1*)(
				    (
				     (TF1*)file->Get("SigmaFit_sysdown_f")
				     )->Clone(
					      ("SigmaFit_sysdown_f_ybin"+std::to_string(i)).c_str()
					      )
				    );

	sigma_fit_10xsysup[i]=(TF1*)(
				  (
				   (TF1*)file->Get("SigmaFit_10xsysup_f")
				   )->Clone(
					    ("SigmaFit_10xsysup_f_ybin"+std::to_string(i)).c_str()
					    )
				  );
	
	sigma_fit_10xsysdown[i]=(TF1*)(
				    (
				     (TF1*)file->Get("SigmaFit_10xsysdown_f")
				     )->Clone(
					      ("SigmaFit_10xsysdown_f_ybin"+std::to_string(i)).c_str()
					      )
				    );
      }
      
      if(drawJohnDSCB&&false)
	sigma_fit_dscb[i]=(TF1*)(
				 (
				  (TF1*)file->Get("SigmaFit_dscb")
				  )->Clone(
					   ("SigmaFit_dscb_ybin"+std::to_string(i)).c_str()
					   )
				 );
      
      if(drawPatDSCB&&false)
	sigma_fit_patdscb[i]=(TF1*)(
				    (
				     (TF1*)file->Get("")
				     )->Clone(
					      ("SigmaFit_patdscb_ybin"+std::to_string(i)).c_str()
					      )
				    );
      
      
      
    }

    sigma_datapoints[i] =  (TH1D*)(
				   (
				    (TH1D*)file->Get(sigma_datapoints_name.c_str()) 
				    )->Clone( 
					     (sigma_datapoints_name+"_ybin"+std::to_string(i)).c_str() 
					      )
				   );        
    sigma_datapoints[i]->SetMarkerStyle(kOpenCircle);
    sigma_datapoints[i]->SetMarkerColor(sigma_datapoints[i]->GetLineColor());
    sigma_datapoints[i]->SetMarkerSize(1);

    if(drawJohnDSCB){
      sigma_datapoints_johndscb[i] =  (TH1D*)(
					  (
					   (TH1D*)patfile->Get("hJohnDSCBSigma") 
					   )->Clone( 
						    ("hJohnDSCBSigma_clone_ybin"+std::to_string(i)).c_str() 
						     )
					  );        
      sigma_datapoints_johndscb[i]->SetMarkerStyle(kOpenSquare);
      sigma_datapoints_johndscb[i]->SetMarkerSize(1);
      for(int p=1;p<=sigma_datapoints_johndscb[i]->GetNbinsX();p++)
	sigma_datapoints_johndscb[i]->SetBinError(p, 1.e-20);//so the marker actually draws.
    }
    
    
    
    if(drawPatDSCB){
      sigma_datapoints_patdscb[i] =  (TH1D*)(
					     (
					      (TH1D*)patfile->Get("hDSCBSigma") 
					      )->Clone( 
						       ("hPatDSCBSigma_clone_ybin"+std::to_string(i)).c_str() 
							)
					     );        
      sigma_datapoints_patdscb[i]->SetMarkerStyle(kOpenCross);
    }
    
    
    //std::cout<<"hello bin #"<<i<<std::endl;
    //if(drawJohnDSCB)
    sigma_datapoints[i]->GetYaxis()->SetTitle("#sigma Extracted From Fit");
    //else sigma_datapoints[i]->GetYaxis()->SetTitle("Gaussian Core Fit #sigma");
    sigma_datapoints[i]->GetYaxis()->SetTitleOffset(0.95);
    sigma_datapoints[i]->GetYaxis()->SetTitleSize(0.05);
    sigma_datapoints[i]->GetYaxis()->SetLabelOffset(0.005);
    sigma_datapoints[i]->GetYaxis()->SetLabelSize(0.045);
    sigma_datapoints[i]->GetYaxis()->CenterTitle(true);

    sigma_datapoints[i]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");    
    //sigma_datapoints[i]->GetXaxis()->SetTitleOffset(1.0);
    sigma_datapoints[i]->GetXaxis()->SetTitleSize(0.045);
    sigma_datapoints[i]->GetXaxis()->SetLabelOffset(0.005);
    sigma_datapoints[i]->GetXaxis()->SetLabelSize(0.045);

    bool minimumset=false;
    for(int j=1; j<=sigma_datapoints[i]->GetNbinsX();j++){
      //std::cout<<"hello pt bin #"<<i<<std::endl;
      float lowedge  =sigma_datapoints[i]->GetBinLowEdge(j);          //std::cout<<"lowedge ="<<lowedge <<std::endl;
      float highedge =sigma_datapoints[i]->GetBinWidth(j)+lowedge;    //std::cout<<"highedge="<<highedge<<std::endl;
      float sigma       =sigma_datapoints[i]->GetBinContent(j);	      //std::cout<<"mu      ="<<mu      <<std::endl;
      float sigmaerr    =sigma_datapoints[i]->GetBinError(j);      	      //std::cout<<"muerr   ="<<muerr   <<std::endl;
      if(!(sigma>0.))continue;
      if(j<8)continue;//cause i expanded the pT range and now this macro plots everything. i don't want that!
      if(!minimumset){
	ptmin[i]=lowedge;
	//lineptmin[i]=lowedge-(highedge-lowedge);//...why do i do this again? guess we will see...
	lineptmin[i]=lowedge;
	std::cout<<"MINIMUM SET: ptmin[i]="<<ptmin[i]<<std::endl;
	minimumset=true;
      }
      if(ptmax[i]<highedge){
	ptmax[i]=highedge;
	lineptmax[i]=highedge+sigma_datapoints[i]->GetBinWidth(j+1);//...why do i do this again? guess we will see...	  
      }
    }
    
    if(drawPatDSCB)
      for(int j=1; j<=sigma_datapoints_patdscb[i]->GetNbinsX();j++)
	if((sigma_datapoints_patdscb[i]->GetBinCenter(j)>ptmax[i]) ||
	   (sigma_datapoints_patdscb[i]->GetBinCenter(j)<ptmin[i]) 	 )//FOR COSMETICS (points drawing off plot... why?!)
	  sigma_datapoints_patdscb[i]->SetBinContent(j,1.e+20);
    
  }
  
  std::string canvname="PY8JERsigma_onePadOneEta";
  if(!drawFit)canvname+="_datapointsOnly";
  TCanvas* canv=makeSMPRatioCanvas(canvname);
  
  TPaveText* SMPtitle=makeSimPaveTextTitleRatio();
  bool dologx=true;

  for(int i=0; i<netabins; i++){
    
    canv->cd(i+1)->SetLogx(dologx);//maybe set to 1?
    canv->cd(i+1)->SetLogy(0);
    canv->cd(i+1)->SetTicky(1);
    canv->cd(i+1);
    
    float xlo=ptmin[i];
    float xhi=ptmax[i];    
    float linexlo=lineptmin[i];    
    float linexhi=lineptmax[i];    
    
    sigma_datapoints[i]->SetTitle("");
    if(dologx)sigma_datapoints[i]->GetXaxis()->SetNoExponent(true);
    if(dologx)sigma_datapoints[i]->GetXaxis()->SetMoreLogLabels(true);
    sigma_datapoints[i]->SetAxisRange(linexlo, xhi, "X");
    sigma_datapoints[i]->SetAxisRange(0.04,0.21, "Y");
    sigma_datapoints[i]->Draw("HIST E ][");    
    
    if(drawFit){
      sigma_fit[i]->Draw("L SAME");
      if(drawFitUnc){	
	sigma_fit_sysup[i]->Draw("L SAME");
	sigma_fit_sysdown[i]->Draw("L SAME");
	sigma_fit_10xsysup[i]->Draw("L SAME");
	sigma_fit_10xsysdown[i]->Draw("L SAME");
      }
    }
    
    
    //TLine* line1        =makeTLine( linexlo, 0.05 ,linexhi , 0.05);    line1    ->Draw();       
    //TLine* line2        =makeTLine( linexlo, 0.15 ,linexhi , 0.15);    line2    ->Draw();       

    sigma_datapoints[i]->Draw("HIST E ][ SAME");   
    
    if(drawJohnDSCB){
      sigma_datapoints_johndscb[i]->Draw("HIST E ][ SAME");
      if(drawFit&&false)sigma_fit_dscb[i]->Draw("L SAME");
    }
    
    if(drawPatDSCB){
      sigma_datapoints_patdscb[i]->Draw("HIST E ][ SAME");
      if(drawFit&&false)sigma_fit_patdscb[i]->Draw("L SAME");
    }
    
    if(i==1){
      TLegend* leg=makeLegend(0.60, 0.60, 0.98, 0.90);
      leg->AddEntry((TObject*)0, "PYTHIA8", "");
      leg->AddEntry(sigma_datapoints[i],"Gaussian Core Fits","lpe");
      if(drawFit){
	leg->AddEntry(sigma_fit[i],"Fit of Gaussian Core #sigma","l");
	if(drawFitUnc){
	  leg->AddEntry(sigma_fit_sysup[i],"Fit Unc.","l");
	  leg->AddEntry(sigma_fit_10xsysup[i],"10#times(Fit Unc.)","l");
	}
      }

      if(drawJohnDSCB){
	leg->AddEntry(sigma_datapoints_johndscb[i],"from DSCB+Extra Fits","lpe");
	if(drawFit&&false)leg->AddEntry(sigma_fit_dscb[i],"Fit of DSCB+Extra #sigma","lpe");
      }
      
      if(drawPatDSCB){
	leg->AddEntry(sigma_datapoints_patdscb[i],"from DSCB Fits","lpe");
	if(drawFit&&false)leg->AddEntry(sigma_fit_patdscb[i],"Fit of DSCB #sigma","lpe");
      }
      
      leg->Draw();
    }
    
    TPaveText* desc=NULL;
    if(i==0)desc=makePaveText( 0.15, 0.67, 0.45, 0.88);
    else    desc=makePaveText( 0.15, 0.74, 0.45, 0.88);
    desc->AddText(etabin_strs[i].c_str());
    std::string ptrange=std::to_string( (int)xlo) +
                            " GeV < Jet p_{T} < " +
                        std::to_string( (int)xhi) + " GeV";
    desc->AddText(ptrange.c_str());
    if(i==0)desc->AddText(jettype.c_str());
    desc->Draw();
    
    SMPtitle->Draw();  
  }
  
  //makeSMPInclJetXsec_PY8JERsigma_onePadOneEta
  saveCanv(outdir, canv, fout);
  for(int i=0; i<netabins;i++){
    sigma_datapoints[i]->Delete();
    if(drawJohnDSCB)sigma_datapoints_johndscb[i]->Delete();
    if(drawPatDSCB)sigma_datapoints_patdscb[i]->Delete();
  }

  return;
}


//--------------------------------------------------------------------------------------------------------------------------------
void  makeSMPInclJetXsec_PY8JERkLR_onePadOneEta (std::string outdir, TFile* fout, bool drawJohnDSCB){
  std::cout<<"running makeSMPInclJetXsec_PY8JERkLR_onePadOneEta"<<std::endl;

  const int netabins=JERDIR_MC_Nfiles;  
  if(netabins!=n_etabin_strs){
    printetabinerrormessage();
    return;
  }
  
  TH1D* sigma_datapoints[netabins]={};//i'm grabbing this sigma plot for the drawing style only

  TH1D* kL_datapoints_johndscb[netabins]={};
  TH1D* kL_datapoints_patdscb[netabins]={};

  TH1D* kR_datapoints_johndscb[netabins]={};
  TH1D* kR_datapoints_patdscb[netabins]={};

  float lineptmax[netabins]={0.};
  float ptmax[netabins]={0.};
  float lineptmin[netabins]={0.};
  float ptmin[netabins]={0.};
  
  //first get the plots, clone + divide accordingly. binning should be set for me already, essentially
  for(int i=0; i<netabins; i++){
    
    //for the plot style only (gauss fits from here)
    std::string filepath= JERDIR_MC + JERDIR_MC_file_array[i] + ".root";    
    TFile* file=TFile::Open(( filepath).c_str(), "READ");
    std::string sigma_datapoints_name="hSigma_fit";

    //for the actual points
    std::string filepath2=JERDIR_MC + JERDIR2_MC_file_array[i] + ".root";
    TFile* patfile=patfile=TFile::Open( (filepath2).c_str(),"READ");
    
    sigma_datapoints[i] =  (TH1D*)(
				   (
				    (TH1D*)file->Get(sigma_datapoints_name.c_str()) 
				    )->Clone( 
					     (sigma_datapoints_name+"_ybin"+std::to_string(i)).c_str() 
					      )
				   );        
    sigma_datapoints[i]->SetMarkerSize(0);
    sigma_datapoints[i]->SetMarkerColorAlpha(kBlack,0.);
    sigma_datapoints[i]->SetLineColorAlpha(kBlack,0.);
    

    //std::cout<<"hello bin #"<<i<<std::endl;
    sigma_datapoints[i]->GetYaxis()->SetTitle("k_{L} and k_{R} Extracted From Fit");
    sigma_datapoints[i]->GetYaxis()->SetTitleOffset(0.95);
    sigma_datapoints[i]->GetYaxis()->SetTitleSize(0.05);
    sigma_datapoints[i]->GetYaxis()->SetLabelOffset(0.005);
    sigma_datapoints[i]->GetYaxis()->SetLabelSize(0.045);
    sigma_datapoints[i]->GetYaxis()->CenterTitle(true);

    sigma_datapoints[i]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");    
    //sigma_datapoints[i]->GetXaxis()->SetTitleOffset(1.0);
    sigma_datapoints[i]->GetXaxis()->SetTitleSize(0.045);
    sigma_datapoints[i]->GetXaxis()->SetLabelOffset(0.005);
    sigma_datapoints[i]->GetXaxis()->SetLabelSize(0.045);
    

    if(drawJohnDSCB){
      kL_datapoints_johndscb[i] =  (TH1D*)(
					  (
					   (TH1D*)patfile->Get("hJohnDSCBkL") 
					   )->Clone( 
						    ("hJohnDSCBkL_clone_ybin"+std::to_string(i)).c_str() 
						     )
					  );        
      kL_datapoints_johndscb[i]->SetMarkerStyle(kOpenTriangleDown);
      kL_datapoints_johndscb[i]->SetMarkerSize(1);
      kL_datapoints_johndscb[i]->SetMarkerColor(kBlue);
      kL_datapoints_johndscb[i]->SetLineColor(kBlue);
      for(int p=1;p<=kL_datapoints_johndscb[i]->GetNbinsX();p++){
	kL_datapoints_johndscb[i]->SetBinError(p, 1.e-20);//so the marker actually draws.
	if( !(kL_datapoints_johndscb[i]->GetBinContent(p)>0.) &&
	    !(kL_datapoints_johndscb[i]->GetBinContent(p)<0.) )//safe check for 0 bin content
	  kL_datapoints_johndscb[i]->SetBinContent(p, 1.e+20);

      }


      kR_datapoints_johndscb[i] =  (TH1D*)(
					  (
					   (TH1D*)patfile->Get("hJohnDSCBkR") 
					   )->Clone( 
						    ("hJohnDSCBkR_clone_ybin"+std::to_string(i)).c_str() 
						     )
					  );        
      kR_datapoints_johndscb[i]->SetMarkerStyle(kOpenTriangleUp);
      kR_datapoints_johndscb[i]->SetMarkerSize(1);
      kR_datapoints_johndscb[i]->SetMarkerColor(kBlue);
      kR_datapoints_johndscb[i]->SetLineColor(kBlue);
      for(int p=1;p<=kR_datapoints_johndscb[i]->GetNbinsX();p++){
	kR_datapoints_johndscb[i]->SetBinError(p, 1.e-20);//so the marker actually draws.
	if( !(kR_datapoints_johndscb[i]->GetBinContent(p)>0.) &&
	    !(kR_datapoints_johndscb[i]->GetBinContent(p)<0.) )//safe check for 0 bin content
	  kR_datapoints_johndscb[i]->SetBinContent(p, 1.e+20);
      }
      

    }            
    
    kL_datapoints_patdscb[i] =  (TH1D*)(
					     (
					      (TH1D*)patfile->Get("hDSCBkL") 
					      )->Clone( 
						       ("hPatDSCBkL_clone_ybin"+std::to_string(i)).c_str() 
							)
					   );        
    kL_datapoints_patdscb[i]->SetMarkerStyle(kOpenTriangleDown);
    kL_datapoints_patdscb[i]->SetMarkerSize(1);
    kL_datapoints_patdscb[i]->SetMarkerColor(kBlack);
    kL_datapoints_patdscb[i]->SetLineColor(kBlack);
    
    kR_datapoints_patdscb[i] =  (TH1D*)(
					     (
					      (TH1D*)patfile->Get("hDSCBkR") 
					      )->Clone( 
						       ("hPatDSCBkR_clone_ybin"+std::to_string(i)).c_str() 
							)
					   );        
    kR_datapoints_patdscb[i]->SetMarkerStyle(kOpenTriangleUp);
    kR_datapoints_patdscb[i]->SetMarkerSize(1);
    kR_datapoints_patdscb[i]->SetMarkerColor(kBlack);
    kR_datapoints_patdscb[i]->SetLineColor(kBlack);


    //for(int j=1; j<=kL_datapoints_patdscb[i]->GetNbinsX();j++){
    //  std::cout<<"___________________________________________________\n";
    //  std::cout<<"kL_datapoints_patdscb["<<i<<"]->GetBinError("<<j<<")="<<kL_datapoints_patdscb[i]->GetBinError(j)<<"\n";
    //  std::cout<<"kR_datapoints_patdscb["<<i<<"]->GetBinError("<<j<<")="<<kR_datapoints_patdscb[i]->GetBinError(j)<<"\n";
    //}
    
    bool minimumset=false;
    for(int j=1; j<=sigma_datapoints[i]->GetNbinsX();j++){
      //std::cout<<"hello pt bin #"<<i<<std::endl;
      float lowedge  =sigma_datapoints[i]->GetBinLowEdge(j);          //std::cout<<"lowedge ="<<lowedge <<std::endl;
      float highedge =sigma_datapoints[i]->GetBinWidth(j)+lowedge;    //std::cout<<"highedge="<<highedge<<std::endl;
      float sigma       =sigma_datapoints[i]->GetBinContent(j);	      //std::cout<<"mu      ="<<mu      <<std::endl;
      //float sigmaerr    =sigma_datapoints[i]->GetBinError(j);      	      //std::cout<<"muerr   ="<<muerr   <<std::endl;
      if(!(sigma>0.))continue;
      if(j<8)continue;//cause i expanded the pT range and now this macro plots everything. i don't want that!
      if(!minimumset){
	ptmin[i]=lowedge;
	//lineptmin[i]=lowedge-(highedge-lowedge);//...why do i do this again? guess we will see...
	lineptmin[i]=lowedge;
	std::cout<<"MINIMUM SET: ptmin[i]="<<ptmin[i]<<std::endl;
	minimumset=true;
      }
      if(ptmax[i]<highedge){
	ptmax[i]=highedge;
	lineptmax[i]=highedge+sigma_datapoints[i]->GetBinWidth(j+1);//...why do i do this again? guess we will see...	  
      }
    }

    sigma_datapoints[i]->Reset("ICES");

    for(int j=1; j<=kR_datapoints_patdscb[i]->GetNbinsX();j++){
      if((kR_datapoints_patdscb[i]->GetBinCenter(j)>ptmax[i]) ||
	 (kR_datapoints_patdscb[i]->GetBinCenter(j)<ptmin[i]) 	 )//FOR COSMETICS (points drawing off plot... why?!)
	kR_datapoints_patdscb[i]->SetBinContent(j,1.e+20);
      
      if((kL_datapoints_patdscb[i]->GetBinCenter(j)>ptmax[i]) ||
	 (kL_datapoints_patdscb[i]->GetBinCenter(j)<ptmin[i]) 	 )//FOR COSMETICS (points drawing off plot... why?!)
	kL_datapoints_patdscb[i]->SetBinContent(j,1.e+20);
    }

    
    
  }
  
  std::string canvname="PY8JERkLR_onePadOneEta";
  TCanvas* canv=makeSMPRatioCanvas(canvname);
  
  TPaveText* SMPtitle=makeSimPaveTextTitleRatio();
  bool dologx=true;

  for(int i=0; i<netabins; i++){
    
    canv->cd(i+1)->SetLogx(dologx);//maybe set to 1?
    canv->cd(i+1)->SetLogy(0);
    canv->cd(i+1)->SetTicky(1);
    canv->cd(i+1);
    
    float xlo=ptmin[i];
    float xhi=ptmax[i];    
    float linexlo=lineptmin[i];    
    float linexhi=lineptmax[i];    
    
    sigma_datapoints[i]->SetTitle("");
    if(dologx)sigma_datapoints[i]->GetXaxis()->SetNoExponent(true);
    if(dologx)sigma_datapoints[i]->GetXaxis()->SetMoreLogLabels(true);
    sigma_datapoints[i]->SetAxisRange(linexlo, xhi, "X");
    sigma_datapoints[i]->SetAxisRange(0., 2.5, "Y");
    sigma_datapoints[i]->Draw("HIST E ][");    
    
    
    TLine* line1        =makeTLine( linexlo, 1. ,linexhi , 1.);    line1    ->Draw();       
    //TLine* line2        =makeTLine( linexlo, 0.15 ,linexhi , 0.15);    line2    ->Draw();       
    
    //sigma_datapoints[i]->Draw("HIST E ][ SAME");   
    
    if(drawJohnDSCB)
      kL_datapoints_johndscb[i]->Draw("HIST E ][ SAME");        
    kL_datapoints_patdscb[i]->Draw("HIST E ][ SAME");

    if(drawJohnDSCB)
      kR_datapoints_johndscb[i]->Draw("HIST E ][ SAME");        
    kR_datapoints_patdscb[i]->Draw("HIST E ][ SAME");
    
    
    if(i==1){
      TLegend* leg=makeLegend(0.7, 0.70, 0.98, 0.90);
      leg->AddEntry((TObject*)0, "PYTHIA8", "");
      //leg->AddEntry(sigma_datapoints[i],"Gaussian Core Fits","lpe");
      
      if(drawJohnDSCB)
	leg->AddEntry(kR_datapoints_johndscb[i],"k_{R} from DSCB+Extra Fits","lpe");            
      leg->AddEntry(kR_datapoints_patdscb[i],"k_{R} from DSCB Fits","lpe");
      if(drawJohnDSCB)
	leg->AddEntry(kL_datapoints_johndscb[i],"k_{L} from DSCB+Extra Fits","lpe");            
      leg->AddEntry(kL_datapoints_patdscb[i],"k_{L} from DSCB Fits","lpe");
      
      leg->Draw();
    }
    
    TPaveText* desc=NULL;
    if(i==0)desc=makePaveText( 0.15, 0.67, 0.45, 0.88);
    else    desc=makePaveText( 0.15, 0.74, 0.45, 0.88);
    desc->AddText(etabin_strs[i].c_str());
    std::string ptrange=std::to_string( (int)xlo) +
                            " GeV < Jet p_{T} < " +
                        std::to_string( (int)xhi) + " GeV";
    desc->AddText(ptrange.c_str());
    if(i==0)desc->AddText(jettype.c_str());
    desc->Draw();
    
    SMPtitle->Draw();  
  }
  
  //makeSMPInclJetXsec_PY8JERkLR_onePadOneEta
  saveCanv(outdir, canv, fout);
  for(int i=0; i<netabins;i++){
    sigma_datapoints[i]->Delete();
    if(drawJohnDSCB)kL_datapoints_johndscb[i]->Delete();
    kL_datapoints_patdscb[i]->Delete();
    if(drawJohnDSCB)kR_datapoints_johndscb[i]->Delete();
    kR_datapoints_patdscb[i]->Delete();
  }

  return;
}



//--------------------------------------------------------------------------------------------------------------------------------
void  makeSMPInclJetXsec_PY8JERalphaLR_onePadOneEta (std::string outdir, TFile* fout, bool drawJohnDSCB){
  std::cout<<"running makeSMPInclJetXsec_PY8JERalphaLR_onePadOneEta"<<std::endl;

  const int netabins=JERDIR_MC_Nfiles;  
  if(netabins!=n_etabin_strs){
    printetabinerrormessage();
    return;
  }
  
  TH1D* sigma_datapoints[netabins]={};//i'm grabbing this sigma plot for the drawing style only

  TH1D* alphaL_datapoints_johndscb[netabins]={};
  TH1D* alphaL_datapoints_patdscb[netabins]={};

  TH1D* alphaR_datapoints_johndscb[netabins]={};
  TH1D* alphaR_datapoints_patdscb[netabins]={};

  float lineptmax[netabins]={0.};
  float ptmax[netabins]={0.};
  float lineptmin[netabins]={0.};
  float ptmin[netabins]={0.};
  
  //first get the plots, clone + divide accordingly. binning should be set for me already, essentially
  for(int i=0; i<netabins; i++){
    
    //for the plot style only (gauss fits from here)
    std::string filepath= JERDIR_MC + JERDIR_MC_file_array[i] + ".root";    
    TFile* file=TFile::Open(( filepath).c_str(), "READ");
    std::string sigma_datapoints_name="hSigma_fit";

    //for the actual points
    std::string filepath2=JERDIR_MC + JERDIR2_MC_file_array[i] + ".root";
    TFile* patfile=patfile=TFile::Open( (filepath2).c_str(),"READ");
    
    sigma_datapoints[i] =  (TH1D*)(
				   (
				    (TH1D*)file->Get(sigma_datapoints_name.c_str()) 
				    )->Clone( 
					     (sigma_datapoints_name+"_ybin"+std::to_string(i)).c_str() 
					      )
				   );        
    sigma_datapoints[i]->SetMarkerSize(0);
    sigma_datapoints[i]->SetMarkerColorAlpha(kBlack,0.);
    sigma_datapoints[i]->SetLineColorAlpha(kBlack,0.);
    

    //std::cout<<"hello bin #"<<i<<std::endl;
    sigma_datapoints[i]->GetYaxis()->SetTitle("#alpha_{R} and -#alpha_{L} From Fit");
    sigma_datapoints[i]->GetYaxis()->SetTitleOffset(0.95);
    sigma_datapoints[i]->GetYaxis()->SetTitleSize(0.05);
    sigma_datapoints[i]->GetYaxis()->SetLabelOffset(0.005);
    sigma_datapoints[i]->GetYaxis()->SetLabelSize(0.045);
    sigma_datapoints[i]->GetYaxis()->CenterTitle(true);

    sigma_datapoints[i]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");    
    //sigma_datapoints[i]->GetXaxis()->SetTitleOffset(1.0);
    sigma_datapoints[i]->GetXaxis()->SetTitleSize(0.045);
    sigma_datapoints[i]->GetXaxis()->SetLabelOffset(0.005);
    sigma_datapoints[i]->GetXaxis()->SetLabelSize(0.045);
    

    if(drawJohnDSCB){
      alphaL_datapoints_johndscb[i] =  (TH1D*)(
					  (
					   (TH1D*)patfile->Get("hJohnDSCBkL_distFromMu_sigma") 
					   )->Clone( 
						    ("hJohnDSCBkL_distFromMu_sigma_clone_ybin"+std::to_string(i)).c_str() 
						     )
					  );        
      alphaL_datapoints_johndscb[i]->SetMarkerStyle(kOpenTriangleDown);
      alphaL_datapoints_johndscb[i]->SetMarkerSize(1);
      alphaL_datapoints_johndscb[i]->SetMarkerColor(kBlue);
      alphaL_datapoints_johndscb[i]->SetLineColor(kBlue);
      alphaL_datapoints_johndscb[i]->Scale(-1.);
      for(int p=1;p<=alphaL_datapoints_johndscb[i]->GetNbinsX();p++){
	alphaL_datapoints_johndscb[i]->SetBinError(p, 1.e-20);//so the marker actually draws.
	if( !(alphaL_datapoints_johndscb[i]->GetBinContent(p)>0.) &&
	    !(alphaL_datapoints_johndscb[i]->GetBinContent(p)<0.) )//safe check for 0 bin content
	  alphaL_datapoints_johndscb[i]->SetBinContent(p, 1.e+20);	    
      }


      alphaR_datapoints_johndscb[i] =  (TH1D*)(
					  (
					   (TH1D*)patfile->Get("hJohnDSCBkR_distFromMu_sigma") 
					   )->Clone( 
						    ("hJohnDSCBkR_distFromMu_sigma_clone_ybin"+std::to_string(i)).c_str() 
						     )
					  );        
      alphaR_datapoints_johndscb[i]->SetMarkerStyle(kOpenTriangleUp);
      alphaR_datapoints_johndscb[i]->SetMarkerSize(1);
      alphaR_datapoints_johndscb[i]->SetMarkerColor(kBlue);
      alphaR_datapoints_johndscb[i]->SetLineColor(kBlue);
      for(int p=1;p<=alphaR_datapoints_johndscb[i]->GetNbinsX();p++){
	alphaR_datapoints_johndscb[i]->SetBinError(p, 1.e-20);//so the marker actually draws.
	if( !(alphaR_datapoints_johndscb[i]->GetBinContent(p)>0.) &&
	    !(alphaR_datapoints_johndscb[i]->GetBinContent(p)<0.) )//safe check for 0 bin content
	  alphaR_datapoints_johndscb[i]->SetBinContent(p, 1.e+20);
      }

    }            
    
    alphaL_datapoints_patdscb[i] =  (TH1D*)(
					     (
					      (TH1D*)patfile->Get("hDSCBkL_distFromMu_sigma") 
					      )->Clone( 
						       ("hPatDSCBkL_distFromMu_sigma_clone_ybin"+std::to_string(i)).c_str() 
							)
					   );        
    alphaL_datapoints_patdscb[i]->SetMarkerStyle(kOpenTriangleDown);
    alphaL_datapoints_patdscb[i]->SetMarkerSize(1);
    alphaL_datapoints_patdscb[i]->SetMarkerColor(kBlack);
    alphaL_datapoints_patdscb[i]->SetLineColor(kBlack);
    alphaL_datapoints_patdscb[i]->Scale(-1.);
    
    alphaR_datapoints_patdscb[i] =  (TH1D*)(
					     (
					      (TH1D*)patfile->Get("hDSCBkR_distFromMu_sigma") 
					      )->Clone( 
						       ("hPatDSCBkR_distFromMu_sigma_clone_ybin"+std::to_string(i)).c_str() 
							)
					   );        
    alphaR_datapoints_patdscb[i]->SetMarkerStyle(kOpenTriangleUp);
    alphaR_datapoints_patdscb[i]->SetMarkerSize(1);
    alphaR_datapoints_patdscb[i]->SetMarkerColor(kBlack);
    alphaR_datapoints_patdscb[i]->SetLineColor(kBlack);


    //for(int j=1; j<=alphaL_datapoints_patdscb[i]->GetNbinsX();j++){
    //  std::cout<<"___________________________________________________\n";
    //  std::cout<<"alphaL_datapoints_patdscb["<<i<<"]->GetBinError("<<j<<")="<<alphaL_datapoints_patdscb[i]->GetBinError(j)<<"\n";
    //  std::cout<<"alphaR_datapoints_patdscb["<<i<<"]->GetBinError("<<j<<")="<<alphaR_datapoints_patdscb[i]->GetBinError(j)<<"\n";
    //}
    
    bool minimumset=false;
    for(int j=1; j<=sigma_datapoints[i]->GetNbinsX();j++){
      //std::cout<<"hello pt bin #"<<i<<std::endl;
      float lowedge  =sigma_datapoints[i]->GetBinLowEdge(j);          //std::cout<<"lowedge ="<<lowedge <<std::endl;
      float highedge =sigma_datapoints[i]->GetBinWidth(j)+lowedge;    //std::cout<<"highedge="<<highedge<<std::endl;
      float sigma       =sigma_datapoints[i]->GetBinContent(j);	      //std::cout<<"mu      ="<<mu      <<std::endl;
      //float sigmaerr    =sigma_datapoints[i]->GetBinError(j);      	      //std::cout<<"muerr   ="<<muerr   <<std::endl;
      if(!(sigma>0.))continue;
      if(j<8)continue;//cause i expanded the pT range and now this macro plots everything. i don't want that!
      if(!minimumset){
	ptmin[i]=lowedge;
	//lineptmin[i]=lowedge-(highedge-lowedge);//...why do i do this again? guess we will see...
	lineptmin[i]=lowedge;
	std::cout<<"MINIMUM SET: ptmin[i]="<<ptmin[i]<<std::endl;
	minimumset=true;
      }
      if(ptmax[i]<highedge){
	ptmax[i]=highedge;
	lineptmax[i]=highedge+sigma_datapoints[i]->GetBinWidth(j+1);//...why do i do this again? guess we will see...	  
      }
    }        

    sigma_datapoints[i]->Reset("ICES");
    
    for(int j=1; j<=alphaR_datapoints_patdscb[i]->GetNbinsX();j++){
      if((alphaR_datapoints_patdscb[i]->GetBinCenter(j)>ptmax[i]) ||
	 (alphaR_datapoints_patdscb[i]->GetBinCenter(j)<ptmin[i]) 	 )//FOR COSMETICS (points drawing off plot... why?!)
	alphaR_datapoints_patdscb[i]->SetBinContent(j,1.e+20);
      
      if((alphaL_datapoints_patdscb[i]->GetBinCenter(j)>ptmax[i]) ||
	 (alphaL_datapoints_patdscb[i]->GetBinCenter(j)<ptmin[i]) 	 )//FOR COSMETICS (points drawing off plot... why?!)
	alphaL_datapoints_patdscb[i]->SetBinContent(j,1.e+20);
    }
    
    
  }
  
  std::string canvname="PY8JERalphaLR_onePadOneEta";
  TCanvas* canv=makeSMPRatioCanvas(canvname);
  
  TPaveText* SMPtitle=makeSimPaveTextTitleRatio();
  bool dologx=true;

  for(int i=0; i<netabins; i++){
    
    canv->cd(i+1)->SetLogx(dologx);//maybe set to 1?
    canv->cd(i+1)->SetLogy(0);
    canv->cd(i+1)->SetTicky(1);
    canv->cd(i+1);
    
    float xlo=ptmin[i];
    float xhi=ptmax[i];    
    float linexlo=lineptmin[i];    
    float linexhi=lineptmax[i];    
    
    sigma_datapoints[i]->SetTitle("");
    if(dologx)sigma_datapoints[i]->GetXaxis()->SetNoExponent(true);
    if(dologx)sigma_datapoints[i]->GetXaxis()->SetMoreLogLabels(true);
    sigma_datapoints[i]->SetAxisRange(linexlo, xhi, "X");
    sigma_datapoints[i]->SetAxisRange(-6., 9., "Y");
    sigma_datapoints[i]->Draw("HIST E ][");    
    
    TLine* lineneg3        =makeTLine( linexlo, -3. ,linexhi , -3.);    lineneg3    ->Draw();       
    TLine* line0        =makeTLine( linexlo, 0. ,linexhi , 0.);    line0    ->Draw();       
    //TLine* line1        =makeTLine( linexlo, 1. ,linexhi , 1.);    line1    ->Draw();       
    //TLine* line2        =makeTLine( linexlo, 2. ,linexhi , 2.);    line2    ->Draw();       
    TLine* line3        =makeTLine( linexlo, 3. ,linexhi , 3.);    line3    ->Draw();       
    //TLine* line4        =makeTLine( linexlo, 4. ,linexhi , 4.);    line4    ->Draw();       

    
    //sigma_datapoints[i]->Draw("HIST E ][ SAME");   
    
    if(drawJohnDSCB)
      alphaL_datapoints_johndscb[i]->Draw("HIST E ][ SAME");        
    alphaL_datapoints_patdscb[i]->Draw("HIST E ][ SAME");

    if(drawJohnDSCB)
      alphaR_datapoints_johndscb[i]->Draw("HIST E ][ SAME");        
    alphaR_datapoints_patdscb[i]->Draw("HIST E ][ SAME");
    
    
    if(i==1){
      TLegend* leg=makeLegend(0.7, 0.70, 0.98, 0.90);
      leg->AddEntry((TObject*)0, "PYTHIA8", "");
      //leg->AddEntry(sigma_datapoints[i],"Gaussian Core Fits","lpe");
      
      if(drawJohnDSCB)
	leg->AddEntry(alphaR_datapoints_johndscb[i],"#alpha_{R} from DSCB+Extra Fits","lpe");            
      leg->AddEntry(alphaR_datapoints_patdscb[i],"#alpha_{R} from DSCB Fits","lpe");
      if(drawJohnDSCB)
	leg->AddEntry(alphaL_datapoints_johndscb[i],"-#alpha_{L} from DSCB+Extra Fits","lpe");            
      leg->AddEntry(alphaL_datapoints_patdscb[i],"-#alpha_{L} from DSCB Fits","lpe");
      
      leg->Draw();
    }
    
    TPaveText* desc=NULL;
    if(i==0)desc=makePaveText( 0.15, 0.67, 0.45, 0.88);
    else    desc=makePaveText( 0.15, 0.74, 0.45, 0.88);
    desc->AddText(etabin_strs[i].c_str());
    std::string ptrange=std::to_string( (int)xlo) +
                            " GeV < Jet p_{T} < " +
                        std::to_string( (int)xhi) + " GeV";
    desc->AddText(ptrange.c_str());
    if(i==0)desc->AddText(jettype.c_str());
    desc->Draw();
    
    SMPtitle->Draw();  
  }
  
  //makeSMPInclJetXsec_PY8JERalphaLR_onePadOneEta
  saveCanv(outdir, canv, fout);
  for(int i=0; i<netabins;i++){
    sigma_datapoints[i]->Delete();
    if(drawJohnDSCB)alphaL_datapoints_johndscb[i]->Delete();
    alphaL_datapoints_patdscb[i]->Delete();
    if(drawJohnDSCB)alphaR_datapoints_johndscb[i]->Delete();
    alphaR_datapoints_patdscb[i]->Delete();
  }

  return;
}



//--------------------------------------------------------------------------------------------------------------------------------
void  makeSMPInclJetXsec_NLOSmearingFits_onePadAllEta (std::string outdir, TFile* fout){
  std::cout<<"running makeSMPInclJetXsec_NLOSmearingFits_onePadAllEta"<<std::endl;
  
  const int netabins=JERDIR_MC_Nfiles;  
  if(netabins!=n_etabin_strs){
    printetabinerrormessage();
    return;
  }

  TF1* JERfit[netabins]={};
  TH1D* JERfithist[netabins]={};//to make it easier to draw stuff
  float maxy=-1., miny=100000000.;//global min/maxy
  
  //first get the plots, scale accordingly, get the min/max y's
  for(int i=0; i<netabins; i++){
    std::string filepath= JERDIR_MC + JERDIR_MC_file_array[i] + ".root";
    TFile* file=TFile::Open(( filepath).c_str(), "READ");
    
    JERfit[i] =  (TF1*)(
			 (
			  (TF1*)file->Get("SigmaFit_f") 
			  )->Clone( 
				   ("SigmaFit_f_ybin"+std::to_string(i)).c_str() 
				    )
			);
    
    JERfithist[i]=(TH1D*)(
			  JERfit[i]->GetHistogram() 
			  );
    JERfithist[i]->SetName(("SigmaFit_f_hist_ybin"+std::to_string(i)).c_str());
    if(maxy<JERfithist[i]->GetMaximum())      maxy=JERfithist[i]->GetMaximum();
    if(miny>JERfithist[i]->GetMinimum())      miny=JERfithist[i]->GetMinimum();
    
    
  }
  
  // now style hists stuff
  JERfit[0]->SetLineColor(kRed);      
  JERfit[1]->SetLineColor(kGreen);    
  JERfit[2]->SetLineColor(kBlue);     
  JERfit[3]->SetLineColor(kMagenta);  
  
  JERfithist[0]=(TH1D*)JERfithist[0]->Clone("forDrawingOnly_JERfithist");
  JERfithist[0]->Reset("MICES");
  JERfithist[0]->SetLineColorAlpha(kBlack, 1.);
  JERfithist[0]->SetMarkerColorAlpha(kBlack, 1.);
  JERfithist[0]->SetMaximum(maxy*1.2);
  JERfithist[0]->SetMinimum(miny/1.4);
  JERfithist[0]->SetTitle("");  
  JERfithist[0]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  JERfithist[0]->GetXaxis()->SetNoExponent(true);
  JERfithist[0]->GetXaxis()->SetMoreLogLabels(true);  
  JERfithist[0]->GetYaxis()->CenterTitle(true);
  //JERfithist[0]->GetYaxis()->SetTitle("#sigma / #mu");
  //JERfithist[0]->GetYaxis()->SetTitle("#sigma");
  JERfithist[0]->GetYaxis()->SetTitle("Gaussian Core Fit #sigma");

    

  //first hist to be drawn, so this gets the max/min/labels/titles set up
  TLegend* leg=makeLegend(0.55, 0.60, 0.90, 0.85);
  //leg->SetHeader( "PYTHIA8 Gaussian Core Resolution            ","C" );
  leg->SetHeader( "PYTHIA8 Gaussian Core Resolution" );
  leg->AddEntry((TObject*)0,  jettype.c_str(), "");
  
  TCanvas* canv=makeSMPSpectraCanvas("NLOSmearingFits_onePadAllEta");
  canv->SetLogy(0);
  canv->SetLogx(0);
  canv->SetTicky(1);
  canv->cd();
  
  JERfithist[0]->Draw("");
  
  for(int i=0; i<netabins; i++){
    JERfit[i]->Draw("SAME");
    std::string legstr=etabin_strs[i] ;
    leg->AddEntry( JERfit[i], 
		   (legstr).c_str() ,"lp");        
    
    
  }
  
  leg->Draw();
  
  TPaveText* SMPtitle=makeSimPaveTextTitle();
  SMPtitle->Draw();  
  
  //makeSMPInclJetXsec_NLOSmearingFits_onePadAllEta
  saveCanv(outdir, canv, fout);

  for(int i=0; i<netabins;i++){
    JERfit[i]->Delete();
  }
  
  return;
}



//--------------------------------------------------------------------------------------------------------------------------------
void  makeSMPInclJetXsec_NLOSmearingFitswSFs_onePadAllEta (std::string outdir, TFile* fout){
  std::cout<<"running makeSMPInclJetXsec_NLOSmearingFitswSFs_onePadAllEta"<<std::endl;
  
  const int netabins=JERDIR_MC_Nfiles;  
  if(netabins!=n_etabin_strs){
    printetabinerrormessage();
    return;
  }

  TF1* JERfit[netabins]={};
  TH1D* JERfithist[netabins]={};//to make it easier to draw stuff
  float maxy=-1., miny=100000000.;//global min/maxy
  
  //first get the plots, scale accordingly, get the min/max y's
  for(int i=0; i<netabins; i++){
    std::string filepath= JERDIR_MC + JERDIR_MC_file_array[i] + ".root";
    TFile* file=TFile::Open(( filepath).c_str(), "READ");
    
    JERfit[i] =  (TF1*)(
			 (
			  (TF1*)file->Get("SigmaFit_f") 
			  )->Clone( 
				   ("SigmaFit_f_ybin"+std::to_string(i)).c_str() 
				    )
			);
    JERfit[i]->SetParameter(0,
			    JERfit[i]->GetParameter(0)*1.1);
    JERfit[i]->SetParameter(1,
			    JERfit[i]->GetParameter(1)*1.1);
    
    JERfithist[i]=(TH1D*)(
			  JERfit[i]->GetHistogram() 
			  );
    JERfithist[i]->SetName(("SigmaFit_f_hist_ybin"+std::to_string(i)).c_str());
    if(maxy<JERfithist[i]->GetMaximum())      maxy=JERfithist[i]->GetMaximum();
    if(miny>JERfithist[i]->GetMinimum())      miny=JERfithist[i]->GetMinimum();
    
    
  }
  
  // now style hists stuff
  JERfit[0]->SetLineColor(kRed);      
  JERfit[1]->SetLineColor(kGreen);    
  JERfit[2]->SetLineColor(kBlue);     
  JERfit[3]->SetLineColor(kMagenta);  
  
  JERfithist[0]=(TH1D*)JERfithist[0]->Clone("forDrawingOnly_JERfithistwSFs");
  JERfithist[0]->Reset("MICES");
  JERfithist[0]->SetLineColorAlpha(kBlack, 1.);
  JERfithist[0]->SetMarkerColorAlpha(kBlack, 1.);
  JERfithist[0]->SetMaximum(maxy*1.2);
  JERfithist[0]->SetMinimum(miny/1.4);
  JERfithist[0]->SetTitle("");  
  JERfithist[0]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  JERfithist[0]->GetXaxis()->SetNoExponent(true);
  JERfithist[0]->GetXaxis()->SetMoreLogLabels(true);  
  //JERfithist[0]->GetYaxis()->SetTitle("#sigma / #mu");
  JERfithist[0]->GetYaxis()->SetTitle("#sigma");
    

  //first hist to be drawn, so this gets the max/min/labels/titles set up
  TLegend* leg=makeLegend(0.55, 0.60, 0.90, 0.85);
  //leg->SetHeader( "PYTHIA8 Gaussian Core Resolution            ","C" );
  leg->SetHeader( "PYTHIA8 JER Smearing Fits" );
  leg->AddEntry((TObject*)0,  "Scale Factor = 1.1 +/- 0.1", "");
  leg->AddEntry((TObject*)0,  jettype.c_str(), "");
  
  TCanvas* canv=makeSMPSpectraCanvas("NLOSmearingFitswSFs_onePadAllEta");
  canv->SetLogy(0);
  canv->SetLogx(0);
  canv->SetTicky(1);
  canv->cd();
  
  JERfithist[0]->Draw("");
  
  for(int i=0; i<netabins; i++){
    JERfit[i]->Draw("SAME");
    std::string legstr=etabin_strs[i] ;
    leg->AddEntry( JERfit[i], 
		   (legstr).c_str() ,"lp");        
    
    
  }
  
  leg->Draw();
  
  TPaveText* SMPtitle=makeSimPaveTextTitle();
  SMPtitle->Draw();  
  
  //makeSMPInclJetXsec_NLOSmearingFitswSFs_onePadAllEta
  saveCanv(outdir, canv, fout);

  for(int i=0; i<netabins;i++){
    JERfit[i]->Delete();
  }    
  
  return;
}

//--------------------------------------------------------------------------------------------------------------------------------
void  makeSMPInclJetXsec_NLOSmearingFits_onePadOneEta (std::string outdir, TFile* fout){
  std::cout<<"running makeSMPInclJetXsec_NLOSmearingFits_onePadOneEta"<<std::endl;

  const int netabins=JERDIR_MC_Nfiles;  
  if(netabins!=n_etabin_strs){
    printetabinerrormessage();
    return;
  }
  
  TF1* JERfit[netabins]={};	       
  TF1* JERfit_sysup[netabins]={};      
  TF1* JERfit_sysdown[netabins]={};    
                                       
  TH1D* sigmu_datapoints[netabins]={}; 
  TH1D* sigma_datapoints[netabins]={}; 
  TH1D* mu_datapoints[netabins]={};    
  float lineptmax[netabins]={0.};
  float ptmax[netabins]={0.};
  float lineptmin[netabins]={0.};
  float ptmin[netabins]={0.};             
  
  //first get the plots, clone + divide accordingly. binning should be set for me already, essentially
  for(int i=0; i<netabins; i++){
    std::string filepath= JERDIR_MC + JERDIR_MC_file_array[i] + ".root";
    TFile* file=TFile::Open(( filepath).c_str(), "READ");
    
    JERfit[i] =  (TF1*)(
			(
			 (TF1*)file->Get("SigmaFit_f") 
			 )->Clone( 
				  ("SigmaFit_f_ybin"+std::to_string(i)).c_str() 
				   )
			);
    
    JERfit_sysup[i] =  (TF1*)(
			      (
			       (TF1*)file->Get("SigmaFit_sysup_f") 
			       )->Clone( 
					("SigmaFit_sysup_f_ybin"+std::to_string(i)).c_str() 
					 )

			      );
    
    JERfit_sysdown[i] =  (TF1*)(
				(
			 (TF1*)file->Get("SigmaFit_sysdown_f") 
			 )->Clone( 
				  ("SigmaFit_sysdown_f_ybin"+std::to_string(i)).c_str() 
				   )
				
				);
    
    sigma_datapoints[i] =  (TH1D*)(
				   (
				    (TH1D*)file->Get("hSigma_fit") 
				    )->Clone( 
					     ("hSigma_datapoints_ybin"+std::to_string(i)).c_str() 
					      )
				   );        

    mu_datapoints[i] =  (TH1D*)(
				(
				 (TH1D*)file->Get("hMean_fit") 
				    )->Clone( 
					     ("hMean_datapoints_ybin"+std::to_string(i)).c_str() 
					      )
				);        
    
    sigmu_datapoints[i]=(TH1D*) sigma_datapoints[i]->Clone(
							   ("hSigmu_datapoints_ybin"+std::to_string(i)).c_str()
							   );
    //std::cout<<"hello bin #"<<i<<std::endl;
    sigmu_datapoints[i]->Reset("MICES");
    sigmu_datapoints[i]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");

    //sigmu_datapoints[i]->GetYaxis()->SetTitle("#sigma / #mu");
    //sigmu_datapoints[i]->GetYaxis()->SetTitle("#sigma");
    sigmu_datapoints[i]->GetYaxis()->SetTitle("Gaussian Core Fit #sigma");
    sigmu_datapoints[i]->GetYaxis()->SetTitleOffset(0.95);
    sigmu_datapoints[i]->GetYaxis()->SetTitleSize(0.05);
    sigmu_datapoints[i]->GetYaxis()->SetLabelOffset(0.005);
    sigmu_datapoints[i]->GetYaxis()->SetLabelSize(0.045);
    sigmu_datapoints[i]->GetYaxis()->CenterTitle(true);

    bool minimumset=false;
    for(int j=1; j<=sigma_datapoints[i]->GetNbinsX();j++){
      //std::cout<<"hello pt bin #"<<i<<std::endl;
      float lowedge  =sigma_datapoints[i]->GetBinLowEdge(j);          //std::cout<<"lowedge ="<<lowedge <<std::endl;
      float highedge =sigma_datapoints[i]->GetBinWidth(j)+lowedge;    //std::cout<<"highedge="<<highedge<<std::endl;
      float sigma    =sigma_datapoints[i]->GetBinContent(j);	      //std::cout<<"sigma   ="<<sigma   <<std::endl;
      float sigmaerr =sigma_datapoints[i]->GetBinError(j);	      //std::cout<<"sigmaerr="<<sigmaerr<<std::endl;
      float mu       =mu_datapoints[i]->GetBinContent(j);	      //std::cout<<"mu      ="<<mu      <<std::endl;
      float muerr    =mu_datapoints[i]->GetBinError(j);      	      //std::cout<<"muerr   ="<<muerr   <<std::endl;      

      sigmu_datapoints[i]->SetBinContent(j,sigma);//thisis happening because i reset the hist for drawing purposes i think. 
      sigmu_datapoints[i]->SetBinError(j,sigmaerr);            

      if(!(mu>0.))continue;
      if(j<7)continue;//cause i expanded the pT range and now this macro plots everything. i don't want that!
      //if(j<8)continue;//cause i expanded the pT range and now this macro plots everything. i don't want that!
      float sigmu=sigma/mu;
      float sigmuerr=(sigma/mu)*sqrt(
				     pow(sigmaerr/sigma, 2)+
				     pow(muerr/mu, 2)
				     );
      
      //sigmu_datapoints[i]->SetBinContent(j,sigmu);
      //sigmu_datapoints[i]->SetBinError(j,sigmuerr);            
      if(!minimumset){
	ptmin[i]=lowedge;
	lineptmin[i]=lowedge-(highedge-lowedge);
	std::cout<<"MINIMUM SET: ptmin[i]="<<ptmin[i]<<std::endl;
	minimumset=true;
      }
      if(ptmax[i]<highedge){
	ptmax[i]=highedge;
	lineptmax[i]=highedge+sigma_datapoints[i]->GetBinWidth(j+1);
	  
      }
    }
    
    

  }
  
  
  
  TCanvas* canv=makeSMPRatioCanvas("NLOSmearingFits_onePadOneEta");

  TPaveText* SMPtitle=makeSimPaveTextTitleRatio();
  bool dologx=true;
  for(int i=0; i<netabins; i++){

    canv->cd(i+1)->SetLogx(dologx);//maybe set to 1?
    canv->cd(i+1)->SetLogy(0);
    canv->cd(i+1)->SetTicky(1);
    canv->cd(i+1);

    float xlo=ptmin[i];
    float xhi=ptmax[i];    
    float linexlo=lineptmin[i];    
    float linexhi=lineptmax[i];    
    
    sigmu_datapoints[i]->SetTitle("");
    if(dologx)sigmu_datapoints[i]->GetXaxis()->SetNoExponent(true);
    if(dologx)sigmu_datapoints[i]->GetXaxis()->SetMoreLogLabels(true);
    sigmu_datapoints[i]->SetAxisRange(linexlo, xhi, "X");
    sigmu_datapoints[i]->SetAxisRange(0.04,0.23, "Y");
    sigmu_datapoints[i]->Draw("HIST E ][");    

    
    //TLine* p05     =makeTLine( linexlo,.05 ,linexhi , .05);    p05     ->Draw();       
    //TLine* p10     =makeTLine( linexlo,.10 ,linexhi , .10);    p10     ->Draw();       
    //TLine* p15     =makeTLine( linexlo,.15 ,linexhi , .15);    p15     ->Draw();       
    sigmu_datapoints[i]->Draw("HIST E ][ SAME");    
    JERfit[i]        ->Draw("SAME");    
    JERfit_sysup[i]  ->Draw("SAME");    
    JERfit_sysdown[i]->Draw("SAME");    
    

    
    
    TPaveText* desc=NULL;
    if(i==0)desc=makePaveText( 0.15, 0.67, 0.45, 0.88);
    else    desc=makePaveText( 0.15, 0.74, 0.45, 0.88);
    desc->AddText(etabin_strs[i].c_str());
    std::string ptrange=std::to_string( (int)xlo)+
      " GeV < Jet p_{T} < "+
      std::to_string( (int)xhi)+" GeV";
    desc->AddText(ptrange.c_str());
    if(i==0)desc->AddText(jettype.c_str());
    desc->Draw();

    
    if(i==1){
      TLegend* leg=makeLegend(0.56, 0.68, 0.87, 0.87);
      leg->AddEntry(sigmu_datapoints[i],"PYTHIA8 Resolution","lpe");
      leg->AddEntry(JERfit[i],"Gauss Core Fit Resolution","lp");
      leg->AddEntry(JERfit_sysdown[i],"Upper/Lower Fit Unc.","l");
      leg->Draw();
    }
    
    
    
    SMPtitle->Draw();  
  }
  
  //makeSMPInclJetXsec_NLOSmearingFits_onePadOneEta
  saveCanv(outdir, canv, fout);

  for(int i=0; i<netabins;i++){
    JERfit[i]->Delete();	           
    JERfit_sysup[i]->Delete();      
    JERfit_sysdown[i]->Delete();    
    
    sigmu_datapoints[i]->Delete(); 
    sigma_datapoints[i]->Delete(); 
    mu_datapoints[i]->Delete();    
  }


  return;
}


//--------------------------------------------------------------------------------------------------------------------------------
void  makeSMPInclJetXsec_NLOSmearingFitswSFs_onePadOneEta (std::string outdir, TFile* fout, bool drawJohnDSCB){
  std::cout<<"running makeSMPInclJetXsec_NLOSmearingFitswSFs_onePadOneEta"<<std::endl;

  const int netabins=JERDIR_MC_Nfiles;  
  if(netabins!=n_etabin_strs){
    printetabinerrormessage();
    return;
  }
  
  TF1* JERfit[netabins]={};
  TF1* JERfit_sysup[netabins]={};
  TF1* JERfit_sysdown[netabins]={};
  
  TH1D* sigmu_datapoints[netabins]={};
  TH1D* sigma_datapoints[netabins]={};
  TH1D* mu_datapoints[netabins]={};
  float lineptmax[netabins]={0.};
  float ptmax[netabins]={0.};
  float lineptmin[netabins]={0.};
  float ptmin[netabins]={0.};

  std::string fitname="SigmaFit_f";
  if(drawJohnDSCB)fitname="SigmaFit_dscb";
  
  std::string sigma_datapoints_name="hSigma_fit";
  if(drawJohnDSCB)sigma_datapoints_name="hJohnDSCBSigma_clone";
  
  std::string mu_datapoints_name="hMean_fit";
  if(drawJohnDSCB)mu_datapoints_name="hJohnDSCBMean_clone";
  
  //first get the plots, clone + divide accordingly. binning should be set for me already, essentially
  for(int i=0; i<netabins; i++){
    std::string filepath= JERDIR_MC + JERDIR_MC_file_array[i] + ".root";
    TFile* file=TFile::Open(( filepath).c_str(), "READ");
    
    JERfit[i] =  (TF1*)(
			(
			 (TF1*)file->Get(fitname.c_str()) 
			 )->Clone( 
				  (fitname+"_ybin"+std::to_string(i)).c_str() 
				   )
			);    
    
    JERfit_sysup[i] =  (TF1*)( JERfit[i]->Clone( 
						(fitname+"_sysup_ybin"+std::to_string(i)).c_str() 
						 )			       
			       );
    
    JERfit_sysdown[i] =  (TF1*)( JERfit[i]->Clone( 
						  (fitname+"_sysdown_ybin"+std::to_string(i)).c_str() 
						   )				 
				 );
    
    JERfit[i]->SetParameter(0,
			    JERfit[i]->GetParameter(0)*(1.1));
    JERfit[i]->SetParameter(1,
			    JERfit[i]->GetParameter(1)*(1.1));

    JERfit_sysup[i]->SetParameter(0,
				  JERfit_sysup[i]->GetParameter(0)*(1.2));
    JERfit_sysup[i]->SetParameter(1,
				  JERfit_sysup[i]->GetParameter(1)*(1.2));
    JERfit_sysup[i]->SetLineColor(kRed);

    JERfit_sysdown[i]->SetParameter(0,
				  JERfit_sysdown[i]->GetParameter(0)*(1.0));
    JERfit_sysdown[i]->SetParameter(1,
				  JERfit_sysdown[i]->GetParameter(1)*(1.0));
    JERfit_sysdown[i]->SetLineColor(kRed);
    
    
    sigma_datapoints[i] =  (TH1D*)(
				   (
				    (TH1D*)file->Get(sigma_datapoints_name.c_str()) 
				    )->Clone( 
					     (sigma_datapoints_name+"_ybin"+std::to_string(i)).c_str() 
					      )
				   );        

    mu_datapoints[i] =  (TH1D*)(
				(
				 (TH1D*)file->Get(mu_datapoints_name.c_str()) 
				    )->Clone( 
					     (mu_datapoints_name+"_ybin"+std::to_string(i)).c_str() 
					      )
				);        
    
    sigmu_datapoints[i]=(TH1D*) sigma_datapoints[i]->Clone(
							   ("hSigmu_datapoints_ybin"+std::to_string(i)).c_str()
							   );
    //std::cout<<"hello bin #"<<i<<std::endl;
    sigmu_datapoints[i]->Reset("MICES");
    sigmu_datapoints[i]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
    //sigmu_datapoints[i]->GetYaxis()->SetTitle("#sigma / #mu");
    sigmu_datapoints[i]->GetYaxis()->SetTitle("#sigma");
    sigmu_datapoints[i]->SetMarkerColorAlpha(kBlack, 0.);
    sigmu_datapoints[i]->SetLineColorAlpha(kBlack, 0.);
    
    bool minimumset=false;
    for(int j=1; j<=sigma_datapoints[i]->GetNbinsX();j++){
      //std::cout<<"hello pt bin #"<<i<<std::endl;
      float lowedge  =sigma_datapoints[i]->GetBinLowEdge(j);          //std::cout<<"lowedge ="<<lowedge <<std::endl;
      float highedge =sigma_datapoints[i]->GetBinWidth(j)+lowedge;    //std::cout<<"highedge="<<highedge<<std::endl;
      float sigma    =sigma_datapoints[i]->GetBinContent(j);	      //std::cout<<"sigma   ="<<sigma   <<std::endl;
      float sigmaerr =sigma_datapoints[i]->GetBinError(j);	      //std::cout<<"sigmaerr="<<sigmaerr<<std::endl;
      float mu       =mu_datapoints[i]->GetBinContent(j);	      //std::cout<<"mu      ="<<mu      <<std::endl;
      float muerr    =mu_datapoints[i]->GetBinError(j);      	      //std::cout<<"muerr   ="<<muerr   <<std::endl;
      sigmu_datapoints[i]->SetBinContent(j,sigma);
      sigmu_datapoints[i]->SetBinError(j,sigmaerr);            

      if(!(mu>0.))continue;
      //if(j<7)continue;//cause i expanded the pT range and now this macro plots everything. i don't want that!
      if(j<8)continue;//cause i expanded the pT range and now this macro plots everything. i don't want that!
      float sigmu=sigma/mu;
      float sigmuerr=(sigma/mu)*sqrt(
				     pow(sigmaerr/sigma, 2)+
				     pow(muerr/mu, 2)
				     );
      
      //sigmu_datapoints[i]->SetBinContent(j,sigmu);
      //sigmu_datapoints[i]->SetBinError(j,sigmuerr);            
      if(!minimumset){
	ptmin[i]=lowedge;
	lineptmin[i]=lowedge;//-(highedge-lowedge);
	std::cout<<"MINIMUM SET: ptmin[i]="<<ptmin[i]<<std::endl;
	minimumset=true;
      }
      if(ptmax[i]<highedge){
	ptmax[i]=highedge;
	lineptmax[i]=highedge+sigma_datapoints[i]->GetBinWidth(j+1);
	  
      }
    }
    
    

  }
  
  
  std::string canvname="NLOSmearingFitswSFs_onePadOneEta";
  if(drawJohnDSCB)canvname="DSCBNLOSmearingFitswSFs_onePadOneEta";
  //TCanvas* canv=makeSMPRatioCanvas("NLOSmearingFitswSFs_onePadOneEta");
  TCanvas* canv=makeSMPRatioCanvas(canvname);

  TPaveText* SMPtitle=makeSimPaveTextTitleRatio();
  bool dologx=true;
  for(int i=0; i<netabins; i++){

    canv->cd(i+1)->SetLogx(dologx);//maybe set to 1?
    canv->cd(i+1)->SetLogy(0);
    canv->cd(i+1)->SetTicky(1);
    canv->cd(i+1);

    float xlo=ptmin[i];
    float xhi=ptmax[i];    
    float linexlo=lineptmin[i];    
    float linexhi=lineptmax[i];    
    
    sigmu_datapoints[i]->SetTitle("");
    if(dologx)sigmu_datapoints[i]->GetXaxis()->SetNoExponent(true);
    if(dologx)sigmu_datapoints[i]->GetXaxis()->SetMoreLogLabels(true);
    sigmu_datapoints[i]->SetAxisRange(linexlo, xhi, "X");
    sigmu_datapoints[i]->SetAxisRange(0.04,0.23, "Y");
    sigmu_datapoints[i]->Draw("HIST E ][");    

    
    //TLine* p05     =makeTLine( linexlo,.05 ,linexhi , .05);    p05     ->Draw();       
    //TLine* p10     =makeTLine( linexlo,.10 ,linexhi , .10);    p10     ->Draw();       
    //TLine* p15     =makeTLine( linexlo,.15 ,linexhi , .15);    p15     ->Draw();       
    sigmu_datapoints[i]->Draw("HIST E ][ SAME");    
    JERfit[i]        ->Draw("SAME");    
    JERfit_sysup[i]  ->Draw("SAME");    
    JERfit_sysdown[i]->Draw("SAME");    
    

    
    
    
    TPaveText* desc=NULL;
    if(i==0)desc=makePaveText( 0.15, 0.67, 0.45, 0.88);
    else    desc=makePaveText( 0.15, 0.74, 0.45, 0.88);
    desc->AddText(etabin_strs[i].c_str());
    std::string ptrange=std::to_string( (int)xlo)+
      " GeV < Jet p_{T} < "+
      std::to_string( (int)xhi)+" GeV";
    desc->AddText(ptrange.c_str());
    if(i==0)desc->AddText(jettype.c_str());
    desc->Draw();
    
    if(i==1){
      TLegend* leg=makeLegend(0.56, 0.68, 0.87, 0.87);
      if(drawJohnDSCB)leg->AddEntry((TObject*)0,"DSCB Smearing Fits","");
      else        leg->AddEntry((TObject*)0,"Gauss Core Smearing Fits","");
      leg->AddEntry((TObject*)0,"Scale Factor = 1.1 +/- 0.1","");
      leg->AddEntry(JERfit[i],"Fit","lp");
      leg->AddEntry(JERfit_sysdown[i],"Upper/Lower Unc.","l");
      leg->Draw();
      
    }
    
    
    SMPtitle->Draw();  
  }
  
  //makeSMPInclJetXsec_NLOSmearingFitswSFs_onePadOneEta
  saveCanv(outdir, canv, fout);

  for(int i=0; i<netabins;i++){
    JERfit[i]->Delete();	           
    JERfit_sysup[i]->Delete();      
    JERfit_sysdown[i]->Delete();    
    
    sigmu_datapoints[i]->Delete(); 
    sigma_datapoints[i]->Delete(); 
    mu_datapoints[i]->Delete();    
  }

  
  return;
}

