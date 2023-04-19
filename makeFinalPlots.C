//what this code is supposed to do, is take in the final root files, and make additional plots/canvases/etc. 
// .png --> for quick reference + adding things to ANs/Theses/slides/etc...
// .C --> for quick/last minute editing figure etc.

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <string>
#include <cstddef>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TKey.h"
#include "TClass.h"
#include "TColor.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TAttLine.h"
#include "TAxis.h"
#include "TF1.h"
#include "TObject.h"
#include "TMemStatShow.h"

#include "printCanvases_methods.h"
//#include "file_strings_03.18.20.h"//Y RESULTS W NLO FROM JOAO
#include "file_strings_01.05.21.h"//Y RESULTS W FULL-STAT NLO FROM KLAUS
#include "makeSMPInclJetXsec.h"


const std::string uncopts[1]={//"",
			      //"datatotUnc",
			      //"NLOtotUnc",
			      "totUnc"};

const std::string targpdfs[3]={ "CT14nnlo",
				"NNPDF31_nnlo_as_0118",
				"NNPDF31_nnlo_as_0120"};


//const std::string scaleopts[2]={"murmufHTp_v4",
//				"murmufpt_v4"};
const std::string scaleopts[2]={"murmufHTp_v5",
				"murmufpt_v5"};
const std::string orders[2]={"NLO","NNLO"};

const std::string outdir="final_plots/";
const bool writeToFile=false;
const bool printFileBaseToScreen=true;

const bool doDataPlots  =true;//upon review of output 9/13/2021; good
const bool doPY8Plots   =false;//upon review of output 9/13/2021; good
const bool doPY8unfPlots=false;//upon review of output 9/13/2021; good
const bool doThyPlots   =false;//
const bool doNLOunfPlots=false;//

void makeFinalPlots(){
  
  
  TFile* fout=NULL;
  if(writeToFile)
    fout=new TFile("final_plots/final_plots.root","RECREATE");

  
  if(printFileBaseToScreen){
    std::cout<<"NLO_UNF_DATA_file_base="<< NLO_UNF_DATA_file_base<<std::endl;
    std::cout<<"NLO_UNF_CLOSURE_file_base="<< NLO_UNF_DATA_file_base<<std::endl;    

    std::cout<<"PY8_UNF_DATA_file_base="<< PY8_UNF_DATA_file_base<<std::endl;
    std::cout<<"PY8_UNF_CLOSURE_file_base="<< PY8_UNF_DATA_file_base<<std::endl;
  }
  
  //////makeSMPInclJetXsec_data.h //good
  if(doDataPlots){
    //////Data/PY8 event counts
    //printEvtCountTable(outdir);
    ////////unfolded data pT bins
    //printpTbins(outdir);
    ////////Data Stat Errors
    //makeSMPInclJetXsec_fracstaterrData_ratios(outdir,fout);
    ////////Data Covariance Matrix
    //makeSMPInclJetXsec_covmatData_onePadOneEta(outdir,fout);
    ////NLO v PY8 unfolding comparison, MY unfolded results
    //makeSMPInclJetXsec_PY8vNLOunfdata_ratios(outdir,fout);    
    makeSMPInclJetXsec_PY8vNLOunfdata_modelsys_ratios(outdir,fout,false,false);    //fine
    makeSMPInclJetXsec_PY8vNLOunfdata_modelsys_ratios(outdir,fout, true,false);    //fine
    makeSMPInclJetXsec_PY8vNLOunfdata_modelsys_ratios(outdir,fout, true,true);    //fine
  }

  
  //////makeSMPInclJetXsec_PY8.h  //good
  if(doPY8Plots){

    bool drawFit=true;
    bool drawJohnDSCB=true;
    bool drawPatDSCB=true;
    bool drawFitUnc=true;//for mu fit only.


    ////////PY8 JER mu
    //makeSMPInclJetXsec_PY8JERmu_onePadOneEta(   outdir, fout, drawFit, drawFitUnc, false, false);//good
    //makeSMPInclJetXsec_PY8JERmu_onePadOneEta(   outdir, fout, false  , false     , drawJohnDSCB,drawPatDSCB);//good
    //
    ////////PY8 JER sigma
    //makeSMPInclJetXsec_PY8JERsigma_onePadOneEta(outdir, fout, drawFit, drawFitUnc,  false,      false);//good
    //makeSMPInclJetXsec_PY8JERsigma_onePadOneEta(outdir, fout, false,  false, drawJohnDSCB,drawPatDSCB);//good
    //
    ////PY8 JER k_L/R
    //makeSMPInclJetXsec_PY8JERkLR_onePadOneEta(   outdir, fout, drawJohnDSCB);//good
    //    
    ////PY8 JER |k_L/R-mu|/sigma (|alpha_L/R|)
    //makeSMPInclJetXsec_PY8JERalphaLR_onePadOneEta(   outdir, fout, drawJohnDSCB);
    //return;

    //////PY8 JER fits
    makeSMPInclJetXsec_NLOSmearingFits_onePadAllEta(    outdir, fout);//meh why bother anymore
    makeSMPInclJetXsec_NLOSmearingFitswSFs_onePadAllEta(outdir, fout);//meh why bother anymore
    makeSMPInclJetXsec_NLOSmearingFitswSFs_onePadOneEta(outdir,fout, false);
    
    return;    
    
    //CRYSTAL BALL FIT VERSIONS HERE??
  }
  
  
  
  
  //////makeSMPInclJetXsec_PY8unf.h //good
  if(doPY8unfPlots){
    ////comparisons with thy/GEN/truth distributions
    makeSMPInclJetXsec_PY8unfdata(outdir,fout);
    makeSMPInclJetXsec_PY8unfdata_ratios(outdir,fout,false);
    makeSMPInclJetXsec_PY8unfdata_ratios(outdir,fout,true);//use nlo unf data instead
    makeSMPInclJetXsec_PY8unfdatasysterr_ratios(outdir,fout);
    //comparisons with measured data distributions
    makeSMPInclJetXsec_PY8unfdata_wdatameas(outdir,fout);
    makeSMPInclJetXsec_PY8unfdata_wdatameas_ratios(outdir,fout);
    ////////response matrices for this unfolding
    makeSMPInclJetXsec_PY8unfrespmat_onePadOneEta(outdir,fout);
    ////////pearson matrices for this unfolding
    makeSMPInclJetXsec_PY8unfpearsonmat_onePadOneEta(outdir,fout);
    //////////closure+folding ratios
    makeSMPInclJetXsec_PY8unf_closure_ratios(outdir,fout);
    makeSMPInclJetXsec_PY8unf_folding_ratios(outdir,fout);
    //////////fake and misses w/ PY8 GEN and RECO
    makeSMPInclJetXsec_PY8unf_missAndFakes(outdir,fout);
    makeSMPInclJetXsec_PY8unf_chi2viter(outdir,fout);    
    
    //return;
  }
  
  
  //////makeSMPInclJetXsec_NLOandNPCs.h 
  if(doThyPlots){
    //NLO spectra
    //makeSMPInclJetXsec_OnlyNLO(outdir,fout);//good

    //////PY8 NPCs
    //makeSMPInclJetXsec_NPCs_onePadAllEta(outdir,fout, "PYTHIA8");
    //makeSMPInclJetXsec_NPCs_onePadOneEta(outdir,fout, "PYTHIA8");
    ////////HERWIG NPCs
    //makeSMPInclJetXsec_NPCs_onePadAllEta(outdir,fout, "HERWIG");
    //makeSMPInclJetXsec_NPCs_onePadOneEta(outdir,fout, "HERWIG");
    ////////AVG NPCs
    //makeSMPInclJetXsec_NPCs_onePadAllEta(outdir,fout, "AVG");
    //makeSMPInclJetXsec_NPCs_onePadOneEta(outdir,fout, "AVG");


    //return;
    
    ////NLOsyst targPDF v order
    for(int l=0; l<2;l++){	//loop over orders
      for(int k=0; k<2;k++){//loop over scale opts
	for(int j=0; j<3;j++){	    //loop over targ pdfs
	  makeSMPInclJetXsec_NLOsyst_targPDF_ratios(outdir,fout, targpdfs[j], scaleopts[k], orders[l]);//add option to use PY8 unfolded data instead
	}}}

    return;
  }
  
  
  //////////makeSMPInclJetXsec_NLOunf.h 
  if(doNLOunfPlots){

    ////////spectra comparisons with thy/GEN/truth distributions
    //makeSMPInclJetXsec_NLOunfdata(outdir,fout);//fine       
    ////////syst err sources on unf data only 
    //makeSMPInclJetXsec_NLOunfdatasysterr_ratios(outdir,fout); //fine
    ////////comparisons with measured data distributions 
    //makeSMPInclJetXsec_NLOunfdata_wdatameas(outdir,fout);//fine
    //makeSMPInclJetXsec_NLOunfdata_wdatameas_ratios(outdir,fout);//fine
    //////////Response Matrices  //fine
    //makeSMPInclJetXsec_NLOunfrespmat_onePadOneEta(outdir,fout);//fine
    //////////closure+folding ratios
    //makeSMPInclJetXsec_NLOunf_closure_ratios(outdir,fout);//fine
    //makeSMPInclJetXsec_NLOunf_folding_ratios(outdir,fout);//fine
    ////fake and misses w/ truth+smeared spectra of toy NLO
    //makeSMPInclJetXsec_NLOunf_missAndFakes(outdir,fout);
    ////chi2viter
    //makeSMPInclJetXsec_NLOunf_chi2viter(outdir,fout);
    ////pearson matrix
    //makeSMPInclJetXsec_NLOunfpearsonmat_onePadOneEta(outdir,fout);
    
    //return;

    //////ratio comparisons with target PDF spectra //fine
    //for(int i=0; i<1;i++)//loop over unc opts
    //  for(int l=0; l<2;l++)	//loop over orders
    //	for(int k=0; k<2;k++){ //loop over scale opts
    //	  //if(k==1)continue;//skip murmufpt
    //	  //if(k==1)continue;//skips the one scale choice
    //	  for(int j=0; j<3;j++){	     //loop over targ pdfs    
    //	    //makeSMPInclJetXsec_NLOunfdata_targPDF_ratios(outdir,fout, uncopts[i], targpdfs[j], scaleopts[k], orders[l]);
    //	    //makeSMPInclJetXsec_NLOunfdata_targPDF_ratios(outdir,fout, uncopts[i], targpdfs[j], scaleopts[k], orders[l], true);//use py8unf data instead
    //	    makeSMPInclJetXsec_NLOunfdata_targPDF_ratios(outdir,fout, uncopts[i], targpdfs[j], scaleopts[k], orders[l], false);
    //	  }
    //	}
    
    
    //return;
    
    
    //unf data spectra+ratio v. array of NNPDF3.1 PDFs using slightly diff alphaS values. last arg is "forJohn)
    for(int k=0; k<2;k++){ //loop over scale opts
      if(k==1)continue;//skip murmufpt
      //makeSMPInclJetXsec_NLOunfdata_NNPDFs_ratios(outdir,fout,  (std::vector<std::string>) NNPDFs_vect, scaleopts[k],  false);
      //makeSMPInclJetXsec_NLOunfdata_NNPDFs_ratios(outdir,fout,  (std::vector<std::string>) NNPDFs_vect, scaleopts[k],  true); 
      //makeSMPInclJetXsec_NLOunfdata_PDFEnsemble_ratios(outdir, fout, (std::vector<std::string>) PDFensemble_as118_vect, PDFensemble_as118_colors_vect, scaleopts[k], false, "0.118", PDFensemble_as118_denom);
      makeSMPInclJetXsec_NLOunfdata_PDFEnsemble_ratios(outdir, fout, (std::vector<std::string>) PDFensemble_as120_vect, PDFensemble_as120_colors_vect, scaleopts[k], false, "0.120", PDFensemble_as120_denom);
      
    } 
    
    
    
    
    
  }
  
  
  


  

 
  
  if(writeToFile){
    fout->Write();
    fout->Close();
  }
  return;
}










//deprecrated below here

    ////ratio comparison with target PDF spectra's NLO and NNLO spectra // fine
    //for(int ybin=0; ybin<4; ybin++){//loop over ybins
    //  for(int k=0; k<2;k++){ //loop over scale opts
    //  if(k==1)continue;
    //	for(int j=0; j<2;j++){	     //loop over targ pdfs
    //	  makeSMPInclJetXsec_allOunfdata_ratios(outdir,fout,  targpdfs[j], scaleopts[k], ybin);
    //	}      }    }     
    //return;
    
      //for(int ybin=0; ybin<4; ybin++){//loop over ybins
	//makeSMPInclJetXsec_NLOunfdata_NNPDFs(outdir,  fout,  (std::vector<std::string>) NNPDFs_vect, scaleopts[k], ybin, false);//fine
	//makeSMPInclJetXsec_NLOunfdata_NNPDFs(outdir,  fout,  (std::vector<std::string>) NNPDFs_vect, scaleopts[k], ybin, true);//fine
	//makeSMPInclJetXsec_NLOunfdata_NNPDFs_ratios_v0(outdir,fout,  (std::vector<std::string>) NNPDFs_vect, scaleopts[k], ybin, false); //good
	//makeSMPInclJetXsec_NLOunfdata_NNPDFs_ratios_v0(outdir,fout,  (std::vector<std::string>) NNPDFs_vect, scaleopts[k], ybin, true); // also presumably good


  
  ////////////QUICK COMPARISONS SECTION
  //spectra ratios to |y|<0.5, for NLO+NP theory, GEN PY8, and unf data
  //std::cout<<"running makeSMPInclJetXsec_NLOunfdataAcrossy_ratios"<<std::endl;
  //makeSMPInclJetXsec_NLOunfdataAcrossy_ratios(outdir,fout);
  //std::cout<<"running makeSMPInclJetXsec_PY8unfdataAcrossy_ratios"<<std::endl;
  //makeSMPInclJetXsec_PY8unfdataAcrossy_ratios(outdir,fout);


  //////comparisons with JOHNS unfolded results
  ////std::cout<<"running makeSMPInclJetXsec_JOHNNLOunfdata"<<std::endl;
  ////makeSMPInclJetXsec_JOHNNLOunfdata(outdir,fout);
  ////std::cout<<"running makeSMPInclJetXsec_JOHNNLOunfdata_ratios"<<std::endl;
  ////makeSMPInclJetXsec_JOHNNLOunfdata_ratios(outdir,fout);
  ////////comparisons with JOHNS measured data results
  ////std::cout<<"running makeSMPInclJetXsec_JOHNmeasdata"<<std::endl;
  ////makeSMPInclJetXsec_JOHNmeasdata(outdir,fout);
  ////std::cout<<"running makeSMPInclJetXsec_JOHNmeasdata_ratios"<<std::endl;
  ////makeSMPInclJetXsec_JOHNmeasdata_ratios(outdir,fout);  
  //////comparisons with JOHNS unfolded results
  ////std::cout<<"running makeSMPInclJetXsec_JOAONLOunfdata"<<std::endl;
  ////makeSMPInclJetXsec_JOAONLOunfdata(outdir,fout);
  ////std::cout<<"running makeSMPInclJetXsec_JOAONLOunfdata_ratios"<<std::endl;
  ////makeSMPInclJetXsec_JOAONLOunfdata_ratios(outdir,fout); 


  ////std::cout<<"running makeSMPInclJetXsec_PY8unfdata_DebugLaptop_ratios"<<std::endl;
  //makeSMPInclJetXsec_PY8unfdata_DebugLaptop_ratios(outdir,fout);
  //makeSMPInclJetXsec_PY8unfdata_DebugLaptop_ratios(outdir);


    //////ratio comparisons with thy/GEN/truth distributions used to unfold // error bands look slightly different than the ones produced by targPDF_ratios using equivalent inputs... look into
    //makeSMPInclJetXsec_NLOunfdata_ratios(outdir,fout,"");             //deprecate me
    //makeSMPInclJetXsec_NLOunfdata_ratios(outdir,fout,"totUnc");       //deprecate me
    //makeSMPInclJetXsec_NLOunfdata_ratios(outdir,fout,"datatotUnc");   //deprecate me
    //makeSMPInclJetXsec_NLOunfdata_ratios(outdir,fout,"NLOtotUnc");    //deprecate me
    //return;




    //NLOsyst v order
    //makeSMPInclJetXsec_NLOsyst_ratios(outdir,fout,"LO");   //deprecate me
    //makeSMPInclJetXsec_NLOsyst_ratios(outdir,fout,"NLO");  //deprecate me
    //makeSMPInclJetXsec_NLOsyst_ratios(outdir,fout,"NNLO"); //deprecate me
