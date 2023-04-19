const std::string etaORy_str="y";
//const bool etaORy_str="#eta";

//// ETA OR YBIN STRINGS FOR FINAL PLOTTING MACRO --> STILL CALLED 'etabin_strs' CAUSE IM LAZY
const std::string etabin_strs[]={
  "0.0 < #||{"+etaORy_str+"} < 0.5",  
  "0.5 < #||{"+etaORy_str+"} < 1.0",
  "1.0 < #||{"+etaORy_str+"} < 1.5",
  "1.5 < #||{"+etaORy_str+"} < 2.0"
};
const int n_etabin_strs=sizeof(etabin_strs)/sizeof(std::string);

//this array always ybin strs
const std::string const_ybin_strs[]={
  "0.0 < #||{y} < 0.5",  
  "0.5 < #||{y} < 1.0",
  "1.0 < #||{y} < 1.5",
  "1.5 < #||{y} < 2.0"
};
const int n_const_ybin_strs=sizeof(const_ybin_strs)/sizeof(std::string);

//this array always eta strs
const std::string const_etabin_strs[]={
  "0.0 < #||{#eta} < 0.5",  
  "0.5 < #||{#eta} < 1.0",
  "1.0 < #||{#eta} < 1.5",
  "1.5 < #||{#eta} < 2.0"
};
const int n_const_etabin_strs=sizeof(const_etabin_strs)/sizeof(std::string);


const std::string prelimpavetitle="#bf{CMS}                                                          #bf{27.4 pb^{-1} (5.02 TeV)}";
const std::string prelimratpavetitle="      #bf{CMS}                                                 #bf{27.4 pb^{-1} (5.02 TeV)}";
const std::string prelimth2pavetitle="      #bf{CMS}                                                 #bf{27.4 pb^{-1} (5.02 TeV)}";

const std::string simpavetitle="#bf{CMS} #it{Simulation}                                                      #bf{(5.02 TeV)}";
const std::string simratpavetitle="     #bf{CMS} #it{Simulation}                                              #bf{(5.02 TeV)}";
const std::string simth2pavetitle="     #bf{CMS} #it{Simulation}                                              #bf{(5.02 TeV)}";


const std::string yaxtitle="#frac{d^{2}#sigma}{d"+etaORy_str+" dp_{T}} [#frac{pb}{GeV}]";
const std::string xaxtitle="Jet p_{T} [GeV]";

const std::string sqrts="#sqrt{s} = 5.02 TeV";
const std::string jettype="anti-k_{T} PF Jets (R = 0.4)";
const std::string ptcuts="56 GeV < Jet p_{T} < 967 GeV";
const std::string ptcuts_lo="56 GeV < Jet p_{T} < ";
const std::string intlumi="L_{int} = 27.4 pb^{-1}";

const std::string scalechoice_murmufpt ="#mu_{r} = #mu_{f} = p_{T}^{Jet}";
const std::string scalechoice_murmufpt1="#mu_{r} = #mu_{f} = p_{T}^{Lead Jet}";
//const std::string scalechoice_murmufHTp="#mu_{r} = #mu_{f} = #Sigma p_{T}^{parton}";
const std::string scalechoice_murmufHTp="#mu_{r} = #mu_{f} = H_{T}";


std::vector<std::string> NNPDFs_vect={
  "NNPDF31_nnlo_as_0108",
  "NNPDF31_nnlo_as_0110",
  "NNPDF31_nnlo_as_0112",
  "NNPDF31_nnlo_as_0114",
  "NNPDF31_nnlo_as_0116",
  "NNPDF31_nnlo_as_0117", 
  "NNPDF31_nnlo_as_0118",
  "NNPDF31_nnlo_as_0119",
  "NNPDF31_nnlo_as_0120",
  "NNPDF31_nnlo_as_0122",
  "NNPDF31_nnlo_as_0124" 
};


std::vector<std::string> PDFensemble_as120_vect={
"ABMP16als120_5_nnlo",
"NNPDF23_nnlo_as_0120",
"NNPDF31_nnlo_as_0120"
//"CT10nnlo_as_0120", //problematic
//"CT14nnlo_as_0120", //problematic
//"CT18NNLO_as_0120", //problematic
//"HERAPDF20_NNLO_ALPHAS_120", //problematic

};

std::vector<Color_t> PDFensemble_as120_colors_vect={
  2,
  3,
  4//,
//  5,
//  6,
//  7,
//  8,
//  9,
};

const int PDFensemble_as120_denom=2;//NNPDF31_nnlo_as_120


std::vector<std::string> PDFensemble_as118_vect={
"ABMP16als118_5_nnlo",
"CT10nnlo",
"CT14nnlo",
"CT18NNLO",
"MSHT20nnlo_as118",
"NNPDF23_nnlo_as_0118",
"NNPDF30_nnlo_as_0118",
"NNPDF31_nnlo_as_0118"
//"CT10nnlo_as_0118", //problematic
//"CT14nnlo_as_0118", //problematic
//"CT18NNLO_as_0118", //problematic
//"HERAPDF20_NNLO_ALPHAS_118",//problematic
};

std::vector<Color_t> PDFensemble_as118_colors_vect={
  2,
  3,
  4,
  5,
  6,
  7,
  8,
  9,
};

const int PDFensemble_as118_denom=7;//NNPDF31_nnlo_as_118
