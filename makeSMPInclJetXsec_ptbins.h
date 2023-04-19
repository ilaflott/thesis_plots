

//-----------------------------------------------------------
//-----------------------------------------------------------
//-----------------------------------------------------------
//-----------------------------------------------------------
//-----------------------------------------------------------

std::vector<std::vector<double>> merged_SMP_ptbins{
  {//00eta05
    //    43., 49.,
    56.,
      64.,       74.,       84.,       97.,       114.,
      133.,      153.,      174.,      196.,
      220.,
      245.,
      272.,
      300.,
      330.,      //362.,
      395.,      //430.,      //468.,
      507.,      //548.,      //592.,
      638.,      //686.,      //737.,      // 790.,
      846.},

    {//05eta10
      //43., 49.,
      56.,
	64.,       74.,       84.,       97.,       114.,
	133.,      153.,      174.,      196.,
	220.,
	245.,
	272.,
	300.,
	330.,	//362.,
	395.,	//430.,	//468.,
	507.,	//548.,	//592.,
	638.,	//686.,	//737.,	// 790.,
	846.}, //WARNING empty bin in data 846. - 905.
      //905.,  //last entry in data spectra
      //967.,


  {//10eta15
    //43., 49.,
    56.,
      64.,       74.,       84.,       97.,       114.,
      133.,      153.,      174.,      196.,
      220.,
      245.,
      272.,      //300.,
      330.,      //362.,
      395.,      //430.,      //468.,
      507.,      //548.,      //592.,
      638.},//,
    //686.,  //last entry in data spectra      686-737
    //737.,
    //790.,

    {//15eta20
      //43., 49.,
      56.,
	64.,       74.,       84.,       97.,       114.,
	133.,      153.,      174.,      196.,      
	220.,  
	245.,
	272.,	//300.,
	330.,	//362.,
	395.,	//430.,	//468.,
	507.}//,

};

std::vector<std::vector<double>> default_SMP_ptbins{
  {//00eta05

    ////if i'm allowed to omit some bins at high pt, then i like this one, because it works for 05y10 as well as this one.
    //    43., 49.,
    56.,
      64.,       74.,       84.,       97.,       114.,
      133.,      153.,      174.,      196.,
      220.,
      245.,
      272.,
      300.,
      330.,
      362.,
      395.,
      430.,
      468.,
      507.,
      548.,
      592.,
      638.,
      686.,
      737.,
      790.,
      846.,
      905.,
      967.,
      1032.},
    //1101.},
    
    {//05eta10
      //if i'm allowed to omit some bins at high pt, and i care about preserving the bin edges used in 00y05, this is what i get
      //43., 49.,
      56.,
	64.,       74.,       84.,       97.,       114.,
	133.,      153.,      174.,      196.,
	220.,      245., 272.,
	300.,
	330.,
	362.,
	395.,
	430.,
	468.,
	507.,
	548.,
	592.,
	638.,
	686.,
	737.,
	790.,
	846.}, //empty bin 846-905
      //905.,     
    //967.},   
    //1032.,
    //1101.},

    
      {//10eta15
	//43., 49.,
	56.,
	  64.,       74.,       84.,       97.,       114.,
	  133.,      153.,      174.,      196.,
	  220.,
	  245.,
	  272.,
	  300.,
	  330.,
	  362.,
	  395.,
	  430.,
	  468.,
	  507.,
	  548.,
	  592.,
	  638.,
	  686.,
	  737.},
      //	  790.},
	//846., 
	//905.,     
	//967.,   
	//1032.,
	//1101.},
      
      
      
	{//15eta20
	  //43., 49.,
	  56.,
	    64.,       74.,       84.,       97.,       114.,
	    133.,      153.,      174.,      196.,      220.,      245.,
	    272.,
	    300.,
	    330.,
	    362.,
	    395.,
	    430.,
	    468.,
	    507.,
	    548.,
	    592.}//,//empty bin, 592-638
	//638.,
	//   686.,
	//  737.}
	  //790.}
	  //846., 
	  //905.,     
	  //967.,   
	  //1032.,
	  //1101.},

};








double defbins_00eta05[]={
  56.,
  64.,
  74.,
  84.,
  97.,
  114.,
  133.,
  153.,
  174.,
  196.,
  220.,
  245.,
  272.,
  300.,
  330.,
  362.,
  395.,
  430.,
  468.,
  507.,
  548.,
  592.,
  638.,
  686.//,//generally cut off here
  //737.,
  //790.,
  //846.,
  //905.,
  //967.,
  //1032.//,
  //  1101.
};
const int defbins_00eta05_nbins=sizeof(defbins_00eta05)/sizeof(double)-1;

double defbins_05eta10[]={
  56.,
  64.,
  74.,
  84.,
  97.,
  114.,
  133.,
  153.,
  174.,
  196.,
  220.,
  245.,
  272.,
  300.,
  330.,
  362.,
  395.,
  430.,
  468.,
  507.,
  548.,
  592.,
  638.,
  686.//,//generally cut off here
  //737.,
  //790.,
  //846.,
  //905.,
  //967.,
  //1032.
};
const int defbins_05eta10_nbins=sizeof(defbins_05eta10)/sizeof(double)-1;

double defbins_10eta15[]={


  56.,
  64.,
  74.,
  84.,
  97.,
  114.,
  133.,
  153.,
  174.,
  196.,
  220.,
  245.,
  272.,
  300.,
  330.,
  362.,
  395.,
  430.,
  468.,
  507.,
  548.,
  592.//,//generally cut off here
//  638.,
//  686.,
//  737.,
//  790.,
//  846.,
//  905.,
//  967.,
//  1032.
};
const int defbins_10eta15_nbins=sizeof(defbins_10eta15)/sizeof(double)-1;

double defbins_15eta20[]={
  56.,
  64.,
  74.,
  84.,
  97.,
  114.,
  133.,
  153.,
  174.,
  196.,
  220.,
  245.,
  272.,
  300.,
  330.,
  362.,
  395.,
  430.,
  468.,
  507.//,//generally cut off here
  //548.,
  //592.,
  //638.,
  //686.,
  //737.,
  //790.,
  //846.,
  //905.,
  //967.,
  //1032.
};
const int defbins_15eta20_nbins=sizeof(defbins_15eta20)/sizeof(double)-1;



double johnsbins_00eta05[]={
  56.,
  64.,  84.,  114.,
  153.,  196.,  245.,  300.,
  362.,  430.,  507.,  592.,
  686.,
  790.,  905.,  1032.
};
const int johnsbins_00eta05_nbins=sizeof(johnsbins_00eta05)/sizeof(double)-1;
double johnsbins_05eta10[]={
  56.,
  64.,  84.,  114.,
  153.,  196.,  245.,  300.,
  362.,  430.,  507.,  592.,
  686.,
  790.,  905.,  1032.
//  //56.,
//  64.,  84.,  114.,
//  153.,  196.,  245.,  300.,
//  362.,  430.,  507.,  592.,
//  686.//  ,
//      //  790.,  905.
};
const int johnsbins_05eta10_nbins=sizeof(johnsbins_05eta10)/sizeof(double)-1;
double johnsbins_10eta15[]={
  56.,
  64.,  84.,  114.,
  153.,  196.,  245.,  300.,
  362.,  430.,  507.,  592.,
  686.,
  790.,  905.,  1032.
//  //56.,
//  64.,  84.,  114.,
//  153.,  196.,  245.,  300.,
//  362.,  430.,  507.,  592. //,
//  //  686.
};
const int johnsbins_10eta15_nbins=sizeof(johnsbins_10eta15)/sizeof(double)-1;
double johnsbins_15eta20[]={
  56.,
  64.,  84.,  114.,
  153.,  196.,  245.,  300.,
  362.,  430.,  507.,  592.,
  686.,
  790.,  905.,  1032.
//  //56.,
//  64.,  84.,  114.,
//  153.,  196.,  245.,  300.,
//  362.,  430.,  507.  //,
//  //592.
};
const int johnsbins_15eta20_nbins=sizeof(johnsbins_15eta20)/sizeof(double)-1;






double johnsbins2_00eta05[]={
  56.,
  64.,  84.,  114.,
  153.,  196.,  245.,  300.,
  362.,  430.,  507.,  592.,
  686.//,
  //790.,  905.,  1032.
};
const int johnsbins2_00eta05_nbins=sizeof(johnsbins2_00eta05)/sizeof(double)-1;
double johnsbins2_05eta10[]={
  56.,
  64.,  84.,  114.,
  153.,  196.,  245.,  300.,
  362.,  430.,  507.,  592.,
  686.//  ,
      //  790.,  905.
};
const int johnsbins2_05eta10_nbins=sizeof(johnsbins2_05eta10)/sizeof(double)-1;
double johnsbins2_10eta15[]={
  56.,
  64.,  84.,  114.,
  153.,  196.,  245.,  300.,
  362.,  430.,  507.,  592. //,
  //  686.
};
const int johnsbins2_10eta15_nbins=sizeof(johnsbins2_10eta15)/sizeof(double)-1;
double johnsbins2_15eta20[]={
  56.,
  64.,  84.,  114.,
  153.,  196.,  245.,  300.,
  362.,  430.,  507.  //,
  //592.
};
const int johnsbins2_15eta20_nbins=sizeof(johnsbins2_15eta20)/sizeof(double)-1;





std::vector<std::vector<double>> johnsbins_semifin={
{
64., 84., 114., 153., 196., 245., 300., 395., 507., 638., 846.
},

{
64., 84., 114., 153., 196., 245., 300., 395., 507., 638., 846.
},

{
64., 84., 114., 153., 196., 245., 330., 395., 507., 638.
},

{
64., 84., 114., 153., 196., 245., 330., 395., 507.
}
};
