// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"


namespace Rivet {

  // This analysis is a derived from the class Analysis:
  class CMS_Dijet_Mass_13TeV_test : public Analysis {

 
  private:
    BinnedHistogram<double> _hist_Mjj;
    BinnedHistogram<double> _hist_Pt;
    Histo1DPtr _h_leading_y;
    Histo1DPtr _h_subleading_y;
    Histo1DPtr _h_mjj_y;

  public:
    // @name Constructors, init, analyze, finalize
    // @{

    // Constructor
    CMS_Dijet_Mass_13TeV_test()
      : Analysis("CMS_Dijet_Mass_13TeV_test") {
      //setNeedsCrossSection(true);
    }

    // Book histograms and initialize projections:
    void init() {
      
      const FinalState fs;

      // Initialize the projectors:
      addProjection(FastJets(fs, FastJets::ANTIKT, 0.4),"Jets");
	  
      int MassBins = 103;
      double MassBinning[104] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693,
	                         740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659,
				 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 
				 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000
				 };
		
      std::vector<double> MassHistobins;

      for(int i=0;i<MassBins+1;i++){
	 MassHistobins.push_back(MassBinning[i]);
      }	
	  
	      
		  
         
      const int neta = 5;
      const int nbins = 65;
      double etabins[neta+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
      double vx[neta][nbins] = 
		{	
		 {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389}, // Eta_0.0-0.5
		 {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3637, 5220, 5492, 0}, // Eta_0.5-1.0
		 {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2941, 3832, 4037, 0, 0, 0, 0, 0}, // Eta_1.0-1.5
		 {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2500, 2640, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // Eta_1.5-2.0
		 {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} // Eta_2.0-2.5
		
		}; 
   
     std::vector<vector<double> > tmpPtbins;
     std::vector<double> Ptbins;
		
     for (int i =0;i<neta;i++){
          
	     tmpPtbins.push_back(vector<double>());
 
		   for (int j=0;j<nbins && vx[i][j]!=0;j++){
 		
			tmpPtbins[i].push_back(vx[i][j]);
		    }
	     
	     Ptbins = tmpPtbins[i];
		 
             // Booking Histograms :
	     char Pt_histos[20];
             sprintf(Pt_histos, "%s%d%s","d0",i+1,"-Pt-AK4");
            _hist_Pt.addHistogram(etabins[i], etabins[i+1], bookHisto1D(Pt_histos,Ptbins));
		
	     char Mjj_histos[20];
             sprintf(Mjj_histos, "%s%d%s","d0",i+1,"-x01-y01-AK4");
	    _hist_Mjj.addHistogram(etabins[i], etabins[i+1], bookHisto1D(Mjj_histos,MassHistobins));
	
	}
    _h_leading_y = bookHisto1D(1, 1, 1);
    _h_subleading_y = bookHisto1D(2, 1, 1);
    _h_mjj_y = bookHisto1D(3, 1, 1);
    
	  
    }
    
    
	
	 // Analysis
    void analyze(const Event &event) {
		
	  
      const double weight = event.weight();      
      const Jets& jets = applyProjection<FastJets>(event, "Jets").jetsByPt(Cuts::pt>30.*GeV);
	  
	  //Calculating dijet mass and filling histograms
	  //CONTROL PLOTS: Leading jet Pt spectrum
	  
	  if ( jets.size() < 2 ) vetoEvent;
       
	   
			double ymaxdj = max(jets[0].momentum().absrapidity(), jets[1].momentum().absrapidity());
			double Mjj = FourMomentum(jets[0].momentum() + jets[1].momentum()).mass();
       
	   
	_hist_Mjj.fill(ymaxdj, Mjj, weight);
	_hist_Pt.fill(ymaxdj,jets[0].momentum().pT() / GeV, weight);
	
	/// Extra control plots...
	_h_leading_y->fill(jets[0].momentum().rapidity(), weight);
	_h_subleading_y->fill(jets[1].momentum().rapidity(), weight);
	_h_mjj_y->fill(ymaxdj, weight);
	
	  
    }
	 
	 
	
	// Finalize
    void finalize() {
       cout<<"cross Section: "<<crossSection()<<endl;
      
      _hist_Mjj.scale(crossSection()/sumOfWeights()/2, this);
      _hist_Pt.scale(crossSection()/sumOfWeights()/2, this);
      
     /* _h_leading_y.scale(crossSection()/sumOfWeights()/2, this);
      _h_subleading_y.scale(crossSection()/sumOfWeights()/2, this);
      _h_mjj_y.scale(crossSection()/sumOfWeights()/2, this);
      */
	  
    }
 };

  // This global object acts as a hook for the plugin system. 
  AnalysisBuilder<CMS_Dijet_Mass_13TeV_test> plugin_CMS_Dijet_Mass_13TeV_test;

}
    

    
