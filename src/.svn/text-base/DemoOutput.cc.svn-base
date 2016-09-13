#include <TMTrackTrigger/TMTrackFinder/interface/DemoOutput.h>
#include <TMTrackTrigger/TMTrackFinder/interface/Settings.h>

#include <iostream>

//=== Create EDM collections of all stubs and of the subset of stubs that are on L1 tracks, stored in the HardwareStub class.

void DemoOutput::getStubCollection(const matrix<HTpair>&  mHtPairs,
	                           std::auto_ptr<HwStubCollection>& hwStubs, std::auto_ptr<HwStubCollection>& hwStubsOnTracks) 
{
  for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
    for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {

      const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);

      // Access Rphi Transform array
      const HTrphi& htRphi = htPair.getRphiHT();
      const matrix<HTcell>& htArray = htRphi.getAllCells();
      // Loop over cells in the array
      for(unsigned int j = 0 ; j< htArray.size2(); ++j ){
	for(unsigned int i = 0 ; i < htArray.size1(); ++i) {
	  const vector < const Stub* > stubs = htArray(i,j).stubs();
	  for(const Stub* st:stubs) {
	    Stub stub = *st;
	    // Digitize stub relative to this phi sector.
	    stub.digitize(iPhiSec);
	    const DigitalStub& digiStub = stub.digitalStub();

	    // Calculate bin in Hough transform array of this stub, in format expect by hardware
	    int mbin = i;
	    int cbin = j;
	    // Store stub in format expected by hardware. Variables needed differ slightly according to firmware version.
	    l1t::HardwareStub lstub;
	    const int junk = 99999;
	    if (settings_->firmwareType() == 1) { // "Daisy chain" firmware in use, so dphi & rho not used in hardware.
	      lstub = l1t::HardwareStub(iPhiSec, iEtaReg, 0u, 1u, junk               , mbin, cbin, 
					digiStub.iDigi_PhiS(), digiStub.iDigi_Rt(), digiStub.iDigi_Z(), junk                 , 
					digiStub.m_min(), digiStub.m_max(), digiStub.iDigi_LayerID() );
	    } else { // "Systolic" or "2-c-bin" firmware in use, so m_min & m_max not used in hardware.
	      lstub = l1t::HardwareStub(iPhiSec, iEtaReg, 0u, 1u, digiStub.iDigi_Rho(), mbin, cbin, 
					digiStub.iDigi_PhiS(), digiStub.iDigi_Rt(), digiStub.iDigi_Z(), digiStub.iDigi_Dphi(), 
					junk            , junk            , digiStub.iDigi_LayerID() );
	    }
	    // Add some floating point numbers for debug purposes.
            lstub.setFphi(digiStub.phi());
            lstub.setFphiS(digiStub.phiS());
            lstub.setFrT(digiStub.rt() ) ;

	    // Store stub.
	    hwStubs->push_back(lstub);
	    // Also store only stubs associated with reconstructed L1 tracks.
	    if (htArray(i,j).trackCandFound()) {
	      // An individual HT cell doesn't know if tracks were killed because the sector was too "busy".
	      // So check if this track was read out within TM period. The r-phi HT array object knows this.
	      bool survived;
	      if (settings_->busySectorKill()) { // Is option to kill tracks not read out in TM period enabled?
		survived = false;
		const vector<L1track2D>& tracks = htRphi.trackCands2D();  
		for (const L1track2D& trk: tracks) {
		  if (i ==  trk.getCellLocation().first && j == trk.getCellLocation().second) survived = true;
		}
	      } else {
		survived = true;
	      }
	      if (survived) hwStubsOnTracks->push_back(lstub);
	    }
	  }
	}
      }
    }     
  }
}

//=== Create EDM collectiosn of all stubs produced by all tracking particles, or by the subset of tracking particles used
//=== for the algorithmic tracking efficiency measurement. Store the results in the HardwareTrack class.

void DemoOutput::getTPstubCollection(const matrix<HTpair>&  mHtPairs, const matrix<Sector>& mSectors, const vector<TP>& vTPs,
		    	             std::auto_ptr<HwTrackCollection>& hwStubsOnTP, std::auto_ptr<HwTrackCollection>& hwStubsOnTPforAlgEff) {

  for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
    for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {

      const Sector& sector = mSectors(iPhiSec, iEtaReg);
      const HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);

      // Fill HardwareTrack trees with a digitize version of TP good for efficiency studies       
      for( const TP tp: vTPs ){
	if(tp.useForEff()){
	  l1t::HardwareTrack ltrack(tp.index(), iPhiSec, iEtaReg, tp.pt(), tp.eta(), tp.phi0(), tp.charge(), tp.pdgId());
	  const vector < const Stub* > tpStubs = tp.assocStubs();
	  std::vector<const Stub*> tpStubsInSector;
	  for(const Stub* st : tpStubs){
	    if (sector.inside( st )) {
	      Stub stub = *st;
	      tpStubsInSector.push_back(st);

  	      // Digitize stub relative to this phi sector. 
	      // N.B. TP includes some stubs that failed front-end electronics window cut, and so can't be digitized. Veto these.
	      if (stub.frontendPass()) {
		stub.digitize(iPhiSec);
		const DigitalStub& digiStub = stub.digitalStub();

		// Store stub in format expected by hardware. Variables needed differ slightly according to firmware version.
		l1t::HardwareStub lstub;
		const int junk = 99999;
		if (settings_->firmwareType() == 1) { // "Daisy chain" firmware in use, so dphi & rho not used in hardware.
		  lstub = l1t::HardwareStub(iPhiSec, iEtaReg, 0u, 1u, junk               , 0, 0, 
					    digiStub.iDigi_PhiS(), digiStub.iDigi_Rt(), digiStub.iDigi_Z(), junk                 , 
					    digiStub.m_min(), digiStub.m_max(), digiStub.iDigi_LayerID() );
		} else { // "Systolic" or "2-c-bin" firmware in use, so m_min & m_max not used in hardware.
		  lstub = l1t::HardwareStub(iPhiSec, iEtaReg, 0u, 1u, digiStub.iDigi_Rho(), 0, 0, 
					    digiStub.iDigi_PhiS(), digiStub.iDigi_Rt(), digiStub.iDigi_Z(), digiStub.iDigi_Dphi(), 
					    junk            , junk            , digiStub.iDigi_LayerID() );
		}

		ltrack.Fill(lstub);
	      }
	    }
	  }
	  if (tpStubsInSector.size() > 0) hwStubsOnTP->push_back(ltrack);
 
	  unsigned int numLayers = Utility::countLayers(settings_, tpStubsInSector);
	  unsigned int minLayers = (tp.pt() > settings_->minPtToReduceLayers()) ? settings_->minStubLayers() : (settings_->minStubLayers() - 1 );

	  if (tp.useForAlgEff() && numLayers >= minLayers) hwStubsOnTPforAlgEff->push_back(ltrack);
	}
      }
    }		
  }
}


