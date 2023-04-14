/// \brief Simple macro for ITS Standalone performances

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <array>

#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TProfile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>

#include <TMath.h>
#include <TMatrixD.h>

#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"
#include "DetectorsBase/Propagator.h"
#include "ITSBase/GeometryTGeo.h"
#include "SimulationDataFormat/TrackReference.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

#endif

Float_t Bz = -5;

extern TROOT* gROOT;

const int l0 = 6;
const int l1 = 5;
const int l2 = 4;

const int l3 = 3;
const int l4 = 2;
const int l5 = 1;
const int l6 = 0;

using GIndex = o2::dataformats::VtxTrackIndex;

void extractITSPerformance(){
  float vr = 0.1;

  //output histogram definition
  auto hd = new TH2D("hd", "Impact parameter XY vs phi;#phi (rad);D (cm)", 100, 0., 6.28, 100, -vr, vr);

  Int_t nb = 100;
  Double_t xbins[nb + 1], ptcutl = 0.1, ptcuth = 10.;
  Double_t a = TMath::Log(ptcuth / ptcutl) / nb;
  for (Int_t i = 0; i <= nb; i++) {
    xbins[i] = ptcutl * TMath::Exp(i * a);
  }
  auto hdpt = new TH2D("hdpt", "Impact parameter XY vs pt;pt (GeV/c);D (cm)", nb, xbins, 100, -vr, vr);
  auto hzpt = new TH2D("hzpt", "Impact parameter Z vs pt;pt (GeV/c);D (cm)", nb, xbins, 100, -vr, vr);

  auto hdT = new TH1D("hdT", "Impact parameter XY Top;D_{XY} (cm)", 100, -vr, vr);
  hdT->SetLineColor(2);
  auto hdB = new TH1D("hdB", "Impact parameter XY Bottom;D (cm)", 100, -vr, vr);

  auto hzT = new TH1D("hzT", "Impact parameter Z Top; D_{Z} (cm)", 100, -vr, vr);
  hzT->SetLineColor(2);
  auto hzB = new TH1D("hzB", "Impact parameter Z Bottom", 100, -vr, vr);

  // open the file for vertices
  TChain* pvTree = new TChain("o2sim");
  pvTree->AddFile("o2_primary_vertex.root");
  if(!pvTree) {
    printf("PV tree not found"); return;
  }

  // Reconstructed tracks
  TFile file("o2trac_its.root", "OPEN");
  TChain* trackTree = new TChain("o2sim");
  trackTree->AddFile("o2trac_its.root");
  if(!trackTree) {
    printf("TRACK tree not found"); return;
  }

  if (pvTree && trackTree) {
    //reading primary vertex file
    auto br = pvTree->GetBranch("PrimaryVertex");
    std::vector<o2::dataformats::PrimaryVertex>* vertices = nullptr;
    br->SetAddress(&vertices);

    // this referes to actual tracks
    auto indexbr = pvTree->GetBranch("PVTrackIndices");
    std::vector<GIndex>* vertexTrackIDs = nullptr;
    indexbr->SetAddress(&vertexTrackIDs);

    // this makes the connection of vertex to track indices
    auto v2totrackrefbr = pvTree->GetBranch("PV2TrackRefs");
    std::vector<o2::dataformats::VtxTrackRef>* v2trackref = nullptr;
    v2totrackrefbr->SetAddress(&v2trackref);

    //defines track branch in o2track file
    auto trackBranch = trackTree->GetBranch("ITSTrack");
    std::vector<o2::its::TrackITS>* itstracks= nullptr; 
    trackBranch->SetAddress(&itstracks);
    
    for(int ientry=0; ientry<pvTree->GetEntries(); ientry++){
      br->GetEntry(ientry);
      indexbr->GetEntry(ientry);
      v2totrackrefbr->GetEntry(ientry);
      trackBranch->GetEntry(ientry);

      if (vertices && vertexTrackIDs) {
        //printf("BLA: vertices size %d\n", vertices->size());
        int index = 0;
        for (auto& v : *vertices) {
          int BCid = 0;
          auto& cov = v.getCov();
          auto& ts = v.getTimeStamp();

          // now go over tracks via the indices
          auto& trackref = (*v2trackref)[index];
          int start = trackref.getFirstEntryOfSource(0);
          int ntracks = trackref.getEntriesOfSource(0);
          //printf("Start: %d, track vertex %d",start,ntracks);
          for (int ti = 0; ti < ntracks; ++ti) {
            auto trackindex = (*vertexTrackIDs)[start + ti];

            //fetching the actual tracks
            const auto source = trackindex.getSource();
            o2::its::TrackITS* trackIts = nullptr;
            trackIts = &((*itstracks)[trackindex.getIndex()]);
            //take track interesting variables
            std::array<float, 3> pxpypz;
            trackIts->getPxPyPzGlo(pxpypz);
            Float_t ip[2]{0., 0.};
            trackIts->getImpactParams(v.getX(), v.getY(), v.getZ(), 5., ip);
            hdpt->Fill(trackIts->getPt(), ip[0]);
            hzpt->Fill(trackIts->getPt(), ip[1]);
            if (pxpypz[1] > 0) {
                hdT->Fill(ip[0]);
                hzT->Fill(ip[1]);
            } else {
                hdB->Fill(ip[0]);
                hzB->Fill(ip[1]);
            }
            double phi = TMath::ATan2((double)pxpypz[1], (double)pxpypz[0]);
            if (phi < 0) phi += 6.28318;
                hd->Fill(phi, ip[0]);
            } //end track per PVvertex 
            index++;
        }//end vertices
      } //end if vertices
    } //loop over tree entries
  } // end tree inspection
  
   //output plots
  TFile *fout=new TFile("itsStandalonePerformances.root", "RECREATE");
  auto c = new TCanvas;
  c->SetGrid();
  hd->Draw();
  auto hd_fpx = hd->ProfileX();
  hd_fpx->SetMarkerColor(2);
  hd_fpx->SetMarkerStyle(20);
  hd_fpx->Draw("same");
  hd_fpx->Write();

  c = new TCanvas();//"c", "hdB", 0.,0.,900,700);;
  c->SetGrid();
  hdT->Draw();
  hdB->Draw("same");
  hdB->Write();

  c = new TCanvas();//"c", "hzB", 0.,0.,900,700);;
  c->SetGrid();
  hzT->Draw();
  hzB->Draw("same");
  hzB->Write();

  c = new TCanvas();//"c", "hdpt", 0.,0.,900,700);;
  c->SetGrid();
  c->SetLogx();
  hdpt->Draw("col");
  hdpt->Write();

  c = new TCanvas();//"c", "hdpt_2", 0.,0.,900,700);
  c->SetGrid();
  c->SetLogx();
  hdpt->FitSlicesY();
  auto hdpt_2 = (TH1D*)gROOT->FindObject("hdpt_2");
  hdpt_2->SetMaximum(0.03);
  hdpt_2->SetMinimum(0.00);
  hdpt_2->Draw();
  hdpt_2->Write();

  c = new TCanvas();//"c", "zpt_2", 0.,0.,900,700);
  c->SetGrid();
  c->SetLogx();
  hzpt->FitSlicesY();
  auto hzpt_2 = (TH1D*)gROOT->FindObject("hzpt_2");
  hzpt_2->SetMaximum(0.03);
  hzpt_2->SetMinimum(0.00);
  hzpt_2->Draw();
  hzpt_2->Write();

  //gDirectory->Write();
  fout->Close();

}
