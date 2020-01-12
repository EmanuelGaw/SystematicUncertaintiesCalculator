#define public private
#include "TSystem.h"
#include "TObject.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TH1.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>

void drawHistAndProfComparsion(TH1* inputHistogram, vector<TH1*> outHistograms, TProfile * resultProfile);

void calculateSystematicUncertainties(string fileName) {
	ifstream file;
	file.open(fileName.c_str(), ios::in);

	if(file.is_open()) {
		vector<string> rootFiles;
		string line, histName;
		getline(file, line);
		bool symetrize = line == "true" ? true : false;
		getline(file, histName);

		while(getline(file, line, '\n')) {
			//printf(":%s:\n", line.c_str());
			stringstream ss(line);
			while(getline(ss, line, ' '))
				rootFiles.push_back(line);
		}
		
		file.close();

		unsigned rootFilesCnt = rootFiles.size();
		if(rootFilesCnt > 0) {
			
			vector<TH1*> systHists;
			vector<TH1*> systUpHists;
			vector<TH1*> systDownHists;
			vector<TFile*> tmpFiles;
			for (unsigned i = 0; i < rootFilesCnt; i++) {
				if(!gSystem->AccessPathName(rootFiles[i].c_str())) {
					//printf("%s\n", rootFiles[i].c_str());
					tmpFiles.push_back(new TFile(rootFiles[i].c_str()));
					
					if(rootFiles[i].find("up.root") != string::npos) systUpHists.push_back((TH1*)gFile->Get(histName.c_str()));
					else if(rootFiles[i].find("down.root") != string::npos) systDownHists.push_back((TH1*)gFile->Get(histName.c_str()));
					else systHists.push_back((TH1*)gFile->Get(histName.c_str()));
					//tmpFiles[i]->Close();
				}
			}
		
			unsigned systHistsCnt = systHists.size(), systUpHistsCnt = systUpHists.size(), systDownHistsCnt = systDownHists.size();

			if (systHistsCnt > 0) {
				TFile * result = new TFile("Output.root", "RECREATE");
				systHists[0]->Write("nominal");

				if(systUpHistsCnt != systDownHistsCnt) printf("Missing systXUp.root or systXDown.root to pair.\n");
				else if((systHistsCnt + systUpHistsCnt + systDownHistsCnt) > 1) {

					double nominal, syst, diff, center;
					double systUpTotal, systUp, systDownTotal, systDown;
					double systMeans[3], systUpMeans[3], systDownMeans[3];
					unsigned binCnt = systHists[0]->GetNbinsX();
					unsigned i, j, k;

					TGraphAsymmErrors * systematics[3];
					for(k = 0; k < 3; k++) {
						systematics[k] = new TGraphAsymmErrors();
						//systematics[k] = new TGraphAsymmErrors(systHists[0]);
					}

					for (i = 0; i <= binCnt; i++) {

						for(k = 0; k < 3; k++) {
							systMeans[k] = 0;
							systUpMeans[k] = 0;
							systDownMeans[k] = 0;
						}

						nominal = systHists[0]->GetBinContent(i);
						center = systHists[0]->GetBinCenter(i);
						
						for (j = 1; j < systHists.size(); j++) {
							syst = systHists[j]->GetBinContent(i);
							diff = nominal - syst;

							systMeans[0] += pow(diff, 2);
							systMeans[1] += pow(diff/nominal, 2);
							systMeans[2] += pow(syst/nominal, 2);
							
							if(diff < 0) {
								systUpMeans[0] += pow(diff, 2);
								systUpMeans[1] += pow(diff/nominal, 2);
								systUpMeans[2] += pow(syst/nominal, 2);
							} else {
								systDownMeans[0] += pow(diff, 2);
								systDownMeans[1] += pow(diff/nominal, 2);
								systDownMeans[2] += pow(syst/nominal, 2);
							}
						}
						
						//Handle pairs (systXUp,systXDown)
						for (j = 0; j < systUpHistsCnt; j++) {
							systUp = nominal - systUpHists[j]->GetBinContent(i);
							systUpTotal = 0;
							systDown = nominal - systDownHists[j]->GetBinContent(i);
							systDownTotal = 0;

							//up/down are negative:
							if(systUp*systDown > 0 && systUp < 0) {
								if(systUp > systDown) systUpTotal = systUp;
								else systUpTotal = systDown;
							}
							//up/down are positive:
							else if(systUp*systDown > 0 && systUp > 0) {
								if(systUp > systDown) systDownTotal = systDown;
								else systDownTotal = systUp;
							}
							//up/down are of different sign:
							else if(systUp*systDown < 0) {
								if(systUp > 0) {
									systUpTotal = systDown;
									systDownTotal = systUp;
								} else {
									systUpTotal = systUp;
									systDownTotal = systDown;
								}
							}

							systUpMeans[0] += pow(systUpTotal, 2);
							systUpMeans[1] += pow(systUpTotal/nominal, 2);
							systUpMeans[2] += pow(systUp/nominal, 2);
							
							systDownMeans[0] += pow(systDownTotal, 2);
							systDownMeans[1] += pow(systDownTotal/nominal, 2);
							systDownMeans[2] += pow(systDown/nominal, 2);

							//symetric case, choose bigger according to val
							if(abs(systUpTotal) >= abs(systDownTotal)){
							  systMeans[0] += pow(systUpTotal, 2);
							  systMeans[1] += pow(systUpTotal/nominal, 2);
							  systMeans[2] += pow(systUp/nominal, 2);
							} else {
							  systMeans[0] += pow(systDownTotal, 2);
							  systMeans[1] += pow(systDownTotal/nominal, 2);
							  systMeans[2] += pow(systDown/nominal, 2);
							}
						}

						for(k = 0; k < 3; k++, nominal = 0) {
							systMeans[k] = sqrt(systMeans[k]);
							systUpMeans[k] = sqrt(systUpMeans[k]);
							systDownMeans[k] = sqrt(systDownMeans[k]);
							
							systematics[k]->SetPoint(i, center, nominal);
							if(symetrize) systematics[k]->SetPointError(i, 0, 0, systMeans[k], systMeans[k]);
							else systematics[k]->SetPointError(i, 0, 0, systDownMeans[k], systUpMeans[k]);
						}
					}
					
					systematics[0]->Write("absSyst");
					systematics[1]->Write("relSyst");
					systematics[2]->Write("relSyst1");
					
					delete result;
				}
			} else printf("No .root files from configuration present in current dir.\n");
		} else printf("No .root files provided in configuration file.\n");
	} else printf("Cannot open the given file.\n");
}

void calculateUncorrelatedSystematicUncertainties(TH1 * inputHistogram, unsigned randCount = 10) {

	//gStyle->SetErrorX(0);
	//clock_t tStart = clock();

	TFile * result = new TFile("Output.root", "RECREATE");
	inputHistogram->SetDirectory(result);
	inputHistogram->SetOption("E1");
	inputHistogram->SetStats();
	inputHistogram->Write("inputHistogram");

	double mean, sigma, randomVal;
	unsigned i, j, binCnt = inputHistogram->GetNbinsX();
	const unsigned xbinxCnt = binCnt + 1;

	vector<TH1*> outHistograms;
	Double_t xbins[xbinxCnt];
	for (i = 0; i <= binCnt; i++) {
		xbins[i] = inputHistogram->GetBinLowEdge(i+1);
	}
	TProfile * resultProfile = new TProfile("Result","Result Profile", binCnt, xbins);

	for (i = 0; i < randCount; i++) { //printf("\n");

		outHistograms.push_back((TH1*) inputHistogram->Clone());
		outHistograms[i]->Reset();

		for (j = 0; j <= binCnt; j++) {
			mean = inputHistogram->GetBinContent(j);
		    sigma = inputHistogram->GetBinError(j);
			randomVal = gRandom->Gaus(mean,sigma);

			outHistograms[i]->SetBinContent(j, randomVal);
			outHistograms[i]->SetBinError(j, sigma);

			resultProfile->Fill(outHistograms[i]->GetBinCenter(j), randomVal);
            //printf("%g :%g: %g / ",resultProfile->GetBinCenter(j), randomVal, outHistograms[i]->GetBinCenter(j));
		}

		ostringstream histTitleStr;
		histTitleStr << "Result histogram no. "; histTitleStr << i; histTitleStr << ";X;Y(X)";
		outHistograms[i]->SetTitle(histTitleStr.str().c_str());
		outHistograms[i]->SetLineColor(51);
		outHistograms[i]->SetMarkerStyle(20);
		outHistograms[i]->SetOption("E2");

		ostringstream objNameStr;
		objNameStr << "out"; objNameStr << i;
		outHistograms[i]->Write(objNameStr.str().c_str());
	}

	resultProfile->SetLineWidth(2);
	resultProfile->SetLineColor(92);
	resultProfile->SetOption("E1");
	resultProfile->Write("resultProfile");
	//outHistograms[i]->Fit("gaus", "Q");

	drawHistAndProfComparsion(inputHistogram, outHistograms, resultProfile);

	//delete resultProfile;
	delete result;

	//printf("Time taken: %.5gs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
}

void run(string fileName, unsigned randCount = 10) {
	TFile *f = new TFile(fileName.c_str());
	TH1 *hist = (TH1*)gFile->Get("h_ratio");
	calculateUncorrelatedSystematicUncertainties(hist, randCount);
	delete f;
}

void drawHistAndProfComparsion(TH1* inputHistogram, vector<TH1*> outHistograms, TProfile * resultProfile) {

	gStyle->SetErrorX(0);
	gStyle->SetOptStat(kFALSE); //remove info box from histogram
	//gPad->SetLogy(kTRUE);

	inputHistogram->SetTitle("InputHistogram + ResultProfile;X;Y(X)");
	inputHistogram->GetXaxis()->SetTitleOffset(inputHistogram->GetYaxis()->GetTitleOffset() + 0.5);
	inputHistogram->SetMarkerStyle(25);
	inputHistogram->SetLineWidth(4);
	inputHistogram->SetLineColor(kBlue - 3);

    TCanvas *presenter = new TCanvas("Panel", "Wynik", 800, 1000);
    presenter->Divide(1, 2);

	presenter->cd(1);
    inputHistogram->Draw("E1");

	resultProfile->SetLineColor(kGreen);
	resultProfile->SetLineWidth(4);
	//resultProfile->SetStats(kTRUE);

	//inputHistogram->SetLineWidth(2);
    resultProfile->Draw("E1 same");

	presenter->cd(2);

	unsigned outHistsCnt = outHistograms.size();
	ostringstream histTitle;
	histTitle << "InputHistogram + Out[0:"; histTitle << (outHistsCnt - 1); histTitle << "];X;Y(X)";
	outHistograms[0]->SetTitle(histTitle.str().c_str());
	outHistograms[0]->GetXaxis()->SetTitleOffset(outHistograms[0]->GetYaxis()->GetTitleOffset() + 0.5);
	outHistograms[0]->SetLineColor(kBlue);
	//outHistograms[0]->SetLineWidth(1);
	outHistograms[0]->Draw("E2");
	outHistograms[0]->SetDirectory(0);
	for(unsigned i = 1; i < outHistsCnt; i++) {
		outHistograms[i]->Draw("E2 same");
		outHistograms[i]->SetDirectory(0);
	}
	inputHistogram->Draw("E1 same");
	//presenter->Print("practice.jpg");
	inputHistogram->SetDirectory(0);
	resultProfile->SetDirectory(0);
	//gPad->Print("practice.jpg");

	//presenter->BuildLegend();
}