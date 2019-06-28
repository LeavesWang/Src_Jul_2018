#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <unistd.h>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TStyle.h"

using namespace std;

struct StrtMesytec
{
	int modID;
	int data[32];
	int modRes;
	int modEC_TS;
};

struct StrtS800
{
	int tS;
	int eC;
	int trig;
	int tof[8]; //[0]-->[1]: ORTEC TAC+Phillips ADC; [2]-->[7]: Phillips TDC
	int crdcCath[2][5][64]; //[2]: two CRDC; [5]: [0]---sample; [1]-->[4]---energy
	int crdcAnode[2][2]; // [0]: energy; [1]: time	
	int hodoEgy[32]; // 32 crystals
	int hodoTime;
	int pin[5];
	int mesyTDC[16];
};

struct StrtAna
{
	double tof[4]; //[4]: [0] TAC+ADC+clock; [1] TAC+ADC; [2]: regular CFD+TDC; [3] MCFD16+TDC
	double tD[4][8][8];
	double egyPMT[2][4]; //energies in 4 PMTs and each plastic
	double egyPla[2]; //[2]: [0] is energy in plastic at S800, [1] is energy in plastic at A1900
	double xPlaT[4][2]; //[2]: [0] for S800 plastic from time info, [1] for A1900 plastic from time info
	double yPlaT[4][2];
	double xPlaQ[2]; //[2]: [0] for S800 plastic from amp info, [1] for A1900 plastic from amp info
	double yPlaQ[2];
	double xMCP[2]; //two gain settings
	double yMCP[2]; //two gain settings
	double delE[5];
	double tke;	
	double beta[4];
	double gamma[4];
	double Z[4];
	double dZ[4];
	double brho[2];
	double AoQ[4][2];
	double Q[4][2];
	double ZmQ[4][2];
	double ZImQ[4][2];
	double A[4][2];
	double Araw[4][2];
	double Am2Q[4][2];
	double Am3Q[4][2];
	double Am2Z[4][2];
	double Am3Z[4][2];
	double dAm2Z[4][2];
	double dAm3Z[4][2];
	int numTime[4][2];
	int sigTime[4][2];
	int Zi[4];
	int Am2Zi[4][2];
	int Am3Zi[4][2];
};

struct StrtPid
{
	double tof;
	double beta;
	double gamma;
	double Z;
	double dZ;
	double AoQ;
	double Q;
	double ZmQ;
	double ZImQ;
	double A;
	double Araw;
	double Am2Q;
	double Am3Q;
	double Am2Z;
	double Am3Z;
	double dAm2Z;
	double dAm3Z;
	int Zi;
	int Am2Zi;
	int Am3Zi;
};

const int AdcPinLow=10, AdcPinUp=4096;
const int AdcTofLow=500, AdcTofUp=7680;
const int TdcTofLow[16]={12000, 12000, 12000, 12000, 35000, 29000, 36000, 36000, 10000, 10000, 10000, 10000, 32000, 12000, 1, 14500};
const int TdcTofUp[16]={22000, 22000, 22000, 22000, 55000, 49000, 56000, 56000, 20000, 20000, 20000, 20000, 52000, 32000, 1, 34500};
const int QdcTofLow[8]={729, 750, 754, 778, 770, 780, 760, 780};
const int QdcMcpLow[8]={759, 765, 802, 791, 752, 764, 765, 760}; //new pedestal values 
const int QdcUp=3840;

const double CALADC[12]={6.46209, 6.59645, 6.56230, 6.57185, 6.44156, 6.58265, 6.64827, 6.52219, 6.45537, 6.42844, 6.65406, 6.43436};  //unit: ps/ch
const double CALTDC=3.90625; //ps/ch
const double CALPIN[5][2]={{0,1}, {0,1}, {0,1}, {0,1}, {0,1}};
const double CALTKE[6]={81.071121, 1.071346, 0.662579, 3.013299, 2.826749, 0};

const double CALXMCP[2][10]={{0.733129,26.744387,-0.091781,1.043661,0.047598,9.192684,2.637526,-0.929438,2.056948,0.576781},{0.802060,26.063777,-0.897100,1.296354,1.163047,11.688516,3.208674,-1.230582,-2.736673,3.004569}}; //[0]+[1]*x+[2]*x*x+[3]*y+[4]*y*y+[5]*x*x*x+[6]*y*y*y+[7]*x*y+[8]*x*x*y+[9]*x*y*y
const double CALYMCP[2][10]={{3.652901,19.180574,1.578795,-1.716251,0.330541,11.410052,-0.641449,-0.958885,0.507911,5.328422}, {3.727687,18.762661,-0.510623,-1.588110,-0.511162,10.227921,-1.138502,0.227536,0.858179,4.114189}}; //[0]+[1]*y+[2]*y*y+[3]*x+[4]*x*x+[5]*y*y*y+[6]*x*x*x+[7]*x*y+[8]*y*y*x+[9]*x*x*y	
// const double CALXMCP[2][10]={{0,1,0,0,0,0,0,0,0,0}, {0,1,0,0,0,0,0,0,0,0}}; //[0]+[1]*x+[2]*x*x+[3]*y+[4]*y*y+[5]*x*x*x+[6]*y*y*y+[7]*x*y+[8]*x*x*y+[9]*x*y*y //for raw pos
// const double CALYMCP[2][10]={{0,1,0,0,0,0,0,0,0,0}, {0,1,0,0,0,0,0,0,0,0}}; //[0]+[1]*y+[2]*y*y+[3]*x+[4]*x*x+[5]*y*y*y+[6]*x*x*x+[7]*x*y+[8]*y*y*x+[9]*x*x*y //for raw pos

const double CALTOF[4][2]={{500,-0.001}, {513.269,-0.001}, {579.248,-0.001}, {500,-0.001}};
const double BRHO0=3.7211; //Tm
const double DISP=106.84; // mm/%
const double LOF=60.763; //m
const double CALZ[4][2]={{0,1},{1.0214,5.9613},{1.0616,5.9556},{0,1}};
// const double CALZ[4][2]={{0,1},{0,1},{0,1},{0,1}};

const double CALTOF_PID[2]={579.248,-0.001};
const double CALZ_PID[2]={1.0616,5.9556};

const int iLow[4]={6, 5, 4, 7};
const string sSet[2]={"PS_270_382", "RS_270_382"};

void Root2Ana()
{	
	gStyle->SetOptStat("nemri");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1);
	gStyle->SetCanvasDefH(1080);
	gStyle->SetCanvasDefW(1920);
	
	string sRoot, sAna;	
	StrtMesytec madc, mtdc, mqdcTOF, mqdcMCP;
	StrtS800 s800;
	StrtAna ana;
	StrtPid pid;
	
	double tPMT[8], timeDet[2];
	double calQdcMcp[2][4];
	
	long long iEnt, nEnt;
	int i, j, k, m, n, p;
	int iAna;
	int nQdcTof[2];
	int nGoodEvt, nGoodPin, nGoodMcp[2];
	bool goodPin[5];
	double b;
	string setting;
	int run;
	double xMcpRaw=0, yMcpRaw=0;
		
	int runMin, runMax, runNum;
	cout<<"Input minimum and maximum numbers of run: ";
	cin>>runMin>>runMax;
	
	sAna="/home/kailong/ExpData/Jul2018/AnaData/ana-run-"+to_string(runMin)+"--"+to_string(runMax)+".root";
	TFile *fAna=new TFile(sAna.c_str(), "RECREATE");
	TTree *tAna=new TTree("tAna", "tree for data analysis");

	tAna->Branch("setting", &setting);
	tAna->Branch("run", &run, "run/I");
	tAna->Branch("ana", &ana, "tof[4]/D:tD[4][8][8]/D:egyPMT[2][4]/D:egyPla[2]/D:xPlaT[4][2]/D:yPlaT[4][2]/D:xPlaQ[2]/D:yPlaQ[2]/D:xMCP[2]/D:yMCP[2]/D:delE[5]/D:tke/D:beta[4]/D:gamma[4]/D:Z[4]/D:dZ[4]/D:brho[2]/D:AoQ[4][2]/D:Q[4][2]/D:ZmQ[4][2]/D:ZImQ[4][2]/D:A[4][2]/D:Araw[4][2]/D:Am2Q[4][2]/D:Am3Q[4][2]/D:Am2Z[4][2]/D:Am3Z[4][2]/D:dAm2Z[4][2]/D:dAm3Z[4][2]/D:numTime[4][2]/I:sigTime[4][2]/I:Zi[4]/I:Am2Zi[4][2]/I:Am3Zi[4][2]/I");
	tAna->Branch("pid", &pid, "tof/D:beta/D:gamma/D:Z/D:dZ/D:AoQ/D:Q/D:ZmQ/D:ZImQ/D:A/D:Araw/D:Am2Q/D:Am3Q/D:Am2Z/D:Am3Z/D:dAm2Z/D:dAm3Z/D:Zi/I:Am2Zi/I:Am3Zi/I");
	
	ostringstream ssRun;
	for(runNum=runMin; runNum<=runMax; runNum++)
	{
		run=runNum;
		
		ssRun.str("");
		ssRun<<setw(4)<<setfill('0')<<runNum;
		
		sRoot="/home/kailong/ExpData/Jul2018/RootData/run-"+ssRun.str()+"-00.root";
		if(access("/home/kailong/ExpData/Jul2018/AnaData", F_OK)!=0)
            system("mkdir /home/kailong/ExpData/Jul2018/AnaData");
		printf("\n**********Now converting %s to %s!**********\n\n", sRoot.c_str(), sAna.c_str());
    	
		TFile *fRoot = new TFile(sRoot.c_str());
		if(fRoot->IsZombie())
		{
			cout<<"Error in opening "<<sRoot<<"!\n";
			continue;
		}
		
		TTree *tData;
		fRoot->GetObject("tData",tData);
		if(!tData)
    	{
    		cout<<"Error read the tree of tData!\n";
    		continue;
    	}
		
		double mcpGainMat[8][2];
		for(i=0; i<8; i++)
		{
			mcpGainMat[i][0]=0;
			mcpGainMat[i][1]=1;
		}
		
		if(runNum>=270)
		{
			tData->SetEstimate(-1);
			nEnt=tData->GetEntries();		
			double *errCh=new double[nEnt];
			for(iEnt=0; iEnt<nEnt; iEnt++)
				errCh[iEnt]=0.5;
			// string strCanv="gainMat_"+to_string(runNum);
			// TCanvas *cGain=new TCanvas(strCanv.c_str(), strCanv.c_str());
			// cGain->Divide(2,2);
			TGraphErrors *gr[4];
			for(i=0; i<4; i++)
			{
				j=iLow[i];
				// cGain->cd(i+1);
				string sCut="mqdcMCP.data["+to_string(i)+"]>"+to_string(QdcMcpLow[i])+"&&mqdcMCP.data["+to_string(i)+"]<3840&&mqdcMCP.data["+to_string(j)+"]>"+to_string(QdcMcpLow[j])+"&&mqdcMCP.data["+to_string(j)+"]<3840";
				
				string sDraw="mqdcMCP.data["+to_string(i)+"]-"+to_string(QdcMcpLow[i])+":mqdcMCP.data["+to_string(j)+"]-"+to_string(QdcMcpLow[j]);
				
				// printf("Now drawing {%s} of {%s} when {%s}\n\n", sDraw.c_str(), sRoot.c_str(), sCut.c_str());
				
				long long nData=tData->Draw(sDraw.c_str(), sCut.c_str(), "goff");
				if(nData<2)
					continue;
				double *highGain=tData->GetV1();
				double *lowGain=tData->GetV2();
				gr[i]=new TGraphErrors(nData, lowGain, highGain, errCh, errCh);
				// gr[i]->Draw("AP");
				// gr[i]->SetTitle(("high"+to_string(i)+"_vs_low"+to_string(iLow[i])).c_str());
				TFitResultPtr fitRes=gr[i]->Fit("pol1","SQ");
				int fitSt=fitRes;
				if(fitSt!=0&&fitSt!=4000)
					continue;
				mcpGainMat[j][0]=fitRes->Parameter(0);
				mcpGainMat[j][1]=fitRes->Parameter(1);
				// printf("%f %f\n",mcpGainMat[j][0],mcpGainMat[j][1]);
			}
			delete []errCh;
			// cGain->SaveAs(("/home/kailong/ExpData/Jul2018/Graphs/Charts/"+strCanv+".png").c_str());
			// cGain->Close();
			for(i=0; i<4; i++)
				if(!gr[i])
					delete gr[i];
			// delete cGain;
		}
		
		memset(&madc, 0, sizeof(madc));
		memset(&mtdc, 0, sizeof(mtdc));
		memset(&mqdcTOF, 0, sizeof(mqdcTOF));
		memset(&mqdcMCP, 0, sizeof(mqdcMCP));
		memset(&s800, 0, sizeof(s800));
		
		tData->SetBranchAddress("madc", &madc);
		tData->SetBranchAddress("mtdc", &mtdc);
		tData->SetBranchAddress("mqdcTOF", &mqdcTOF);
		tData->SetBranchAddress("mqdcMCP", &mqdcMCP);
		tData->SetBranchAddress("s800", &s800);	
		
		for(i=0; i<2; i++)
		{
			string sfSet="/home/kailong/ExpData/Jul2018/Src/runNum"+sSet[i]+".dat";
			ifstream fSet(sfSet.c_str());
			string sRead;
			bool isFound=false;
			while(getline(fSet, sRead))
				if( sRead.find(to_string(runNum)) != string::npos )
				{
					isFound=true;
					break;
				}
			if(isFound)
			{
				setting=sSet[i];
				break;
			}
			if(i==1&&!isFound)
				setting="other_270_382";
		}
		
		for(iEnt=0; iEnt<tData->GetEntries(); iEnt++)
		// for(iEnt=0; iEnt<1000; iEnt++)
		{
			tData->GetEntry(iEnt);
			if(s800.trig==1||s800.trig==16)
			{
				TRandom3 r(0);
				
				memset(&ana, 0, sizeof(ana));
				memset(&pid, 0, sizeof(pid));
				
				memset(nQdcTof, 0, sizeof(nQdcTof));
				for(i=0; i<8; i++)
					if(mqdcTOF.data[2*i+1]>QdcTofLow[i]&&mqdcTOF.data[2*i+1]<QdcUp)
					{
						j=i-i/4*4;
						nQdcTof[i/4]++;
						ana.egyPMT[i/4][j]=mqdcTOF.data[2*i+1]+r.Uniform(-0.5,0.5)-QdcTofLow[i];
					}
				if(nQdcTof[0]==4)
				{
					// ana.xPlaQ[2]=log(ana.egyPMT[0][0]*ana.egyPMT[0][1]/ana.egyPMT[0][2]/ana.egyPMT[0][3]);
					// ana.yPlaQ[2]=log(ana.egyPMT[0][1]*ana.egyPMT[0][2]/ana.egyPMT[0][0]/ana.egyPMT[0][3]);
					
					ana.xPlaQ[0]=log(ana.egyPMT[0][1]/ana.egyPMT[0][3]);
					ana.yPlaQ[0]=log(ana.egyPMT[0][2]/ana.egyPMT[0][0]);
					ana.egyPla[0]=(ana.egyPMT[0][0]+ana.egyPMT[0][1]+ana.egyPMT[0][2]+ana.egyPMT[0][3])/4;
				}
				if(nQdcTof[1]==4)
				{					
					ana.xPlaQ[1]=log(ana.egyPMT[1][1]/ana.egyPMT[1][3]);
					ana.yPlaQ[1]=log(ana.egyPMT[1][0]/ana.egyPMT[1][2]);
					ana.egyPla[1]=(ana.egyPMT[1][0]+ana.egyPMT[1][1]+ana.egyPMT[1][2]+ana.egyPMT[1][3])/4;
				}

				ana.tke=CALTKE[0];
				nGoodPin=0;
				memset(goodPin, 0, sizeof(goodPin));
				for(i=0; i<5; i++)
					if(s800.pin[i]>AdcPinLow&&s800.pin[i]<AdcPinUp)
					{
						goodPin[i]=true;
						nGoodPin++;
						
						ana.delE[i]=CALPIN[i][0]+CALPIN[i][1]*(s800.pin[i]+r.Uniform(-0.5,0.5));
						ana.tke+=CALTKE[i+1]*(s800.pin[i]+r.Uniform(-0.5,0.5));
					}
				if(nGoodPin>1&&goodPin[0])
				{
					memset(calQdcMcp, 0, sizeof(calQdcMcp));
					memset(nGoodMcp, 0, sizeof(nGoodMcp));
					if(runNum>=270)
					{
						for(i=0; i<4; i++)
						{
							if(mqdcMCP.data[i]>QdcMcpLow[i]&&mqdcMCP.data[i]<QdcUp)
							{
								nGoodMcp[0]++;
								calQdcMcp[0][i]=mqdcMCP.data[i]-QdcMcpLow[i]+r.Uniform(-0.5,0.5);
								
								nGoodMcp[1]++;
								calQdcMcp[1][i]=calQdcMcp[0][i];
							}
							m=iLow[i];
							if(mqdcMCP.data[i]>=QdcUp&&mqdcMCP.data[m]>QdcMcpLow[m]&&mqdcMCP.data[m]<QdcUp)
							{
								
								calQdcMcp[1][i]=mcpGainMat[m][0]+mcpGainMat[m][1]*(mqdcMCP.data[m]-QdcMcpLow[m]+r.Uniform(-0.5,0.5));
								if(calQdcMcp[1][i]>QdcUp-QdcMcpLow[i])
									nGoodMcp[1]++;
							}
						}
						
						for(i=0; i<2; i++)
							if(nGoodMcp[i]==4)
							{
								xMcpRaw=(calQdcMcp[i][0]+calQdcMcp[i][3]-calQdcMcp[i][1]-calQdcMcp[i][2])/(calQdcMcp[i][0]+calQdcMcp[i][1]+calQdcMcp[i][2]+calQdcMcp[i][3]);
								
								yMcpRaw=(calQdcMcp[i][1]+calQdcMcp[i][3]-calQdcMcp[i][0]-calQdcMcp[i][2])/(calQdcMcp[i][0]+calQdcMcp[i][1]+calQdcMcp[i][2]+calQdcMcp[i][3]);
								
								ana.xMCP[i]=CALXMCP[i][0]+CALXMCP[i][1]*xMcpRaw+CALXMCP[i][2]*pow(xMcpRaw,2)+CALXMCP[i][3]*yMcpRaw+CALXMCP[i][4]*pow(yMcpRaw,2)+CALXMCP[i][5]*pow(xMcpRaw,3)+CALXMCP[i][6]*pow(yMcpRaw,3)+CALXMCP[i][7]*xMcpRaw*yMcpRaw+CALXMCP[i][8]*pow(xMcpRaw,2)*yMcpRaw+CALXMCP[i][9]*xMcpRaw*pow(yMcpRaw,2);
								
								ana.yMCP[i]=CALYMCP[i][0]+CALYMCP[i][1]*yMcpRaw+CALYMCP[i][2]*pow(yMcpRaw,2)+CALYMCP[i][3]*xMcpRaw+CALYMCP[i][4]*pow(xMcpRaw,2)+CALYMCP[i][5]*pow(yMcpRaw,3)+CALYMCP[i][6]*pow(xMcpRaw,3)+CALYMCP[i][7]*yMcpRaw*xMcpRaw+CALYMCP[i][8]*pow(yMcpRaw,2)*xMcpRaw+CALYMCP[i][9]*yMcpRaw*pow(xMcpRaw,2);
							}
					}
					for(i=0; i<2; i++)
						ana.brho[i]=BRHO0*(1+ana.xMCP[i]/DISP/100);
					
					nGoodEvt=0;
					for(iAna=0; iAna<4; iAna++)
					{
						memset(tPMT, 0, sizeof(tPMT));
						memset(timeDet, 0, sizeof(timeDet));
						
						if(iAna==0) //For TAC+ADC+clock
						{
							for(i=0; i<8; i++)
							{
								k=i/4;
								if(madc.data[i]>AdcTofLow&&madc.data[i]<AdcTofUp)
								{
									ana.numTime[0][k]++;
									ana.sigTime[0][k]=10*ana.numTime[0][k]+(i+1);
									tPMT[i]=CALADC[i]*(madc.data[i]+r.Uniform(-0.5,0.5));
									timeDet[k]+=tPMT[i];
								}
							}
							for(i=0; i<7; i++)
								for(j=i+1; j<8; j++)
									if(abs(tPMT[i])>0&&abs(tPMT[j])>0)
									{
										ana.tD[0][i][j]=tPMT[i]-tPMT[j];
										ana.tD[0][j][i]=-ana.tD[0][i][j];
									}
							if(ana.numTime[0][0]>0&&ana.numTime[0][1]>0)
								ana.tof[0]=timeDet[1]/ana.numTime[0][0]-timeDet[0]/ana.numTime[0][1];
						}
						
						if(iAna==1) //For TAC+ADC
						{
							for(i=8; i<12; i++)
							{
								m=i-8;
								p=i-4;
								if(madc.data[i]>AdcTofLow&&madc.data[i]<AdcTofUp)
								{
									ana.numTime[1][0]++;
									ana.numTime[1][1]++;
									ana.sigTime[1][0]=10*ana.numTime[1][0]+(m+1);
									ana.sigTime[1][1]=10*ana.numTime[1][1]+(p+1);
									ana.tD[1][p][m]=CALADC[i]*(madc.data[i]+r.Uniform(-0.5, 0.5));
									ana.tD[1][m][p]=-ana.tD[1][p][m];
									ana.tof[1]+=ana.tD[1][p][m];
								}
							}
							if(ana.numTime[1][0]>0)
								ana.tof[1]/=ana.numTime[1][0];
						}
						
						if(iAna==2)
						{
							for(j=0; j<8; j++)
							{
								if(iAna==2)
								{
									k=j+1;
									p=j;
								}
								if(iAna==3)
								{
									k=2*j+17;
									p=j+8;
								}
								n=j/4;
								if(mtdc.data[k]>TdcTofLow[p]&&mtdc.data[k]<TdcTofUp[p])
								{
									ana.numTime[iAna][n]++;
									ana.sigTime[iAna][n]=10*ana.numTime[iAna][n]+(j+1);
									tPMT[j]=CALTDC*(mtdc.data[k]+r.Uniform(-0.5, 0.5));
									timeDet[n]+=tPMT[j];
								}
							}
							for(i=0; i<7; i++)
								for(j=i+1; j<8; j++)
									if(abs(tPMT[i])>0&&abs(tPMT[j])>0)
									{
										ana.tD[iAna][i][j]=tPMT[i]-tPMT[j];
										ana.tD[iAna][j][i]=-ana.tD[iAna][i][j];
									}
							if(ana.numTime[iAna][0]>0&&ana.numTime[iAna][1]>0)
								ana.tof[iAna]=timeDet[1]/ana.numTime[iAna][1]-timeDet[0]/ana.numTime[iAna][0];
						}
						
						if(ana.numTime[iAna][0]>0&&ana.numTime[iAna][1]>0)
						{
							ana.tof[iAna]=CALTOF[iAna][0]+CALTOF[iAna][1]*ana.tof[iAna];
							
							if(iAna!=1&&ana.numTime[iAna][0]==4)
							{
								// ana.xPlaT[iAna][0]=(tPMT[2]+tPMT[3]-tPMT[0]-tPMT[1]);
								// ana.yPlaT[iAna][0]=(tPMT[0]+tPMT[3]-tPMT[1]-tPMT[2]);	
								ana.xPlaT[iAna][0]=(tPMT[3]-tPMT[1]);
								ana.yPlaT[iAna][0]=(tPMT[0]-tPMT[2]);
							}
							if(iAna!=1&&ana.numTime[iAna][1]==4)
							{
								// ana.xPlaT[iAna][1]=(tPMT[4]+tPMT[7]-tPMT[5]-tPMT[6]);
								// ana.yPlaT[iAna][1]=(tPMT[6]+tPMT[7]-tPMT[4]-tPMT[5]);
								ana.xPlaT[iAna][1]=(tPMT[7]-tPMT[5]);
								ana.yPlaT[iAna][1]=(tPMT[6]-tPMT[4]);
							}

							b=LOF/ana.tof[iAna]/0.299792458;
							if(b>0&&b<1)
							{
								ana.beta[iAna]=b;
								ana.gamma[iAna]=1/sqrt(1-b*b);
														
								ana.Z[iAna]=sqrt( ana.delE[0] / (log(5930/(1/b/b-1))/b/b-1) );
								ana.Z[iAna]=CALZ[iAna][0]+CALZ[iAna][1]*ana.Z[iAna];
								ana.Zi[iAna]=TMath::Nint(ana.Z[iAna]);
								ana.dZ[iAna]=ana.Z[iAna]-ana.Zi[iAna];
								
								if(runNum==150||runNum==152||runNum==153||(nGoodMcp[1]==4&&nQdcTof[0]>0&&nQdcTof[1]>0))
								{
									for(k=0; k<2; k++)
									{
										ana.AoQ[iAna][k]=ana.brho[k]/ana.beta[iAna]/ana.gamma[iAna]*0.32184043;
										ana.Q[iAna][k]=ana.tke/(931.4940954*(ana.gamma[iAna]-1)*ana.AoQ[iAna][k]);
										ana.ZmQ[iAna][k]=ana.Z[iAna]-ana.Q[iAna][k];
										ana.ZImQ[iAna][k]=ana.Zi[iAna]-ana.Q[iAna][k];
										ana.A[iAna][k]=ana.AoQ[iAna][k]*ana.Q[iAna][k];
										ana.Araw[iAna][k]=ana.AoQ[iAna][k]*ana.Zi[iAna];
										ana.Am2Q[iAna][k]=ana.A[iAna][k]-2*ana.Q[iAna][k];
										ana.Am3Q[iAna][k]=ana.A[iAna][k]-3*ana.Q[iAna][k];
										ana.Am2Z[iAna][k]=ana.Araw[iAna][k]-2*ana.Zi[iAna];
										ana.Am3Z[iAna][k]=ana.Araw[iAna][k]-3*ana.Zi[iAna];
										ana.Am2Zi[iAna][k]=TMath::Nint(ana.Am2Z[iAna][k]);
										ana.Am3Zi[iAna][k]=TMath::Nint(ana.Am3Z[iAna][k]);
										ana.dAm2Z[iAna][k]=ana.Am2Z[iAna][k]-ana.Am2Zi[iAna][k];
										ana.dAm3Z[iAna][k]=ana.Am3Z[iAna][k]-ana.Am3Zi[iAna][k];
										// if(ana.Am3Z[iAna][k]<0)
											// ana.dAm3Z[iAna][k]+=1;
									}
									nGoodEvt++;
								}
							}
						}
					}
					
					if(nGoodEvt>1&&ana.numTime[2][0]==4&&ana.numTime[2][1]==4&&(runNum==150||runNum==152||runNum==153||(nGoodMcp[1]==4&&nQdcTof[0]==4&&nQdcTof[1]==4))) //standard condition  //begin to fill branch of pid						
					{
						pid.tof=CALTOF_PID[0]+CALTOF_PID[1]*(ana.tD[2][4][0]+ana.tD[2][5][1]+ana.tD[2][6][2]+ana.tD[2][7][3])/4;
						
						b=LOF/pid.tof/0.299792458;
						if(b>0&&b<1)
						{
							pid.beta=b;
							pid.gamma=1/sqrt(1-b*b);

							pid.Z=sqrt( ana.delE[0]/ (log(5930/(1/b/b-1))/b/b-1) );
							pid.Z=CALZ_PID[0]+CALZ_PID[1]*pid.Z;
							pid.Zi=TMath::Nint(pid.Z);
							pid.dZ=pid.Z-pid.Zi;
							pid.AoQ=ana.brho[1]/pid.beta/pid.gamma*0.32184;
							pid.Q=ana.tke/(931.4940954*(pid.gamma-1)*pid.AoQ);
							pid.ZmQ=pid.Z-pid.Q;
							pid.ZImQ=pid.Zi-pid.Q;
							pid.A=pid.AoQ*pid.Q;
							pid.Araw=pid.AoQ*pid.Zi;
							pid.Am2Q=pid.A-2*pid.Q;
							pid.Am3Q=pid.A-3*pid.Q;
							pid.Am2Z=pid.Araw-2*pid.Zi;
							pid.Am3Z=pid.Araw-3*pid.Zi;
							pid.Am2Zi=TMath::Nint(pid.Am2Z);
							pid.Am3Zi=TMath::Nint(pid.Am3Z);
							pid.dAm2Z=pid.Am2Z-pid.Am2Zi;
							pid.dAm3Z=pid.Am3Z-pid.Am3Zi;
							// if(pid.Am3Z<0)
								// pid.dAm3Z+=1;
							tAna->Fill();
							// cout<<"yes"<<endl;
						}
					}
				}
			}
		}//end of whole tree
		fRoot->Close();
		delete fRoot;
	}
	fAna->cd();
	tAna->Write();
	fAna->Close();
	delete fAna;
}//end of whole function


void StandaloneApplication(int argc, char** argv)
{
	Root2Ana();
}

int main(int argc, char** argv)
{
	TApplication app("ROOT Application", &argc, argv);
	StandaloneApplication(app.Argc(), app.Argv());
	// app.Run();
	return 0;
}
