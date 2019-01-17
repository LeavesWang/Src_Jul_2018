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
	double tof[4];
	double amp[4];
	double tD[4][8][8];
	double egy[4][8];
	double xPla[4][2]; //two plastic
	double yPla[4][2]; //two plastic
	double xMCP[4][2]; //two gain settings
	double yMCP[4][2]; //two gain settings
	double delE[4][5];
	double tke[4];	
	double beta[4];
	double gamma[4];
	double Z[4];
	double dZ[4];
	double brho[4][2];
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
	int Zi[4];
	int sig[4][2][2];
	int run;
};

struct StrtPid
{
	double tof;
	double xMCP[2]; //two gain settings
	double yMCP[2]; //two gain settings
	double delE[5];
	double tke;	
	double beta;
	double gamma;
	double Z;
	double dZ;
	double brho[2];
	double AoQ[2];
	double Q[2];
	double ZmQ[2];
	double ZImQ[2];
	double A[2];
	double Araw[2];
	double Am2Q[2];
	double Am3Q[2];
	double Am2Z[2];
	double Am3Z[2];
	double dAm2Z[2];
	double dAm3Z[2];
	int Zi;
	int run;
};

void Root2Ana()
{	
	gStyle->SetOptStat("nemri");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1);
	gStyle->SetCanvasDefH(1080);
	gStyle->SetCanvasDefW(1920);
	
	const int AdcPinLow=10, AdcPinUp=4096;
	const int AdcTofLow=500, AdcTofUp=7680;
	const int TdcTofLow[16]={12000, 12000, 12000, 12000, 1, 1, 1, 1, 10000, 10000, 10000, 10000, 1, 1, 1, 1};
	const int TdcTofUp[16]={22000, 22000, 22000, 22000, 65536, 65536, 65536, 65536, 20000, 20000, 20000, 20000, 65536, 65536, 65536, 65536};
	const int QdcTofLow[8]={800, 800, 800, 800, 780, 800, 770, 800};
	// const int QdcMcpLow[8]={1000, 760, 810, 850, 790, 790, 780, 790}; //old pedestal values
	const int QdcMcpLow[8]={759, 765, 802, 791, 752, 764, 765, 760}; //new pedestal values 
	const int QdcUp=3840;
	
	const double CALADC[12]={6.46209, 6.59645, 6.56230, 6.57185, 6.44156, 6.58265, 6.64827, 6.52219, 6.45537, 6.42844, 6.65406, 6.43436};  //unit: ps/ch
	const double CALTDC=3.90625; //ps/ch
	const double CALPIN[6][2]={{0,0.6951}, {0,0.6558}, {0,2.9832}, {0,2.7269}, {0,2.9703}, {0,0.4886}}; //Mev/ch  //0.6951 is the original slope of PIN0 and 0.4886 is related to the material in front of Si detectors.
	// const double CALXMCP[2][4]={{-2.62024, 27.6449, -2.11865, 1.01504}, {0,1,0,0}}; //mm, mm/chx, mm/chx^2, mm/chy
	// const double CALYMCP[2][4]={{2.99088, 18.4523, -4.75117, -1.07836}, {0,1,0,0}}; //mm, mm/chy, mm/chy^2, mm/chx
	const double CALXMCP[2][4]={{0,1,0,0}, {0,1,0,0}}; //mm, mm/chx, mm/chx^2, mm/chy
	const double CALYMCP[2][4]={{0,1,0,0}, {0,1,0,0}}; //mm, mm/chy, mm/chy^2, mm/chx
	const double CALTOF[4][2]={{500,-0.001}, {500,-0.001}, {577.648,-0.001}, {500,-0.001}}; //ns, ns/ps
	const double BRHO0=3.7211; //Tm
	const double DISP=112; // mm/%
	const double LOF=60.740; //m
	const double CALZ[2]={0, 1};
	
	const int QdcMcpLow_Pid[8]={726, 730, 745, 742, 760,760,760,760};
	const double CALPIN_PID[6][2]={{0,0.6951}, {0,0.6558}, {0,2.9832}, {0,2.7269}, {0,2.9703}, {0,0.4886}}; //Mev/ch  //0.6951 is the original slope of PIN0 and 0.4886 is related to the material in front of Si detectors.
	const double CALXMCP_PID[2][4]={{-4.11869, -26.62532, -3.38656, -19.399}, {0,1,0,0}}; //mm, mm/ch, mm/ch^2, mm/ch^3
	const double CALYMCP_PID[2][4]={{0,1,0,0}, {0,1,0,0}};
	const double CALTOF_PID[2]={570.4559, -0.001}; //ns, ns/ps
	const double BRHO0_PID=3.7211; //Tm
	const double DISP_PID=112; // mm/%
	const double LOF_PID=60.74; //m
	const double CALZ_PID[2]={1.191, 5.9377};
	
	string sSet[2]={"PS_270_382", "RS_270_382"};
	
	string sRoot, sAna;
	StrtMesytec madc, mtdc, mqdcTOF, mqdcMCP;
	StrtS800 s800;
	StrtAna ana;
	StrtPid pid;

	double tPMT[8];
	double timeDet[2], egyDet[2];
	double calQdcMcp[2][4];
	
	long long iEnt, nEnt;
	int i, j, k, m, n, p, q;
	int iAna;
	int coin[3];
	int nGoodEvt, nGoodPin, nGoodMcp[2];
	bool goodPin[5], goodPid;
	double b;
	string setting;
	double xMcpRaw, yMcpRaw;
	
	int iLow[4]={6,5,4,7};
	
	int runMin, runMax, runNum;
	cout<<"Input minimum and maximum numbers of run: ";
	cin>>runMin>>runMax;
	
	sAna="/home/kailong/ExpData/Jul2018/AnaData/ana-run-"+to_string(runMin)+"--"+to_string(runMax)+".root";
	TFile *fAna=new TFile(sAna.c_str(), "RECREATE");
	TTree *tAna=new TTree("tAna", "tree for data analysis");

	tAna->Branch("setting", &setting);
	tAna->Branch("ana", &ana, "tof[4]/D:amp[4]/D:tD[4][8][8]/D:egy[4][8]/D:xPla[4][2]/D:yPla[4][2]/D:xMCP[4][2]/D:yMCP[4][2]/D:delE[4][5]/D:tke[4]/D:beta[4]/D:gamma[4]/D:Z[4]/D:dZ[4]/D:brho[4][2]/D:AoQ[4][2]/D:Q[4][2]/D:ZmQ[4][2]/D:ZImQ[4][2]/D:A[4][2]/D:Araw[4][2]/D:Am2Q[4][2]/D:Am3Q[4][2]/D:Am2Z[4][2]/D:Am3Z[4][2]/D:dAm2Z[4][2]/D:dAm3Z[4][2]/D:Zi[4]/I:sig[4][2][2]/I:run/I");
	tAna->Branch("pid", &pid, "tof/D:xMCP[2]/D:yMCP[2]/D:delE[5]/D:tke/D:beta/D:gamma/D:Z/D:dZ/D:brho[2]/D:AoQ[2]/D:Q[2]/D:ZmQ[2]/D:ZImQ[2]/D:A[2]/D:Araw[2]/D:Am2Q[2]/D:Am3Q[2]/D:Am2Z[2]/D:Am3Z[2]/D:dAm2Z[2]/D:dAm3Z[2]/D:Zi/I:run/I");
	
	ostringstream ssRun;
	for(runNum=runMin; runNum<=runMax; runNum++)
	{
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
		
		tData->SetEstimate(-1);
		nEnt=tData->GetEntries();
		
		double *errCh=new double[nEnt];
		for(iEnt=0; iEnt<nEnt; iEnt++)
			errCh[iEnt]=0.5;
		string strCanv="gainMat_"+to_string(runNum);
		TCanvas *cGain=new TCanvas(strCanv.c_str(), strCanv.c_str());
		cGain->Divide(2,2);
		bool goodRun=true;
		TGraphErrors *gr[4];
		for(i=0; i<4; i++)
		{
			j=iLow[i];
			cGain->cd(i+1);
			string sCut="mqdcMCP.data["+to_string(i)+"]>"+to_string(QdcMcpLow[i])+"&&mqdcMCP.data["+to_string(i)+"]<3840&&mqdcMCP.data["+to_string(j)+"]>"+to_string(QdcMcpLow[j])+"&&mqdcMCP.data["+to_string(j)+"]<3840";
			
			string sDraw="mqdcMCP.data["+to_string(i)+"]-"+to_string(QdcMcpLow[i])+":mqdcMCP.data["+to_string(j)+"]-"+to_string(QdcMcpLow[j])+">>h"+to_string(i)+"_"+to_string(j);
			
			printf("Now drawing {%s} of {%s} when {%s}\n\n", sDraw.c_str(), sRoot.c_str(), sCut.c_str());
			
			long long nData=tData->Draw(sDraw.c_str(), sCut.c_str(), "goff");
			if(nData<2)
			{
				goodRun=false;
				break;
			}
			double *highGain=tData->GetV1();
			double *lowGain=tData->GetV2();
			gr[i]=new TGraphErrors(nData, lowGain, highGain, errCh, errCh);
			gr[i]->Draw("AP");
			gr[i]->SetTitle(("high"+to_string(i)+"_vs_low"+to_string(iLow[i])).c_str());
			TFitResultPtr fitRes=gr[i]->Fit("pol1","MS");
			mcpGainMat[j][0]=fitRes->Parameter(0);
			mcpGainMat[j][1]=fitRes->Parameter(1);
		}
		delete []errCh;
		if(!goodRun)
			continue;
		cGain->SaveAs(("/home/kailong/ExpData/Jul2018/Graphs/Charts/"+strCanv+".png").c_str());
		cGain->Close();
		
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
		// for(iEnt=0; iEnt<10000; iEnt++)
		{
			tData->GetEntry(iEnt);
			if(s800.trig==1||s800.trig==16)
			{
				TRandom3 r(0);
				
				memset(&ana, 0, sizeof(ana));
				ana.run=runNum;
				
				nGoodEvt=0;				
				for(iAna=0; iAna<4; iAna++)
				{
					memset(tPMT, 0, sizeof(tPMT));
					memset(timeDet, 0, sizeof(timeDet));
					memset(egyDet, 0, sizeof(egyDet));
					
					if(iAna==0) //For TAC+ADC+clock
						for(i=0; i<8; i++)
						{
							j=2*i+1;
							k=i/4;
							if(madc.data[i]>AdcTofLow&&madc.data[i]<AdcTofUp&&mqdcTOF.data[j]>QdcTofLow[i]&&mqdcTOF.data[j]<QdcUp)
							{
								ana.sig[0][k][0]++;
								ana.sig[0][k][1]=10*ana.sig[0][k][1]+(i+1);
								tPMT[i]=CALADC[i]*(madc.data[i]+r.Uniform(-0.5,0.5));
								ana.egy[0][i]=mqdcTOF.data[j]+r.Uniform(-0.5,0.5);
							}
						}
					
					if(iAna==1) //For TAC+ADC
					{
						for(i=8; i<12; i++)
						{
							j=2*(i-8)+1;
							k=2*(i-8)+9;
							m=i-8;
							p=i-4;
							if(madc.data[i]>AdcTofLow&&madc.data[i]<AdcTofUp&&mqdcTOF.data[j]>QdcTofLow[m]&&mqdcTOF.data[j]<QdcUp&&mqdcTOF.data[k]>QdcTofLow[p]&&mqdcTOF.data[k]<QdcUp)
							{
								ana.sig[1][0][0]++;
								ana.sig[1][1][0]++;
								ana.sig[1][0][1]=10*ana.sig[1][0][1]+(m+1);
								ana.sig[1][1][1]=10*ana.sig[1][1][1]+(p+1);
								ana.egy[1][m]=mqdcTOF.data[j]+r.Uniform(-0.5, 0.5);
								ana.egy[1][p]=mqdcTOF.data[k]+r.Uniform(-0.5, 0.5);
								ana.tD[1][p][m]=CALADC[i]*(madc.data[i]+r.Uniform(-0.5, 0.5));
								ana.tD[1][m][p]=-ana.tD[1][p][m];
							}
						}
						
						if(ana.sig[1][0][0]==ana.sig[1][1][0]&&ana.sig[1][0][0]>0)
						{
							coin[0]=0;
							for(i=0; i<4; i++)
								for(j=4; j<8; j++)
									if(abs(ana.tD[1][i][j])>0)
									{
										ana.tof[1]+=ana.tD[1][j][i];
										ana.amp[1]+=ana.egy[1][j]-ana.egy[1][i];
										coin[0]++;
									}
							if(coin[0]==ana.sig[1][0][0])
							{
								ana.tof[1]/=coin[0];
								ana.amp[1]/=coin[0];
							}
						}
					}
					
					if(iAna>=2)
						for(j=0; j<8; j++)
						{
							if(iAna==2)
							{
								k=j+1;
								p=j;
							}
							else
							{
								k=2*j+17;
								p=j+8;
							}
							m=2*j+1;
							n=j/4;
				
							if(mtdc.data[k]>TdcTofLow[p]&&mtdc.data[k]<TdcTofUp[p]&&mqdcTOF.data[m]>QdcTofLow[j]&&mqdcTOF.data[m]<QdcUp)
							{
								ana.sig[iAna][n][0]++;
								ana.sig[iAna][n][1]=10*ana.sig[iAna][n][1]+(j+1);
								tPMT[j]=CALTDC*(mtdc.data[k]+r.Uniform(-0.5, 0.5));
								ana.egy[iAna][j]=mqdcTOF.data[m]+r.Uniform(-0.5, 0.5);
							}
						}
					
					if(iAna!=1)
						if(ana.sig[iAna][0][0]>0&&ana.sig[iAna][1][0]>0)
						{
							for(i=0; i<8; i++)
								for(j=i+1; j<8; j++)
									if(abs(tPMT[i])>0&&abs(tPMT[j])>0)
									{
										ana.tD[iAna][i][j]=tPMT[i]-tPMT[j];
										ana.tD[iAna][j][i]=-ana.tD[iAna][i][j];
									}
									
							for(k=0; k<2; k++)
							{
								coin[k+1]=0;
								for(i=k*4; i<4+k*4; i++)
									if(abs(tPMT[i])>0)
									{
										timeDet[k]+=tPMT[i];
										egyDet[k]+=ana.egy[iAna][i];
										coin[k+1]++;
									}
								if(coin[k+1]==ana.sig[iAna][k][0])
								{
									timeDet[k]/=coin[k+1];
									egyDet[k]/=coin[k+1];
								}
							}

							ana.tof[iAna]=timeDet[1]-timeDet[0];
							ana.amp[iAna]=egyDet[1]-egyDet[0];
						}
					
					if(ana.sig[iAna][0][0]>0&&ana.sig[iAna][1][0]>0)
					{
						ana.tof[iAna]=CALTOF[iAna][0]+CALTOF[iAna][1]*ana.tof[iAna];
						for(j=0; j<2; j++)
							if(ana.sig[iAna][j][0]==4)
							{
								k=4*j;
								ana.xPla[iAna][j]=(tPMT[2+k]+tPMT[3+k]-tPMT[0+k]-tPMT[1+k])/2;
								ana.yPla[iAna][j]=(tPMT[0+k]+tPMT[3+k]-tPMT[1+k]-tPMT[2+k])/2;
							}
						b=LOF/ana.tof[iAna]/0.299792458;
						if(b>0&&b<1)
						{
							ana.beta[iAna]=b;
							ana.gamma[iAna]=1/sqrt(1-b*b);
							
							ana.tke[iAna]=0;
							nGoodPin=0;
							memset(goodPin, 0, sizeof(goodPin));
							for(q=0; q<5; q++)
								if(s800.pin[q]>AdcPinLow&&s800.pin[q]<AdcPinUp)
								{
									goodPin[q]=true;
									nGoodPin++;
									
									ana.delE[iAna][q]=CALPIN[q][0]+CALPIN[q][1]*(s800.pin[q]+r.Uniform(-0.5,0.5));
									ana.tke[iAna]+=ana.delE[iAna][q];
								}
								
							if(nGoodPin>1&&goodPin[0]) //standard condition
							// if(goodPin[0]) //loose condition
							{
								ana.tke[iAna]+=(CALPIN[5][0]+CALPIN[5][1]*(s800.pin[0]+r.Uniform(-0.5,0.5))); //consider the absorption effect of material in front of Si detectors
								
								ana.Z[iAna]=sqrt( ana.delE[iAna][0] / (log(5930/(1/b/b-1))/b/b-1) );
								ana.Z[iAna]=CALZ[0]+CALZ[1]*ana.Z[iAna];
								ana.Zi[iAna]=TMath::Nint(ana.Z[iAna]);
								ana.dZ[iAna]=ana.Z[iAna]-ana.Zi[iAna];
								
								memset(calQdcMcp, 0, sizeof(calQdcMcp));
								memset(nGoodMcp, 0, sizeof(nGoodMcp));
								for(p=0; p<4; p++)
									if(mqdcMCP.data[p]>QdcMcpLow[p]&&mqdcMCP.data[p]<QdcUp)
									{
										nGoodMcp[0]++;
										calQdcMcp[0][p]=mqdcMCP.data[p]-QdcMcpLow[p]+r.Uniform(-0.5,0.5);
									}
									
								for(p=0; p<4; p++)
								{
									if(mqdcMCP.data[p]>QdcMcpLow[p]&&mqdcMCP.data[p]<QdcUp)
									{
										nGoodMcp[1]++;
										calQdcMcp[1][p]=mqdcMCP.data[p]-QdcMcpLow[p]+r.Uniform(-0.5,0.5);
									}
									m=iLow[p];
									if(mqdcMCP.data[p]>=QdcUp&&mqdcMCP.data[m]>QdcMcpLow[m]&&mqdcMCP.data[m]<QdcUp&&mcpGainMat[m][0]+mcpGainMat[m][1]*(mqdcMCP.data[m]-QdcMcpLow[m])>QdcUp)
									{
										nGoodMcp[1]++;
										calQdcMcp[1][p]=mcpGainMat[m][0]+mcpGainMat[m][1]*(mqdcMCP.data[m]-QdcMcpLow[m]+r.Uniform(-0.5,0.5));
									}
								}
									
								for(k=0; k<2; k++)
									if(nGoodMcp[k]==4)
									{
										xMcpRaw=0;
										yMcpRaw=0;
										
										xMcpRaw=(calQdcMcp[k][0]+calQdcMcp[k][3]-calQdcMcp[k][1]-calQdcMcp[k][2])/(calQdcMcp[k][0]+calQdcMcp[k][1]+calQdcMcp[k][2]+calQdcMcp[k][3]);
										
										yMcpRaw=(calQdcMcp[k][1]+calQdcMcp[k][3]-calQdcMcp[k][0]-calQdcMcp[k][2])/(calQdcMcp[k][0]+calQdcMcp[k][1]+calQdcMcp[k][2]+calQdcMcp[k][3]);
										
										ana.xMCP[iAna][k]=CALXMCP[k][0]+CALXMCP[k][1]*xMcpRaw+CALXMCP[k][2]*pow(xMcpRaw,2)+CALXMCP[k][3]*yMcpRaw;
										
										ana.yMCP[iAna][k]=CALYMCP[k][0]+CALYMCP[k][1]*yMcpRaw+CALYMCP[k][2]*pow(yMcpRaw,2)+CALYMCP[k][3]*xMcpRaw;
										
										ana.brho[iAna][k]=BRHO0*(1+ana.xMCP[iAna][k]/DISP/100);
										ana.AoQ[iAna][k]=ana.brho[iAna][k]/ana.beta[iAna]/ana.gamma[iAna]*0.32184;
										ana.Q[iAna][k]=ana.tke[iAna]/(931.4940954*(ana.gamma[iAna]-1)*ana.AoQ[iAna][k]);
										ana.ZmQ[iAna][k]=ana.Z[iAna]-ana.Q[iAna][k];
										ana.ZImQ[iAna][k]=ana.Zi[iAna]-ana.Q[iAna][k];
										ana.A[iAna][k]=ana.AoQ[iAna][k]*ana.Q[iAna][k];
										ana.Araw[iAna][k]=ana.AoQ[iAna][k]*ana.Z[iAna];
										ana.Am2Q[iAna][k]=ana.A[iAna][k]-2*ana.Q[iAna][k];
										ana.Am3Q[iAna][k]=ana.A[iAna][k]-3*ana.Q[iAna][k];
										ana.Am2Z[iAna][k]=(ana.Q[iAna][k]-2)*ana.Zi[iAna];
										ana.Am3Z[iAna][k]=(ana.Q[iAna][k]-3)*ana.Zi[iAna];
										ana.dAm2Z[iAna][k]=ana.Am2Z[iAna][k]-TMath::Nint(ana.Am2Z[iAna][k]);
										ana.dAm3Z[iAna][k]=ana.Am3Z[iAna][k]-TMath::Nint(ana.Am3Z[iAna][k]);
										if(ana.Am3Z[iAna][k]<0)
											ana.dAm3Z[iAna][k]+=1;
									}
								
								if(nGoodMcp[0]==4||nGoodMcp[1]==4) //standard condition
								// if(nGoodMcp[0]>=0) //loose condition
									nGoodEvt++;
							}
						}
					}
				}
				
				if(nGoodEvt>1) //standard condition
				// if(nGoodEvt>0) //loose condition
				{
					//begin to fill branch of pid
					goodPid=false;
					memset(&pid, 0, sizeof(pid));
					pid.run=runNum;

					if(0.0625*(s800.mesyTDC[2]-s800.mesyTDC[0])>=60&&0.0625*(s800.mesyTDC[2]-s800.mesyTDC[0])<=150)
					{
						pid.tof=CALTOF_PID[0]+CALTOF_PID[1]*62.5*((s800.mesyTDC[2]+r.Uniform(-0.5,0.5))-(s800.mesyTDC[0]+r.Uniform(-0.5,0.5)));
						
						b=LOF_PID/pid.tof/0.299792458;
						if(b>0&&b<1)
						{
							pid.beta=b;
							pid.gamma=1/sqrt(1-b*b);
							nGoodPin=0;
							memset(goodPin, 0, sizeof(goodPin));
							pid.tke=0;
							for(i=0; i<5; i++)
								if(s800.pin[i]>AdcPinLow&&s800.pin[i]<AdcPinUp)
								{
									goodPin[i]=true;
									nGoodPin++;
									pid.delE[i]=CALPIN_PID[i][0]+CALPIN_PID[i][1]*(s800.pin[0]+r.Uniform(-0.5,0.5));
									pid.tke+=pid.delE[i];
								}
							if(nGoodPin>1&&goodPin[0]) //standard condition
							// if(goodPin[0]) //loose condition
							{
								pid.tke+=CALPIN_PID[5][0]+CALPIN_PID[5][1]*(s800.pin[0]+r.Uniform(-0.5,0.5));
								pid.Z=sqrt( s800.pin[0] / (log(5930/(1/b/b-1))/b/b-1) );
								pid.Z=CALZ_PID[0]+CALZ_PID[1]*pid.Z;
								pid.Zi=TMath::Nint(pid.Z);
								pid.dZ=pid.Z-pid.Zi;
								
								memset(calQdcMcp, 0, sizeof(calQdcMcp));
								memset(nGoodMcp, 0, sizeof(nGoodMcp));
								for(i=0; i<4; i++)
									if(mqdcMCP.data[i]>QdcMcpLow_Pid[i]&&mqdcMCP.data[i]<QdcUp)
									{
										nGoodMcp[0]++;
										calQdcMcp[0][i]=mqdcMCP.data[i]-QdcMcpLow_Pid[i]+r.Uniform(-0.5,0.5);
									}
									
								for(p=0; p<4; p++)
								{
									if(mqdcMCP.data[p]>QdcMcpLow_Pid[p]&&mqdcMCP.data[p]<QdcUp)
									{
										nGoodMcp[1]++;
										calQdcMcp[1][p]=mqdcMCP.data[p]-QdcMcpLow_Pid[p]+r.Uniform(-0.5,0.5);
									}
									m=iLow[p];
									if(mqdcMCP.data[p]>=QdcUp&&mqdcMCP.data[m]>QdcMcpLow_Pid[m]&&mqdcMCP.data[m]<QdcUp&&mcpGainMat[m][0]+mcpGainMat[m][1]*(mqdcMCP.data[m]-QdcMcpLow[m])>QdcUp)
									{
										nGoodMcp[1]++;
										calQdcMcp[1][p]=mcpGainMat[m][0]+mcpGainMat[m][1]*(mqdcMCP.data[m]-QdcMcpLow_Pid[m]+r.Uniform(-0.5,0.5));
									}
								}
									
								for(i=0; i<2; i++)
									if(nGoodMcp[i]==4)
									{
										xMcpRaw=0;
										yMcpRaw=0;
										
										xMcpRaw=(calQdcMcp[i][1]+calQdcMcp[i][2]-calQdcMcp[i][0]-calQdcMcp[i][3])/(calQdcMcp[i][0]+calQdcMcp[i][1]+calQdcMcp[i][2]+calQdcMcp[i][3]);
										
										yMcpRaw=(calQdcMcp[i][1]+calQdcMcp[i][3]-calQdcMcp[i][0]-calQdcMcp[i][2])/(calQdcMcp[i][0]+calQdcMcp[i][1]+calQdcMcp[i][2]+calQdcMcp[i][3]);
										
										pid.xMCP[i]=CALXMCP_PID[i][0]+CALXMCP_PID[i][1]*xMcpRaw+CALXMCP_PID[i][2]*pow(xMcpRaw,2)+CALXMCP_PID[i][3]*pow(xMcpRaw,3);
										
										pid.yMCP[i]=CALYMCP_PID[i][0]+CALYMCP_PID[i][1]*yMcpRaw+CALYMCP_PID[i][2]*pow(yMcpRaw,2)+CALYMCP_PID[i][3]*pow(yMcpRaw,3);
										
										pid.brho[i]=BRHO0_PID*(1+pid.xMCP[i]/DISP_PID/100);
										pid.AoQ[i]=pid.brho[i]/pid.beta/pid.gamma*0.32184;
										pid.Q[i]=pid.tke/(931.4940954*(pid.gamma-1)*pid.AoQ[i]);
										pid.ZmQ[i]=pid.Z-pid.Q[i];
										pid.ZImQ[i]=pid.Zi-pid.Q[i];
										pid.A[i]=pid.AoQ[i]*pid.Q[i];
										pid.Araw[i]=pid.AoQ[i]*pid.Z;
										pid.Am2Q[i]=pid.A[i]-2*pid.Q[i];
										pid.Am3Q[i]=pid.A[i]-3*pid.Q[i];
										pid.Am2Z[i]=(pid.Q[i]-2)*pid.Zi;
										pid.Am3Z[i]=(pid.Q[i]-3)*pid.Zi;
										pid.dAm2Z[i]=pid.Am2Z[i]-TMath::Nint(pid.Am2Z[i]);
										pid.dAm3Z[i]=pid.Am3Z[i]-TMath::Nint(pid.Am3Z[i]);
										if(pid.Am3Z[i]<0)
											pid.dAm3Z[i]+=1;
									}
								if(nGoodMcp[0]==4||nGoodMcp[1]==4) //standard condition
								// if(nGoodMcp[0]>=0) //loose condition
									goodPid=true;
							}
						}
					}
					// if(goodPid)
					tAna->Fill();
				}	
			}	
		}//end of whole tree
		fRoot->Close();
	}
	fAna->cd();
	tAna->Write();
	fAna->Close();	
}//end of whole function

#ifndef __CINT__
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
#endif
