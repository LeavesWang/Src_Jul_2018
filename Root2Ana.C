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
	int sig[2][2];
	double tof;
	double amp;
	double tD[8][8];
	double egy[8];
	double xMCP[2];
	double yMCP[2];
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
};

void Root2Ana()
{		
	const int LADC=100, HADC=7680, LTDC=1, HTDC=65536, LQDCTOF=800, HQDC=3840;
	const int LQDCMCP[8]={726, 730, 745, 742, 700,700,700,700};
	const double CALADC[12]={6.46209, 6.59645, 6.56230, 6.57185, 6.44156, 6.58265, 6.64827, 6.52219, 6.45537, 6.42844, 6.65406, 6.43436};  //unit: ps/ch
	const double CALTDC=3.90625; //ps/ch
	const double CALPIN[6][2]={{0,0.6951}, {0,0.6558}, {0,2.9832}, {0,2.7269}, {0,2.9703}, {0,0.4881}}; //Mev/ch  //0.6951 is the original slope and 0.4886 is related to the material in front of Si detectors.

	const double CALXMCP[2][4]={{-4.11869, -26.6253, -3.38656, -19.399}, {0,1,0,0}}; //mm, mm/ch, mm/ch^2, mm/ch^3
	const double CALYMCP[2][4]={{0,1,0,0}, {0,1,0,0}};
	const double CALTOF[2]={570.456, -0.001}; //ns, ns/ps

	const double BRHO0=3.7221; //Tm
	const double DISP=112; //m/100
	const double LOF=60.74; //m

	const double CALZ[3]={1.191, 5.9377, 0};
	
	string sRoot, sAna;
	StrtMesytec madc, mtdc, mqdcTOF, mqdcMCP;
	StrtS800 s800;
	StrtAna anaADC[2], anaTDC[2];

	double tPMT[8];
	double tAdcDet[2], tTdcDet[2], egyDet[2];
	
	int iEntry;
	int i, j, k, m, n, p, q, u;
	int coin, coin1, coin2;
	int goodSi, goodMCP[2];
	double b;

	int runMin, runMax, runNum;
	cout<<"Input minimum and maximum numbers of run: ";
	cin>>runMin>>runMax;
	ostringstream ssRun;	
	for(runNum=runMin; runNum<=runMax; runNum++)
	{
		ssRun.str("");
		ssRun<<setw(4)<<setfill('0')<<runNum;
		
		sRoot="/home/kailong/ExpData/Jul2018/RootData/run-"+ssRun.str()+"-00.root";
		if(access("/home/kailong/ExpData/Jul2018/AnaData", F_OK)!=0)
            system("mkdir /home/kailong/ExpData/Jul2018/AnaData");
		sAna="/home/kailong/ExpData/Jul2018/AnaData/ana-run-"+ssRun.str()+"-00.root";
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
		
		TFile *fAna=new TFile(sAna.c_str(), "RECREATE");
		TTree *tAnaADC[2], *tAnaTDC[2];
		tAnaADC[0]=new TTree("tAnaADC1", "tree for clk+TAC+ADC analysis");
		tAnaADC[1]=new TTree("tAnaADC2", "tree for TAC+ADC analysis");
		tAnaTDC[0]=new TTree("tAnaTDC1", "tree for CFD+TDC analysis");
		tAnaTDC[1]=new TTree("tAnaTDC2", "tree for MCFD16+TDC analysis");
		tAnaADC[0]->Branch("anaADC1", &anaADC[0], "sig[2][2]/I:tof/D:amp/D:tD[8][8]/D:egy[8]/D:xMCP[2]/D:yMCP[2]/D:delE[5]/D:tke/D:beta/D:gamma/D:Z/D:dZ/D:brho[2]/D:AoQ[2]/D:Q[2]/D:ZmQ[2]/D:ZImQ[2]/D:A[2]/D:Araw[2]/D:Am2Q[2]/D:Am3Q[2]/D:Am2Z[2]/D:Am3Z[2]/D:dAm2Z[2]/D:dAm3Z[2]/D:Zi/I");
		tAnaADC[1]->Branch("anaADC2", &anaADC[1], "sig[2][2]/I:tof/D:amp/D:tD[8][8]/D:egy[8]/D:xMCP[2]/D:yMCP[2]/D:delE[5]/D:tke/D:beta/D:gamma/D:Z/D:dZ/D:brho[2]/D:AoQ[2]/D:Q[2]/D:ZmQ[2]/D:ZImQ[2]/D:A[2]/D:Araw[2]/D:Am2Q[2]/D:Am3Q[2]/D:Am2Z[2]/D:Am3Z[2]/D:dAm2Z[2]/D:dAm3Z[2]/D:Zi/I");
		tAnaTDC[0]->Branch("anaTDC1", &anaTDC[0], "sig[2][2]/I:tof/D:amp/D:tD[8][8]/D:egy[8]/D:xMCP[2]/D:yMCP[2]/D:delE[5]/D:tke/D:beta/D:gamma/D:Z/D:dZ/D:brho[2]/D:AoQ[2]/D:Q[2]/D:ZmQ[2]/D:ZImQ[2]/D:A[2]/D:Araw[2]/D:Am2Q[2]/D:Am3Q[2]/D:Am2Z[2]/D:Am3Z[2]/D:dAm2Z[2]/D:dAm3Z[2]/D:Zi/I");
		tAnaTDC[1]->Branch("anaTDC2", &anaTDC[1], "sig[2][2]/I:tof/D:amp/D:tD[8][8]/D:egy[8]/D:xMCP[2]/D:yMCP[2]/D:delE[5]/D:tke/D:beta/D:gamma/D:Z/D:dZ/D:brho[2]/D:AoQ[2]/D:Q[2]/D:ZmQ[2]/D:ZImQ[2]/D:A[2]/D:Araw[2]/D:Am2Q[2]/D:Am3Q[2]/D:Am2Z[2]/D:Am3Z[2]/D:dAm2Z[2]/D:dAm3Z[2]/D:Zi/I");
		
		for(iEntry=0; iEntry<tData->GetEntries(); iEntry++)
		// for(iEntry=0; iEntry<100; iEntry++)
		{
			tData->GetEntry(iEntry);
			if(s800.trig==1)
			{
				TRandom3 r(0);
			
				memset(anaADC, 0, sizeof(anaADC));
				memset(anaTDC, 0, sizeof(anaTDC));
				
				//Fill tree of anaADC[0]
				memset(tPMT, 0, sizeof(tPMT));
				for(i=0; i<8; i++)
				{
					j=2*i+1;
					k=i/4;
					if(mtdc.data[i]>LTDC&&mtdc.data[i]<HTDC&&mqdcTOF.data[j]>LQDCTOF&&mqdcTOF.data[j]<HQDC)
					{
						anaADC[0].sig[k][0]++;
						anaADC[0].sig[k][1]=10*anaADC[0].sig[k][1]+(i+1);
						tPMT[i]=CALADC[i]*(madc.data[i]+r.Uniform(-0.5,0.5));
						anaADC[0].egy[i]=mqdcTOF.data[j]+r.Uniform(-0.5,0.5);
					}
				}

				if(anaADC[0].sig[0][0]>0&&anaADC[0].sig[1][0]>0)
				{
					for(i=0; i<8; i++)
						for(j=i+1; j<8; j++)
							if(abs(tPMT[i])>0&&abs(tPMT[j])>0)
							{
								anaADC[0].tD[i][j]=tPMT[i]-tPMT[j];
								anaADC[0].tD[j][i]=-anaADC[0].tD[i][j];
							}
					memset(tAdcDet, 0, sizeof(tAdcDet));
					memset(egyDet, 0, sizeof(egyDet));
					
					coin1=0;
					for(i=0; i<4; i++)
						if(abs(tPMT[i])>0)
						{
							tAdcDet[0]+=tPMT[j];
							egyDet[0]+=anaADC[0].egy[j];
							coin1++;
						}
					if(coin1==anaADC[0].sig[0][0])
					{
						tAdcDet[0]/=coin1;
						egyDet[0]/=coin1;
					}
					coin2=0;
					for(i=4; i<8; i++)
						if(abs(tPMT[i])>0)
						{
							tAdcDet[1]+=tPMT[i];
							egyDet[1]+=anaADC[0].egy[j];
							coin2++;
						}
					if(coin2==anaADC[0].sig[1][0])
					{
						tAdcDet[1]/=coin2;
						egyDet[1]/=coin2;
					}
					if(coin1==anaADC[0].sig[0][0]&&coin2==anaADC[0].sig[1][0])
					{
						anaADC[0].tof=tAdcDet[1]-tAdcDet[0];
						anaADC[0].tof=CALTOF[0]+CALTOF[1]*anaADC[0].tof;
						anaADC[0].amp=egyDet[1]-egyDet[0];
						
						memset(goodMCP, 0, sizeof(goodMCP));
						for(p=0; p<8; p++)
							if(mqdcMCP.data[p]>LQDCMCP[p]&&mqdcMCP.data[p]<HQDC)
							{
								m=p/4;
								goodMCP[m]++;
							}
							
						for(k=0; k<2; k++)
							if(goodMCP[k]==4)
							{
								anaADC[0].xMCP[k]=1.0*(mqdcMCP.data[0+4*k]-LQDCMCP[0+4*k]+mqdcMCP.data[3+4*k]-LQDCMCP[3+4*k]-mqdcMCP.data[1+4*k]-LQDCMCP[1+4*k]-mqdcMCP.data[2+4*k]-LQDCMCP[2+4*k])/1.0/(mqdcMCP.data[0+4*k]-LQDCMCP[0+4*k]+mqdcMCP.data[3+4*k]-LQDCMCP[3+4*k]+mqdcMCP.data[1+4*k]-LQDCMCP[1+4*k]+mqdcMCP.data[2+4*k]-LQDCMCP[2+4*k]);
								anaADC[0].xMCP[k]=CALXMCP[k][0]+CALXMCP[k][1]*anaADC[0].xMCP[k]+CALXMCP[k][2]*pow(anaADC[0].xMCP[k],2)+CALXMCP[k][3]*pow(anaADC[0].xMCP[k],3);
								
								anaADC[0].yMCP[k]=1.0*(mqdcMCP.data[1+4*k]-LQDCMCP[1+4*k]+mqdcMCP.data[3+4*k]-LQDCMCP[3+4*k]-mqdcMCP.data[0+4*k]-LQDCMCP[0+4*k]-mqdcMCP.data[2+4*k]-LQDCMCP[2+4*k])/1.0/(mqdcMCP.data[0+4*k]-LQDCMCP[0+4*k]+mqdcMCP.data[3+4*k]-LQDCMCP[3+4*k]+mqdcMCP.data[1+4*k]-LQDCMCP[1+4*k]+mqdcMCP.data[2+4*k]-LQDCMCP[2+4*k]);
								anaADC[0].yMCP[k]=CALYMCP[k][0]+CALYMCP[k][1]*anaADC[0].yMCP[k]+CALYMCP[k][2]*pow(anaADC[0].yMCP[k],2)+CALYMCP[k][3]*pow(anaADC[0].yMCP[k],3);
							}
						
						anaADC[0].tke=0;
						goodSi=0;
						for(q=0; q<5; q++)
							if(s800.pin[q]>100&&s800.pin[q]<4000)
							{
								anaADC[0].delE[q]=CALPIN[q][0]+CALPIN[q][1]*s800.pin[q];
								anaADC[0].tke+=anaADC[0].delE[q];
								goodSi++;
							}
						if(s800.pin[0]>100&&s800.pin[0]<4000)
							anaADC[0].tke+=(CALPIN[5][0]+CALPIN[5][1]*s800.pin[0]); //consider the absorption effect of material in front of Si detectors
							
						b=LOF/anaADC[0].tof/0.299792458;
						if(b>0&&b<1)
						{
							anaADC[0].beta=b;
							anaADC[0].gamma=1/sqrt(1-b*b);
							if(s800.pin[0]>100&&s800.pin[0]<4000&&goodSi>1)
							{
								anaADC[0].Z=anaADC[0].delE[0]/sqrt(1/b/b*log(5930.0/(1/b/b-1))-1);
								anaADC[0].Z=CALZ[0]+CALZ[1]*anaADC[0].Z+CALZ[2]*pow(anaADC[0].Z,2);
								anaADC[0].Zi=TMath::Nint(anaADC[0].Z);
								anaADC[0].dZ=anaADC[0].Z-anaADC[0].Zi;
								for(k=0; k<2; k++)
									if(goodMCP[k]==4)
									{
										anaADC[0].brho[k]=BRHO0*(1+anaADC[0].xMCP[k]/DISP/100);
										anaADC[0].AoQ[k]=anaADC[0].brho[k]/anaADC[0].beta/anaADC[0].gamma*0.32184;
										anaADC[0].Q[k]=anaADC[0].tke/(931.4940954*(anaADC[0].gamma-1)*anaADC[0].AoQ[k]);
										anaADC[0].ZmQ[k]=anaADC[0].Z-anaADC[0].Q[k];
										anaADC[0].ZImQ[k]=anaADC[0].Zi-anaADC[0].Q[k];
										anaADC[0].A[k]=anaADC[0].AoQ[k]*anaADC[0].Q[k];
										anaADC[0].Araw[k]=anaADC[0].AoQ[k]*anaADC[0].Z;
										anaADC[0].Am2Q[k]=anaADC[0].A[k]-2*anaADC[0].Q[k];
										anaADC[0].Am3Q[k]=anaADC[0].A[k]-3*anaADC[0].Q[k];
										anaADC[0].Am2Z[k]=(anaADC[0].Q[k]-2)*anaADC[0].Zi;
										anaADC[0].Am3Z[k]=(anaADC[0].Q[k]-3)*anaADC[0].Zi;
										anaADC[0].dAm2Z[k]=anaADC[0].Am2Z[k]-TMath::Nint(anaADC[0].Am2Z[k]);
										anaADC[0].dAm3Z[k]=anaADC[0].Am3Z[k]-TMath::Nint(anaADC[0].Am3Z[k]);
										if(anaADC[0].Am3Z[k]<0)
											anaADC[0].dAm3Z[k]+=1;
									}
							}
						}
						if(b>0&&b<1&&s800.pin[0]>100&&s800.pin[0]<4000&&goodSi>1&&goodMCP[0]==4)
						{
							// cout<<"anaADC[0].Z: "<<anaADC[0].Z<<endl;
							// cout<<"anaADC[0].brho[0]: "<<anaADC[0].brho[0]<<endl;
							// cout<<"anaADC[0].AoQ[0]: "<<anaADC[0].AoQ[k]<<endl;
							tAnaADC[0]->Fill();
						}
					}
				}				
				
				//Fill tree of tAnaADC[1]
				for(i=8; i<12; i++)
				{
					j=2*(i-8)+1;
					k=2*(i-8)+9;
					m=i-8;
					p=i-4;
					if(madc.data[i]>LADC&&madc.data[i]<HADC&&mqdcTOF.data[j]>LQDCTOF&&mqdcTOF.data[j]<HQDC&&mqdcTOF.data[k]>LQDCTOF&&mqdcTOF.data[k]<HQDC)
					{					
						anaADC[1].sig[0][0]++;
						anaADC[1].sig[1][0]++;
						anaADC[1].sig[0][1]=10*anaADC[1].sig[0][1]+(m+1);
						anaADC[1].sig[1][1]=10*anaADC[1].sig[1][1]+(p+1);
						anaADC[1].egy[m]=mqdcTOF.data[j]+r.Uniform(-0.5, 0.5);
						anaADC[1].egy[p]=mqdcTOF.data[k]+r.Uniform(-0.5, 0.5);
						anaADC[1].tD[p][m]=CALADC[i]*(madc.data[i]+r.Uniform(-0.5, 0.5));
						anaADC[1].tD[m][p]=-anaADC[1].tD[p][m];
					}
				}
				
				if(anaADC[1].sig[0][0]==anaADC[1].sig[1][0]&&anaADC[1].sig[0][0]>0)
				{
					coin=0;
					for(i=0; i<4; i++)
						for(j=4; j<8; j++)
							if(abs(anaADC[1].tD[i][j])>0)
							{
								anaADC[1].tof+=anaADC[1].tD[j][i];
								anaADC[1].amp+=anaADC[1].egy[j]-anaADC[1].egy[i];
								coin++;
							}
					if(coin==anaADC[1].sig[0][0])
					{
						anaADC[1].tof/=coin;
						anaADC[1].tof=CALTOF[0]+CALTOF[1]*anaADC[1].tof;
						anaADC[1].amp/=coin;
						
						memset(goodMCP, 0, sizeof(goodMCP));
						for(p=0; p<8; p++)
							if(mqdcMCP.data[p]>LQDCMCP[p]&&mqdcMCP.data[p]<HQDC)
							{
								m=p/4;
								goodMCP[m]++;
							}
						for(k=0; k<2; k++)
							if(goodMCP[k]==4)
							{
								anaADC[1].xMCP[k]=1.0*(mqdcMCP.data[0+4*k]-LQDCMCP[0+4*k]+mqdcMCP.data[3+4*k]-LQDCMCP[3+4*k]-mqdcMCP.data[1+4*k]-LQDCMCP[1+4*k]-mqdcMCP.data[2+4*k]-LQDCMCP[2+4*k])/1.0/(mqdcMCP.data[0+4*k]-LQDCMCP[0+4*k]+mqdcMCP.data[3+4*k]-LQDCMCP[3+4*k]+mqdcMCP.data[1+4*k]-LQDCMCP[1+4*k]+mqdcMCP.data[2+4*k]-LQDCMCP[2+4*k]);
								anaADC[1].xMCP[k]=CALXMCP[k][0]+CALXMCP[k][1]*anaADC[1].xMCP[k]+CALXMCP[k][2]*pow(anaADC[1].xMCP[k],2)+CALXMCP[k][3]*pow(anaADC[1].xMCP[k],3);
								
								anaADC[1].yMCP[k]=1.0*(mqdcMCP.data[1+4*k]-LQDCMCP[1+4*k]+mqdcMCP.data[3+4*k]-LQDCMCP[3+4*k]-mqdcMCP.data[0+4*k]-LQDCMCP[0+4*k]-mqdcMCP.data[2+4*k]-LQDCMCP[2+4*k])/1.0/(mqdcMCP.data[0+4*k]-LQDCMCP[0+4*k]+mqdcMCP.data[3+4*k]-LQDCMCP[3+4*k]+mqdcMCP.data[1+4*k]-LQDCMCP[1+4*k]+mqdcMCP.data[2+4*k]-LQDCMCP[2+4*k]);
								anaADC[1].yMCP[k]=CALYMCP[k][0]+CALYMCP[k][1]*anaADC[1].yMCP[k]+CALYMCP[k][2]*pow(anaADC[1].yMCP[k],2)+CALYMCP[k][3]*pow(anaADC[1].yMCP[k],3);
							}
						anaADC[1].tke=0;
						goodSi=0;
						for(q=0; q<5; q++)
							if(s800.pin[q]>100&&s800.pin[q]<4000)
							{
								anaADC[1].delE[q]=CALPIN[q][0]+CALPIN[q][1]*s800.pin[q];
								anaADC[1].tke+=anaADC[0].delE[q];
								goodSi++;
							}
						if(s800.pin[0]>100&&s800.pin[0]<4000)
							anaADC[1].tke+=(CALPIN[5][0]+CALPIN[5][1]*s800.pin[0]); //consider the absorption effect of material in front of Si detectors
							
						b=LOF/anaADC[1].tof/0.299792458;
						if(b>0&&b<1)
						{
							anaADC[1].beta=b;
							anaADC[1].gamma=1/sqrt(1-b*b);
							if(s800.pin[0]>100&&s800.pin[0]<4000&&goodSi>1)
							{
								anaADC[1].Z=anaADC[1].delE[0]/sqrt(1/b/b*log(5930.0/(1/b/b-1))-1);
								anaADC[1].Z=CALZ[0]+CALZ[1]*anaADC[1].Z+CALZ[2]*pow(anaADC[1].Z,2);
								anaADC[1].Zi=TMath::Nint(anaADC[1].Z);
								anaADC[1].dZ=anaADC[1].Z-anaADC[1].Zi;
								for(k=0; k<2; k++)
									if(goodMCP[k]==4)
									{
										anaADC[1].brho[k]=BRHO0*(1+(anaADC[1].xMCP[k]/DISP)/100);
										anaADC[1].AoQ[k]=anaADC[1].brho[k]/anaADC[1].beta/anaADC[1].gamma*0.32184;
										anaADC[1].Q[k]=anaADC[1].tke/(931.4940954*(anaADC[1].gamma-1)*anaADC[1].AoQ[k]);
										anaADC[1].ZmQ[k]=anaADC[1].Z-anaADC[1].Q[k];
										anaADC[1].ZImQ[k]=anaADC[1].Zi-anaADC[1].Q[k];
										anaADC[1].A[k]=anaADC[1].AoQ[k]*anaADC[1].Q[k];
										anaADC[1].Araw[k]=anaADC[1].AoQ[k]*anaADC[1].Z;
										anaADC[1].Am2Q[k]=anaADC[1].A[k]-2*anaADC[1].Q[k];
										anaADC[1].Am3Q[k]=anaADC[1].A[k]-3*anaADC[1].Q[k];
										anaADC[1].Am2Z[k]=(anaADC[1].Q[k]-2)*anaADC[1].Zi;
										anaADC[1].Am3Z[k]=(anaADC[1].Q[k]-3)*anaADC[1].Zi;
										anaADC[1].dAm2Z[k]=anaADC[1].Am2Z[k]-TMath::Nint(anaADC[1].Am2Z[k]);
										anaADC[1].dAm3Z[k]=anaADC[1].Am3Z[k]-TMath::Nint(anaADC[1].Am3Z[k]);
										if(anaADC[1].Am3Z[k]<0)
											anaADC[1].dAm3Z[k]+=1;
									}
							}
						}
						if(b>0&&b<1&&s800.pin[0]>100&&s800.pin[0]<4000&&goodSi>1&&goodMCP[0]==4)
							tAnaADC[1]->Fill();
					}	
				}
				
				//Fill tree of anaTDC[0] and anaTDC[1]
				for(i=0; i<2; i++)
				{				
					memset(tPMT, 0, sizeof(tPMT));
					for(j=0; j<8; j++)
					{ 
						if(i==0)
							k=j+1;
						else
							k=2*j+17;
						m=2*j+1;
						n=j/4;
			
						if(mtdc.data[k]>LTDC&&mtdc.data[k]<HTDC&&mqdcTOF.data[m]>LQDCTOF&&mqdcTOF.data[m]<HQDC)
						{
							anaTDC[i].sig[n][0]++;
							anaTDC[i].sig[n][1]=10*anaTDC[i].sig[n][1]+(j+1);
							tPMT[j]=CALTDC*(mtdc.data[k]+r.Uniform(-0.5, 0.5));
							anaTDC[i].egy[j]=mqdcTOF.data[m]+r.Uniform(-0.5, 0.5);
						}
					}

					if(anaTDC[i].sig[0][0]>0&&anaTDC[i].sig[1][0]>0)
					{
						for(j=0; j<8; j++)
							for(k=j+1; k<8; k++)
								if(abs(tPMT[j])>0&&abs(tPMT[k])>0)
								{
									anaTDC[i].tD[j][k]=tPMT[j]-tPMT[k];
									anaTDC[i].tD[k][j]=-anaTDC[i].tD[j][k];
								}
								
						memset(tTdcDet, 0, sizeof(tTdcDet));
						memset(egyDet, 0, sizeof(egyDet));
						
						coin1=0;
						for(j=0; j<4; j++)
							if(abs(tPMT[j])>0)
							{
								tTdcDet[0]+=tPMT[j];
								egyDet[0]+=anaTDC[i].egy[j];
								coin1++;
							}
						if(coin1==anaTDC[i].sig[0][0])
						{
							tTdcDet[0]/=coin1;
							egyDet[0]/=coin1;
						}
						
						coin2=0;
						for(j=4; j<8; j++)
							if(abs(tPMT[j])>0)
							{
								tTdcDet[1]+=tPMT[j];
								egyDet[1]+=anaTDC[i].egy[j];
								coin2++;
							}
						if(coin2==anaTDC[i].sig[1][0])
						{
							tTdcDet[1]/=coin2;
							egyDet[1]/=coin2;
						}
						
						if(coin1==anaTDC[i].sig[0][0]&&coin2==anaTDC[i].sig[1][0])
						{
							anaTDC[i].tof=tTdcDet[1]-tTdcDet[0];
							anaTDC[i].tof=CALTOF[0]+CALTOF[1]*anaTDC[i].tof;
							anaTDC[i].amp=egyDet[1]-egyDet[0];
							memset(goodMCP, 0, sizeof(goodMCP));
							for(p=0; p<8; p++)
								if(mqdcMCP.data[p]>LQDCMCP[p]&&mqdcMCP.data[p]<HQDC)
								{
									m=p/4;
									goodMCP[m]++;
								}
							for(k=0; k<2; k++)
								if(goodMCP[k]==4)
								{
									anaTDC[i].xMCP[k]=1.0*(mqdcMCP.data[0+4*k]-LQDCMCP[0+4*k]+mqdcMCP.data[3+4*k]-LQDCMCP[3+4*k]-mqdcMCP.data[1+4*k]-LQDCMCP[1+4*k]-mqdcMCP.data[2+4*k]-LQDCMCP[2+4*k])/1.0/(mqdcMCP.data[0+4*k]-LQDCMCP[0+4*k]+mqdcMCP.data[3+4*k]-LQDCMCP[3+4*k]+mqdcMCP.data[1+4*k]-LQDCMCP[1+4*k]+mqdcMCP.data[2+4*k]-LQDCMCP[2+4*k]);
									anaTDC[i].xMCP[k]=CALXMCP[k][0]+CALXMCP[k][1]*anaTDC[i].xMCP[k]+CALXMCP[k][2]*pow(anaTDC[i].xMCP[k],2)+CALXMCP[k][3]*pow(anaTDC[i].xMCP[k],3);
									
									anaTDC[i].yMCP[k]=1.0*(mqdcMCP.data[1+4*k]-LQDCMCP[1+4*k]+mqdcMCP.data[3+4*k]-LQDCMCP[3+4*k]-mqdcMCP.data[0+4*k]-LQDCMCP[0+4*k]-mqdcMCP.data[2+4*k]-LQDCMCP[2+4*k])/1.0/(mqdcMCP.data[0+4*k]-LQDCMCP[0+4*k]+mqdcMCP.data[3+4*k]-LQDCMCP[3+4*k]+mqdcMCP.data[1+4*k]-LQDCMCP[1+4*k]+mqdcMCP.data[2+4*k]-LQDCMCP[2+4*k]);
									anaTDC[i].yMCP[k]=CALYMCP[k][0]+CALYMCP[k][1]*anaTDC[i].yMCP[k]+CALYMCP[k][2]*pow(anaTDC[i].yMCP[k],2)+CALYMCP[k][3]*pow(anaTDC[i].yMCP[k],3);
								}
							anaTDC[i].tke=0;
							goodSi=0;
							for(q=0; q<5; q++)
								if(s800.pin[q]>100&&s800.pin[q]<4000)
								{
									anaTDC[i].delE[q]=CALPIN[q][0]+CALPIN[q][1]*s800.pin[q];
									anaTDC[i].tke+=anaADC[0].delE[q];
									goodSi++;
								}
							if(s800.pin[0]>100&&s800.pin[0]<4000)
								anaTDC[i].tke+=(CALPIN[5][0]+CALPIN[5][1]*s800.pin[0]); //consider the absorption effect of material in front of Si detectors
								
							b=LOF/anaTDC[i].tof/0.299792458;
							if(b>0&&b<1)
							{
								anaTDC[i].beta=b;
								anaTDC[i].gamma=1/sqrt(1-b*b);
								if(s800.pin[0]>100&&s800.pin[0]<4000&&goodSi>1)
								{
									anaTDC[i].Z=anaTDC[i].delE[0]/sqrt(1/b/b*log(5930.0/(1/b/b-1))-1);
									anaTDC[i].Z=CALZ[0]+CALZ[1]*anaTDC[i].Z+CALZ[2]*pow(anaTDC[i].Z,2);
									anaTDC[i].Zi=TMath::Nint(anaTDC[i].Z);
									anaTDC[i].dZ=anaTDC[i].Z-anaTDC[i].Zi;
									for(k=0; k<2; k++)
										if(goodMCP[k]==4)
										{
											anaTDC[i].brho[k]=BRHO0*(1+(anaTDC[i].xMCP[k]/DISP)/100);
											anaTDC[i].AoQ[k]=anaTDC[i].brho[k]/anaTDC[i].beta/anaTDC[i].gamma*0.32184;
											anaTDC[i].Q[k]=anaTDC[i].tke/(931.4940954*(anaTDC[i].gamma-1)*anaTDC[i].AoQ[k]);
											anaTDC[i].ZmQ[k]=anaTDC[i].Z-anaTDC[i].Q[k];
											anaTDC[i].ZImQ[k]=anaTDC[i].Zi-anaTDC[i].Q[k];
											anaTDC[i].A[k]=anaTDC[i].AoQ[k]*anaTDC[i].Q[k];
											anaTDC[i].Araw[k]=anaTDC[i].AoQ[k]*anaTDC[i].Z;
											anaTDC[i].Am2Q[k]=anaTDC[i].A[k]-2*anaTDC[i].Q[k];
											anaTDC[i].Am3Q[k]=anaTDC[i].A[k]-3*anaTDC[i].Q[k];
											anaTDC[i].Am2Z[k]=(anaTDC[i].Q[k]-2)*anaTDC[i].Zi;
											anaTDC[i].Am3Z[k]=(anaTDC[i].Q[k]-3)*anaTDC[i].Zi;
											anaTDC[i].dAm2Z[k]=anaTDC[i].Am2Z[k]-TMath::Nint(anaTDC[i].Am2Z[k]);
											anaTDC[i].dAm3Z[k]=anaTDC[i].Am3Z[k]-TMath::Nint(anaTDC[i].Am3Z[k]);
											if(anaTDC[i].Am3Z[k]<0)
												anaTDC[i].dAm3Z[k]+=1;
										}
								}
							}
							if(b>0&&b<1&&s800.pin[0]>100&&s800.pin[0]<4000&&goodSi>1&&goodMCP[0]==4)
								tAnaTDC[i]->Fill();
						}
					}
				}
			}	
		}//end of whole tree
		fAna->Write();
		fAna->Close();
		fRoot->Close();
	}//end of runs
}//end of whole function

#ifndef __CINT__
/*
void StandaloneApplication(int argc, char** argv)
{
	Root2Ana();
}

int main(int argc, char** argv)
{
	TApplication app("ROOT Application", &argc, argv);
	StandaloneApplication(app.Argc(), app.Argv());
	app.Run();
	return 0;
}
*/
int main()
{
	Root2Ana();
	return 0;
}
#endif
