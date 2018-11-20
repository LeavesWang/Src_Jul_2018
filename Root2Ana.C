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
	const double DISP=112; // unit??
	const double LOF=60.74; //m

	const double CALZ[3]={0, 1, 0};
	
	string sSet[2]={"PS_270_382", "RS_270_382"};
	
	string sRoot, sAna;
	StrtMesytec madc, mtdc, mqdcTOF, mqdcMCP;
	StrtS800 s800;
	StrtAna ana;

	double tPMT[8];
	double timeDet[2], egyDet[2];
	
	long long iEntry;
	int i, j, k, m, n, p, q;
	int iAna;
	int coin[3];
	int goodEvt, goodSi, goodMCP[2];
	double b;
	string setting;
	
	int runMin, runMax, runNum;
	cout<<"Input minimum and maximum numbers of run: ";
	cin>>runMin>>runMax;
	
	sAna="/home/kailong/ExpData/Jul2018/AnaData/ana-run-"+to_string(runMin)+"--"+to_string(runMax)+".root";
	TFile *fAna=new TFile(sAna.c_str(), "RECREATE");
	TTree *tAna=new TTree("tAna", "tree for data analysis");

	tAna->Branch("setting", &setting);
	tAna->Branch("ana", &ana, "tof[4]/D:amp[4]/D:tD[4][8][8]/D:egy[4][8]/D:xPla[4][2]/D:yPla[4][2]/D:xMCP[4][2]/D:yMCP[4][2]/D:delE[4][5]/D:tke[4]/D:beta[4]/D:gamma[4]/D:Z[4]/D:dZ[4]/D:brho[4][2]/D:AoQ[4][2]/D:Q[4][2]/D:ZmQ[4][2]/D:ZImQ[4][2]/D:A[4][2]/D:Araw[4][2]/D:Am2Q[4][2]/D:Am3Q[4][2]/D:Am2Z[4][2]/D:Am3Z[4][2]/D:dAm2Z[4][2]/D:dAm3Z[4][2]/D:Zi[4]/I:sig[4][2][2]/I:run/I");
	
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
		
		for(iEntry=0; iEntry<tData->GetEntries(); iEntry++)
		// for(iEntry=0; iEntry<10000; iEntry++)
		{
			tData->GetEntry(iEntry);
			if(s800.trig==1)
			{
				TRandom3 r(0);
				memset(&ana, 0, sizeof(ana));
				ana.run=runNum;
				goodEvt=0;
				
				for(iAna=0; iAna<4; iAna++)
				{
					memset(tPMT, 0, sizeof(tPMT));
					memset(timeDet, 0, sizeof(timeDet));
					memset(egyDet, 0, sizeof(egyDet));
					
					if(iAna==0) //For TAC+ADC+clock
					{
						for(i=0; i<8; i++)
						{
							j=2*i+1;
							k=i/4;
							if(madc.data[i]>LADC&&madc.data[i]<HADC&&mqdcTOF.data[j]>LQDCTOF&&mqdcTOF.data[j]<HQDC)
							{
								ana.sig[0][k][0]++;
								ana.sig[0][k][1]=10*ana.sig[0][k][1]+(i+1);
								tPMT[i]=CALADC[i]*(madc.data[i]+r.Uniform(-0.5,0.5));
								ana.egy[0][i]=mqdcTOF.data[j]+r.Uniform(-0.5,0.5);
							}
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
							if(madc.data[i]>LADC&&madc.data[i]<HADC&&mqdcTOF.data[j]>LQDCTOF&&mqdcTOF.data[j]<HQDC&&mqdcTOF.data[k]>LQDCTOF&&mqdcTOF.data[k]<HQDC)
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
								k=j+1;
							else
								k=2*j+17;
							m=2*j+1;
							n=j/4;
				
							if(mtdc.data[k]>LTDC&&mtdc.data[k]<HTDC&&mqdcTOF.data[m]>LQDCTOF&&mqdcTOF.data[m]<HQDC)
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
						ana.tof[iAna]=CALTOF[0]+CALTOF[1]*ana.tof[iAna];
						for(j=0; j<2; j++)
							if(ana.sig[iAna][j][0]==4)
							{
								k=4*j;
								ana.xPla[iAna][j]=(tPMT[2+k]+tPMT[3+k]-tPMT[0+k]-tPMT[1+k])/2;
								ana.yPla[iAna][j]=(tPMT[0+k]+tPMT[3+k]-tPMT[1+k]-tPMT[2+k])/2;
							}
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
								ana.xMCP[iAna][k]=1.0*(mqdcMCP.data[0+4*k]-LQDCMCP[0+4*k]+mqdcMCP.data[3+4*k]-LQDCMCP[3+4*k]-mqdcMCP.data[1+4*k]-LQDCMCP[1+4*k]-mqdcMCP.data[2+4*k]-LQDCMCP[2+4*k])/1.0/(mqdcMCP.data[0+4*k]-LQDCMCP[0+4*k]+mqdcMCP.data[3+4*k]-LQDCMCP[3+4*k]+mqdcMCP.data[1+4*k]-LQDCMCP[1+4*k]+mqdcMCP.data[2+4*k]-LQDCMCP[2+4*k]);
								
								ana.xMCP[iAna][k]=CALXMCP[k][0]+CALXMCP[k][1]*ana.xMCP[0][k]+CALXMCP[k][2]*pow(ana.xMCP[0][k],2)+CALXMCP[k][3]*pow(ana.xMCP[0][k],3);
								
								ana.yMCP[iAna][k]=1.0*(mqdcMCP.data[1+4*k]-LQDCMCP[1+4*k]+mqdcMCP.data[3+4*k]-LQDCMCP[3+4*k]-mqdcMCP.data[0+4*k]-LQDCMCP[0+4*k]-mqdcMCP.data[2+4*k]-LQDCMCP[2+4*k])/1.0/(mqdcMCP.data[0+4*k]-LQDCMCP[0+4*k]+mqdcMCP.data[3+4*k]-LQDCMCP[3+4*k]+mqdcMCP.data[1+4*k]-LQDCMCP[1+4*k]+mqdcMCP.data[2+4*k]-LQDCMCP[2+4*k]);
								
								ana.yMCP[iAna][k]=CALYMCP[k][0]+CALYMCP[k][1]*ana.yMCP[0][k]+CALYMCP[k][2]*pow(ana.yMCP[0][k],2)+CALYMCP[k][3]*pow(ana.yMCP[0][k],3);
							}

						ana.tke[iAna]=0;
						goodSi=0;
						for(q=0; q<5; q++)
							if(s800.pin[q]>100&&s800.pin[q]<4000)
							{
								ana.delE[iAna][q]=CALPIN[q][0]+CALPIN[q][1]*s800.pin[q];
								ana.tke[iAna]+=ana.delE[iAna][q];
								goodSi++;
							}
						if(s800.pin[0]>100&&s800.pin[0]<4000)
							ana.tke[iAna]+=(CALPIN[5][0]+CALPIN[5][1]*s800.pin[iAna]); //consider the absorption effect of material in front of Si detectors
							
						b=LOF/ana.tof[iAna]/0.299792458;
						if(b>0&&b<1)
						{
							ana.beta[iAna]=b;
							ana.gamma[iAna]=1/sqrt(1-b*b);
							if(s800.pin[0]>100&&s800.pin[0]<4000&&goodSi>1)
							{
								ana.Z[iAna]=ana.delE[iAna][0]/sqrt(1/b/b*log(5930.0/(1/b/b-1))-1);
								ana.Z[iAna]=CALZ[0]+CALZ[1]*ana.Z[iAna]+CALZ[2]*pow(ana.Z[iAna],2);
								ana.Zi[iAna]=TMath::Nint(ana.Z[iAna]);
								ana.dZ[iAna]=ana.Z[iAna]-ana.Zi[iAna];
								for(k=0; k<2; k++)
									if(goodMCP[k]==4)
									{
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
							}
						}
						
						if(b>0&&b<1&&s800.pin[0]>100&&s800.pin[0]<4000&&goodSi>1&&goodMCP[0]==4)
							goodEvt++;
					}
				}
				
				if(goodEvt>1)
					tAna->Fill();
			}	
		}//end of whole tree
		fRoot->Close();
	}
	fAna->cd();
	fAna->Write();
	fAna->Close();	
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