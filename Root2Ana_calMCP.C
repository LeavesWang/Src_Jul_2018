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
	double xMCP[2];
	double yMCP[2];
};

void Root2Ana_calMCP()
{	
	const int LQDC=700, HQDC=3839;
	string sRoot, sAna;
	StrtMesytec mqdcMCP;
	double xMCP[2], yMCP[2];
	
	int iEntry;
	int i, j, k, m, n, p, q, u, nMCP[2];

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
		sAna="/home/kailong/ExpData/Jul2018/AnaData/calMCP-run-"+ssRun.str()+"-00.root";
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

		memset(&mqdcMCP, 0, sizeof(mqdcMCP));
		
		tData->SetBranchAddress("mqdcMCP", &mqdcMCP);		
		
		TFile *fAna=new TFile(sAna.c_str(), "RECREATE");
		TTree *tAna=new TTree("tAna", "tree for calibrating MCP");
		tAna->Branch("xMCP", xMCP, "xMCP[2]/D");
		tAna->Branch("yMCP", yMCP, "yMCP[2]/D");
		
		for(iEntry=0; iEntry<tData->GetEntries(); iEntry++)
		{
			tData->GetEntry(iEntry);
			TRandom3 r(0);
		
			memset(xMCP, 0, sizeof(xMCP));
			memset(yMCP, 0, sizeof(yMCP));				
			memset(nMCP, 0, sizeof(nMCP));
			for(p=0; p<8; p++)
				if(mqdcMCP.data[p]>LQDC&&mqdcMCP.data[p]<HQDC)
				{
					m=p/4;
					nMCP[m]++;
				}
			for(k=0; k<2; k++)
				if(nMCP[k]==4)
				{
					xMCP[k]=1.0*(mqdcMCP.data[0+4*k]+mqdcMCP.data[3+4*k]-mqdcMCP.data[1+4*k]-mqdcMCP.data[2+4*k])/1.0/(mqdcMCP.data[0+4*k]+mqdcMCP.data[3+4*k]+mqdcMCP.data[1+4*k]+mqdcMCP.data[2+4*k]);
					yMCP[k]=1.0*(mqdcMCP.data[1+4*k]+mqdcMCP.data[3+4*k]-mqdcMCP.data[0+4*k]-mqdcMCP.data[2+4*k])/1.0/(mqdcMCP.data[0+4*k]+mqdcMCP.data[3+4*k]+mqdcMCP.data[1+4*k]+mqdcMCP.data[2+4*k]);
				}
			if(nMCP[0]==4 || nMCP[1]==4)
				tAna->Fill();	
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
	Root2Ana_calMCP();
	return 0;
}
#endif
