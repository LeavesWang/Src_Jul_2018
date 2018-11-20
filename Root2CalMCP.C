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

void Root2CalMCP()
{	
	const int HQDC=3840;
	const int LQDC[8]={726, 730, 745, 742, 752,758,764,761};
	string sRoot, sCalMcp;
	StrtMesytec mqdcMCP;
	double xMCP[2], yMCP[2];
	
	int iEntry;
	int i, j, nMCP[2];

	int runMin, runMax, runNum;
	cout<<"Input minimum and maximum numbers of run: ";
	cin>>runMin>>runMax;
	ostringstream ssRun;	
	for(runNum=runMin; runNum<=runMax; runNum++)
	{
		ssRun.str("");
		ssRun<<setw(4)<<setfill('0')<<runNum;
		
		sRoot="/home/kailong/ExpData/Jul2018/RootData/run-"+ssRun.str()+"-00.root";
		sCalMcp="/home/kailong/ExpData/Jul2018/CalMCP/calMCP-run-"+ssRun.str()+"-00.root";
		printf("\n**********Now converting %s to %s!**********\n\n", sRoot.c_str(), sCalMcp.c_str());
    	
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
		
		TFile *fCalMcp=new TFile(sCalMcp.c_str(), "RECREATE");
		TTree *tCalMcp=new TTree("tCalMcp", "tree for calibrating MCP");
		tCalMcp->Branch("xMCP", xMCP, "xMCP[2]/D");
		tCalMcp->Branch("yMCP", yMCP, "yMCP[2]/D");
		
		for(iEntry=0; iEntry<tData->GetEntries(); iEntry++)
		{
			tData->GetEntry(iEntry);
			TRandom3 r(0);
		
			memset(xMCP, 0, sizeof(xMCP));
			memset(yMCP, 0, sizeof(yMCP));				
			memset(nMCP, 0, sizeof(nMCP));
			for(i=0; i<8; i++)
				if(mqdcMCP.data[i]>LQDC[i]&&mqdcMCP.data[i]<HQDC)
				{
					j=i/4;
					nMCP[j]++;
				}
			for(i=0; i<2; i++)
				if(nMCP[i]==4)
				{
					xMCP[i]=1.0*(mqdcMCP.data[0+4*i]+mqdcMCP.data[3+4*i]-mqdcMCP.data[1+4*i]-mqdcMCP.data[2+4*i])/1.0/(mqdcMCP.data[0+4*i]+mqdcMCP.data[3+4*i]+mqdcMCP.data[1+4*i]+mqdcMCP.data[2+4*i]);
					yMCP[i]=1.0*(mqdcMCP.data[1+4*i]+mqdcMCP.data[3+4*i]-mqdcMCP.data[0+4*i]-mqdcMCP.data[2+4*i])/1.0/(mqdcMCP.data[0+4*i]+mqdcMCP.data[3+4*i]+mqdcMCP.data[1+4*i]+mqdcMCP.data[2+4*i]);
				}
			if(nMCP[0]==4 || nMCP[1]==4)
				tCalMcp->Fill();	
		}//end of whole tree
		fCalMcp->Write();
		fCalMcp->Close();
		fRoot->Close();
	}//end of runs
}//end of whole function

#ifndef __CINT__
/*
void StandaloneApplication(int argc, char** argv)
{
	Root2CalMCP();
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
	Root2CalMCP();
	return 0;
}
#endif
