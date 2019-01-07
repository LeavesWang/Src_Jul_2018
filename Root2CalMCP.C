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

struct StrtAna
{
	double xMCP[2];
	double yMCP[2];
};

void Root2CalMCP()
{	
	const int QdcUp=3840;
	const int QdcMcpLow[8]={1000, 760, 810, 850, 790, 790, 780, 790};
	string sRoot, sCalMcp;
	StrtMesytec mqdcMCP;
	int run;
	double xMCP[2], yMCP[2];
	double calQdcMcp[8];
	
	int iEntry;
	int i, j, nMCP[2];

	int runMin, runMax, runNum;
	cout<<"Input minimum and maximum numbers of run: ";
	cin>>runMin>>runMax;
	
	sCalMcp="/home/kailong/ExpData/Jul2018/CalMCP/calMCP-run-"+to_string(runMin)+"--"+to_string(runMax)+".root";
	TFile *fCalMcp=new TFile(sCalMcp.c_str(), "RECREATE");
	TTree *tCalMcp=new TTree("tCalMcp", "tree for calibrating MCP");
	tCalMcp->Branch("run", &run, "run/I");
	tCalMcp->Branch("xMCP", xMCP, "xMCP[2]/D");
	tCalMcp->Branch("yMCP", yMCP, "yMCP[2]/D");
		
	ostringstream ssRun;	
	for(runNum=runMin; runNum<=runMax; runNum++)
	{
		ssRun.str("");
		ssRun<<setw(4)<<setfill('0')<<runNum;
		
		sRoot="/home/kailong/ExpData/Jul2018/RootData/run-"+ssRun.str()+"-00.root";
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
		
		for(iEntry=0; iEntry<tData->GetEntries(); iEntry++)
		{
			tData->GetEntry(iEntry);
			TRandom3 r(0);
			run=runNum;
			memset(xMCP, 0, sizeof(xMCP));
			memset(yMCP, 0, sizeof(yMCP));				
			memset(nMCP, 0, sizeof(nMCP));
			for(i=0; i<8; i++)
				if(mqdcMCP.data[i]>QdcMcpLow[i]&&mqdcMCP.data[i]<QdcUp)
				{
					j=i/4;
					nMCP[j]++;
					calQdcMcp[i]=mqdcMCP.data[i]-QdcMcpLow[i]+r.Uniform(-0.5,0.5);
				}
			for(i=0; i<2; i++)
				if(nMCP[i]==4)
				{
					xMCP[i]=(calQdcMcp[0+4*i]+calQdcMcp[3+4*i]-calQdcMcp[1+4*i]-calQdcMcp[2+4*i])/(calQdcMcp[0+4*i]+calQdcMcp[1+4*i]+calQdcMcp[2+4*i]+calQdcMcp[3+4*i]);
					yMCP[i]=(calQdcMcp[1+4*i]+calQdcMcp[3+4*i]-calQdcMcp[0+4*i]-calQdcMcp[2+4*i])/(calQdcMcp[0+4*i]+calQdcMcp[1+4*i]+calQdcMcp[2+4*i]+calQdcMcp[3+4*i]);
				}
			if(nMCP[0]==4 || nMCP[1]==4)
				tCalMcp->Fill();	
		}//end of whole tree
		fRoot->Close();
	}//end of runs
	fCalMcp->cd();
	tCalMcp->Write();
	fCalMcp->Close();
}//end of whole function

#ifndef __CINT__
void StandaloneApplication(int argc, char** argv)
{
	Root2CalMCP();
}

int main(int argc, char** argv)
{
	TApplication app("ROOT Application", &argc, argv);
	StandaloneApplication(app.Argc(), app.Argv());
	// app.Run();
	return 0;
}
#endif
