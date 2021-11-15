
//base
#include "hades.h"
#include "hspectrometer.h"
#include "hdetector.h"
#include "hruntimedb.h"
#include "htask.h"
#include "hevent.h"
#include "hcategory.h"
#include "hdst.h"
#include "htime.h"
#include "hsrckeeper.h"

//tasksets
#include "hstarttaskset.h"
#include "hrichtaskset.h"
#include "hrich700taskset.h"
#include "hmdctaskset.h"
#include "htoftaskset.h"
#include "hrpctaskset.h"
#include "hshowertaskset.h"
#include "hemctaskset.h"
#include "hwalltaskset.h"
#include "hsplinetaskset.h"

//finders, fillers, cleaners
#include "hparticlevertexfind.h"
#include "hparticlecandfiller.h"
#include "hparticletrackcleaner.h"
#include "hparticleevtinfofiller.h"
#include "hparticlebt.h"
#include "hparticlet0reco.h"
#include "hqamaker.h"

//defs
#include "haddef.h"
#include "richdef.h"
#include "hmdcdef.h"
#include "hmdctrackddef.h"
#include "hmdctrackgdef.h"
#include "showerdef.h"
#include "emcdef.h"
#include "rpcdef.h"
#include "tofdef.h"
#include "walldef.h"

#include "hrich700digipar.h"
#include "hrich700ringfinderpar.h"

//containers
#include "hmdcsetup.h"
#include "hmagnetpar.h"
#include "hmdclayercorrpar.h"
#include "hmdcdigitpar.h"
#include "hmdcdedx2maker.h"
#include "hmdcdigitizer.h"
#include "hrich700digitizer.h"      // rich700
#include "hmdctrackdset.h"
#include "hmdc12fit.h"
#include "hmetamatchF2.h"
#include "hstart2hitfsim.h"

//ROOT
#include "TSystem.h"
#include "TROOT.h"
#include "TString.h"
#include "TStopwatch.h"
#include "TDatime.h"
#include "TRegexp.h"

//standard c++
#include <iostream>
#include <cstdio>

using namespace std;





Int_t analysisDST_protons(TString inFile, TString outdir,Int_t nEvents=1000000, Int_t startEvt=0) //Jochen: nEvents=1 ...
{
    new Hades;
    TStopwatch timer;
    gHades->setTreeBufferSize(8000);
    gHades->makeCounter(10);
    gHades->setBeamTimeID(Particle::kMar19);
    gHades->getSrcKeeper()->addSourceFile("analysisDST_protons.cc");
    gHades->getSrcKeeper()->addSourceFile("treeFilter.h");
    gHades->getSrcKeeper()->addSourceFile("sendScript_SL.sh");


    //####################################################################
    //######################## CONFIGURATION #############################
    printf("Setting configuration...+++\n");

    HRuntimeDb *rtdb = gHades -> getRuntimeDb();
    //Int_t refId = 16000; // mar19 highfield;
    Int_t refId = 16001; // mar19 mediumfield;
    //Int_t refId = 16002; // mar19 lowfield;
    Bool_t HighMult = kTRUE;  // kTRUE : ag1650ag,ag4500ag , kFALSE : pp,pAg


    TString beamtime     ="mar19"; // used for spectrometer setup and track params
    TString beamtimeTrack="mar19"; // highmult = apr12, lowmult = aug14 used for  track params


    //-------------- Default Settings for the File names -----------------

    TString baseDir = outdir;
    if(!baseDir.EndsWith("/")) baseDir+="/";
    TString outDir   = baseDir;
    TString outDirQA = outDir+"qa/";

    TString outFileSuffix = "_dst_" + beamtime + ".root";

    TString asciiParFile = "./params/Geometry_MAR19SIM_v4.txt";
    //TString asciiParFile = "./params/dEdx2param_vladimir_03042021.txt";
    //TString asciiParFile = "./params/emc_par.txt";
    //TString rootParFile  = "./params/allParam_mar19_sim_run_16000_gen2_22102019.root";
    //TString rootParFile  = "./params/allParam_mar19_sim_run_16000_gen2_19022020.root";
    //TString rootParFile  = "./params/allParam_mar19_sim_run_16000_gen4_25112020.root";
    TString rootParFile  = "./params/allParam_mar19_sim_run_16000_gen5_06052021.root";
    //TString paramSource  = "ascii,oracle"; // root, ascii, oracle
    TString paramSource  = "root"; // root, ascii, oracle
    //TString paramRelease = "MAR19SIM_v2"; // 22102019
    //TString paramRelease = "MAR19SIM_v3"; // 19022020
    //TString paramRelease = "MAR19SIM_v4"; // 25112020
    TString paramRelease = "MAR19SIM_v5"; // 03042021
    //TString paramRelease = "now"; //

    Float_t scaleDeltaElectrons=0.69; // Au+Au to Ag+Ag (47*47 / 79*79)
    Bool_t kParamFile        = kFALSE;
    Bool_t doExtendedFit     = kTRUE; // switch on/off fit for initial params of segment fitter (10 x slower!)
    Bool_t doMetaMatch       = kFALSE;  // default : kTRUE, kFALSE switch off metamatch in clusterfinder
    Bool_t useOffVertex      = kTRUE;  // default : kTRUE,  kTRUE=use off vertex procedure  (apr12 gen8:kFALSE, gen9:kTRUE)
    Bool_t doMetaMatchScale  = kTRUE;
    Bool_t useWireStat       = kTRUE;
    Float_t metaScale        = 2.0;
    Bool_t doTree            = kFALSE; //default:: kTRUE (Jochen), kFALSE = use the Robert's version
    Bool_t doMDCDeltaElectron = kTRUE;
    Bool_t doRICHDeltaElectron= kFALSE; // kFALSE : suppress inserted delta electrons
    Bool_t doReverseField     = kFALSE; // kTRUE : switch magnet polarity

    //--------------------------------------------------------------------
    //------standard output for dst production---------------------------- 
    //--------------------------------------------------------------------
    // switch off unwanted categories in the poutput
    Cat_t notPersistentCatAll[] =
    {
	//catRich, catMdc, catShower, catTof, catTracks,
	catRichRaw,
	//catRichHitHdr,  // rich700
    catRichTrack,
	//catRichDirClus,
	//catRichHit,
	//catRichCal,
	catMdcRaw,
	catMdcCal1,   // changed
	catMdcCal2, catMdcHit,
	catMdcClusInf,
	catMdcSeg,
	catMdcTrkCand, catMdcRawEventHeader,

	catEmcRaw,
	//catEmcCal,
    //catEmcCluster

	catTofRaw,
	catTofHit,     // changed
	catTofCluster, // changed

	catRpcRaw,
	catRpcHit,      // changed
	catRpcCluster,  // changed
	catRpcCal,

	catRKTrackB, catSplineTrack,
	catMetaMatch,
	//catParticleCandidate, catParticleEvtInfo,
    catParticleMdc,
	catWallRaw, catWallOneHit, catWallCal,

	catMdcGeantRaw,catTofGeantRaw,catRpcGeantRaw,catEmcGeantRaw,catWallGeantRaw,
	catStartGeantRaw,catRichGeantRaw,catRichGeantRaw+1,catRichGeantRaw+2

    };

    //catMdcHit commented, same way as for the Robert case

        Cat_t notPersistentCat[] =
    {
    //catRich, catMdc, catShower, catTof, catTracks,
    catRichRaw,
    //catRichHitHdr,  // rich700
    catRichTrack,
    //catRichDirClus,
    //catRichHit,
    //catRichCal,
    catMdcRaw,
    catMdcCal1,   // changed
    catMdcCal2, //catMdcHit,
    catMdcClusInf,
    catMdcSeg,
    catMdcTrkCand, catMdcRawEventHeader,

    catEmcRaw,
    //catEmcCal,
    //catEmcCluster

    catTofRaw,
    catTofHit,     // changed
    catTofCluster, // changed

    catRpcRaw,
    catRpcHit,      // changed
    catRpcCluster,  // changed
    catRpcCal,

    catRKTrackB, catSplineTrack,
    catMetaMatch,
    //catParticleCandidate, catParticleEvtInfo,
    catParticleMdc,
    catWallRaw, catWallOneHit, catWallCal,

    catMdcGeantRaw,catTofGeantRaw,catRpcGeantRaw,catEmcGeantRaw,catWallGeantRaw,
    catStartGeantRaw,catRichGeantRaw,catRichGeantRaw+1,catRichGeantRaw+2

    };

    //--------------------------------------------------------------------

    //####################################################################
    //####################################################################


    //------------- Operations on the filenames --------------------------
    TString rootSuffix =".root";
    TString nFile;    // needed to build outputfile
    TString dirname;  // needed for datasource
    TString filename; // needed for datasource
    TString outFile;

    Int_t sourcetype = 3; // 3 root source, 4 merge source


    cout<<"infile: "<<inFile<<endl;

    if(inFile.Contains(",") || sourcetype == 4 ){ // comma seperated list for geant merge source
	sourcetype = 4;
	inFile.ReplaceAll(" ","");
	TObjArray* ar = inFile.Tokenize(",");
	TString firstfile;
        TString f;
	/*
	if(ar){
	    if(ar->GetEntries()>0) {
		firstfile = ((TObjString*)ar->At(0))->GetString();
	    }
	    delete ar;
	}
	*/

	//------------------------------------------------------------
	// build out file name. In case of mult inputs
	// files can be recycled. In this case taking the
	// first file to derive the outfile name is not sufficient.
	// Instead we will take a concat of file numbers of all inputs and
        // attach it to the first file name to prevent over writing
	TString fileNum   = "";
        TString firstFile = "";
	for(Int_t i=0;i<ar->GetEntries();i++)
	{

             f = ((TObjString*)ar->At(i))->GetString();


	     TString f1 = gSystem->BaseName(f);
	     TRegexp reg         ("_[01234456789]*_1.root");
	     TString result = f1(reg);

	     if(i==0)
	     {
		 firstfile  = f;
                 firstfile.ReplaceAll(result.Data(),"");

	     }

	     result.ReplaceAll("_1.root","");
             result.ReplaceAll("_","");

             fileNum=Form("%s_%s",fileNum.Data(),result.Data());
             cout<<firstfile<<" "<<fileNum<<endl;
	}

	firstfile = Form("%s%s_1.root",firstfile.Data(),fileNum.Data());
	//------------------------------------------------------------

	nFile     = gSystem->BaseName(firstfile.Data());
	filename  = inFile;
	dirname   = "";
    //hijacked from Robert, not sure what it does... (evntBuild should be 1 or not 1 for some criteria....)
    Int_t evtBuild = HTime::getEvtBuilderFileName(nFile,kFALSE);

    }  else {  // root source

	nFile     = gSystem->BaseName(inFile.Data());
	filename  = gSystem->BaseName(inFile.Data());
	dirname   = gSystem->DirName(inFile.Data());
    }
    if (nFile.EndsWith(rootSuffix)) nFile.ReplaceAll(rootSuffix,"");
    outFile  = outDir+nFile+outFileSuffix;
    outFile.ReplaceAll("//", "/");


    if(gSystem->AccessPathName(outDir.Data()) != 0){
	cout<<"Creating output dir :"<<outDir.Data()<<endl;
	gSystem->Exec(Form("mkdir -p %s",outDir.Data()));
    }
    if(gSystem->AccessPathName(outDirQA.Data()) != 0){
	cout<<"Creating output qadir :"<<outDirQA.Data()<<endl;
	gSystem->Exec(Form("mkdir -p %s",outDirQA.Data()));
    }
    //--------------------------------------------------------------------




    Int_t mdcMods[6][4]=
    { {1,1,1,1},
    {1,1,1,1},
    {1,1,1,1},
    {1,1,1,1},
    {1,1,1,1},
    {1,1,1,1} };

    // recommendations from Vladimir+Olga
    // according to params from 28.04.2011
    Int_t nLayers[6][4] = {
	{6,6,5,6},
	{6,6,5,6},
	{6,6,5,6},
	{6,6,5,6},
	{6,6,5,6},
	{6,6,5,6} };
    Int_t nLevel[4] = {10,50000,10,5000};

    HDst::setupSpectrometer(beamtime,mdcMods,"start,rich,mdc,tof,rpc,emc,wall");
    // beamtime mdcMods_apr12, mdcMods_full
    // Int_t mdcset[6][4] setup mdc. If not used put NULL (default).
    // if not NULL it will overwrite settings given by beamtime
    // detectors (default)= rich,mdc,tof,rpc,shower,wall,tbox,start


    HDst::setupParameterSources(paramSource,asciiParFile,rootParFile,paramRelease);
    //HDst::setupParameterSources("oracle",asciiParFile,rootParFile,"now"); // use to create param file
    // parsource = oracle,ascii,root (order matters)
    // if source is "ascii" a ascii param file has to provided
    // if source is "root" a root param file has to provided
    // The histDate paramter (default "now") is used wit the oracle source

    cout<<"dir="<<dirname<<", file="<<filename<<", refid="<<refId<<endl;

    HDst::setDataSource(sourcetype,dirname,filename,refId); // Int_t sourceType,TString inDir,TString inFile,Int_t refId, TString eventbuilder"
    // datasource 0 = hld, 1 = hldgrep 2 = hldremote, 3 root, 4 geantmerge
    // like "lxhadeb02.gsi.de"  needed by dataosoure = 2
    // inputDir needed by datasoure = 1,2
    // inputFile needed by datasoure = 1,3
    // for datasource 4 inputFile is a comma seprated list
    // "file1_with_path,file2_with_path,file3_with_path"


    if(kParamFile) {

	TDatime time;
        TString paramfilename= Form("allParam_mar19_sim_run_%i_gen5_%02i%02i%i",refId,time.GetDay(),time.GetMonth(),time.GetYear());  // without .root

	if(gSystem->AccessPathName(Form("%s.root",paramfilename.Data())) == 0){
	    gSystem->Exec(Form("rm -f %s.root",paramfilename.Data()));
	}
	if(gSystem->AccessPathName(Form("%s.log",paramfilename.Data())) == 0){
	    gSystem->Exec(Form("rm -f %s.log" ,paramfilename.Data()));
	}

	if (!rtdb->makeParamFile(Form("%s.root",paramfilename.Data()),"mar19sim","10-MAR-2019 07:00:00","31-MAR-2019 20:00:01")) {
	    delete gHades;
	    exit(1);
	}
    }

    //--------------------------------------------------------------------
    // ----------- Build TASK SETS (using H***TaskSet::make) -------------
    HStartTaskSet        *startTaskSet        = new HStartTaskSet();
    HRich700TaskSet      *richTaskSet         = new HRich700TaskSet();
    HRpcTaskSet          *rpcTaskSet          = new HRpcTaskSet();
    HEmcTaskSet          *emcTaskSet          = new HEmcTaskSet();
    HTofTaskSet          *tofTaskSet          = new HTofTaskSet();
    HWallTaskSet         *wallTaskSet         = new HWallTaskSet();
    HMdcTaskSet          *mdcTaskSet          = new HMdcTaskSet();
    //    mdcTaskSet->setVersionDeDx(1); // 0 = no dEdx, 1 = HMdcDeDx2

    //HMagnetPar* magnet = (HMagnetPar*)rtdb->getContainer("MagnetPar");
    /*
    rtdb->initContainers(refId);
   // magnet->setStatic();
    //magnet->setCurrent(3200);

    if(doReverseField){
        Int_t current = magnet->getCurrent();
	magnet->setCurrent(-1*current);
    }
    magnet->printParams();
    */
    HMdcSetup* mysetup    = (HMdcSetup*)rtdb->getContainer("MdcSetup");
    //HMdcDigitPar* digipar = (HMdcDigitPar*)rtdb->getContainer("MdcDigitPar"); // simulate impact of 1 sparking layer switched off
    rtdb->initContainers(refId);
    mysetup->setStatic();

    //mysetup->getMdcCommonSet()->setIsSimulation(1);                 // fit
    //mysetup->getMdcCommonSet()->setAnalysisLevel(1);                // fit
    //mysetup->getMdcTrackFinderSet()->setIsCoilOff(kFALSE);          // field is on
    //mysetup->getMdcTrackFinderSet()->setNLayers(nLayers[0]);
    //mysetup->getMdcTrackFinderSet()->setNLevel(nLevel);
    //mysetup->getMdc12FitSet()->setMdc12FitSet(2,1,0,kFALSE,kFALSE); // tuned fitter, seg

    HTask *startTasks         = startTaskSet       ->make("simulation");
    HTask *richTasks          = richTaskSet        ->make("simulation"); //"NORINGFINDER"
    HTask *tofTasks           = tofTaskSet         ->make("simulation");
    HTask *wallTasks          = wallTaskSet        ->make("simulation");
    HTask *rpcTasks           = rpcTaskSet         ->make("simulation");
    HTask *emcTasks           = emcTaskSet         ->make("simulation");
    HTask *mdcTasks           = mdcTaskSet         ->make("rtdb","");

    HMdcDigitizer* digi = mdcTaskSet->getDigitizer();
    if(digi){
	//digi->setDeltaElectronUse(doMDCDeltaElectron,kFALSE,109,-750.,600.,20.,scaleDeltaElectrons);
	digi->setDeltaElectronUse(doMDCDeltaElectron,kFALSE,109,-250.,700.,20.,scaleDeltaElectrons); // (Bool_t use, Bool_t useDeltaMomSel=kFALSE, Int_t ionId=109,Float_t t1min=-950.,Float_t t1max=400.,Float_t momCut=20.,Float_t probDelta=2
	digi->setDeltaElectronMinMomCut(2.,2.,4.5,2.,2.,4.5);  // take care of glass mirrors in sec 2+5
	digi->setTimeCutUse(kTRUE);
    }
    HRich700Digitizer* richdigi = HRich700Digitizer::getDigitizer();
    if(richdigi){
	richdigi->setDeltaElectronUse(doRICHDeltaElectron,kFALSE,109,20.,0.66*scaleDeltaElectrons); // 1 - prob  0.5 Mus RICH / 1.5 mus MDC
	richdigi->setDeltaElectronMinMomCut(0.,0.,0.,0.,0.,0.);
    }


    //----------------SPLINE and RUNGE TACKING----------------------------------------
    HSplineTaskSet         *splineTaskSet       = new HSplineTaskSet("","");
    HTask *splineTasks     = splineTaskSet      ->make("","spline,runge");

    HParticleCandFiller    *pParticleCandFiller = new HParticleCandFiller   ("particlecandfiller","particlecandfiller","");
    pParticleCandFiller->setFillMdc(kFALSE); // new : testing close pair
    //pParticleCandFiller->setDoMDCdEdxCorr(kFALSE); // new dev

    HParticleTrackCleaner  *pParticleCleaner    = new HParticleTrackCleaner ("particlecleaner"   ,"particlecleaner");
    HParticleVertexFind    *pParticleVertexFind = new HParticleVertexFind   ("particlevertexfind","particlevertexfind",kTRUE);
    HParticleEvtInfoFiller *pParticleEvtInfo    = new HParticleEvtInfoFiller("particleevtinfo"   ,"particleevtinfo",beamtime);



    //----------------------- Quality Assessment -------------------------
    HQAMaker *qaMaker =0;
    if (!outDirQA.IsNull())
    {
	qaMaker = new HQAMaker("qamaker","qamaker");
	//qaMaker->setUseSlowPar(kFALSE);
	qaMaker->setOutputDir((Text_t *)outDirQA.Data());
	//qaMaker->setPSFileName((Text_t *)hldFile.Data());
	qaMaker->setSamplingRate(1);
	qaMaker->setIntervalSize(50);
    }



    //------------------------ Master task set ---------------------------
    HTaskSet *masterTaskSet = gHades->getTaskSet("all");

    masterTaskSet->add(startTasks);
    masterTaskSet->add(tofTasks);
    masterTaskSet->add(wallTasks);
    masterTaskSet->add(rpcTasks);

    masterTaskSet->add(richTasks);

    masterTaskSet->add(emcTasks);
    masterTaskSet->add(mdcTasks);
    masterTaskSet->add(splineTasks);

    masterTaskSet->add(pParticleCandFiller);
    masterTaskSet->add(pParticleCleaner);
    masterTaskSet->add(pParticleVertexFind); // run after track cleaning
    masterTaskSet->add(pParticleEvtInfo);

    masterTaskSet->add(new HParticleT0Reco("T0","T0",beamtime));

    //--------------------------------------------------------------------
    // write in parallel filtered event files
    addFilter(masterTaskSet,inFile,outdir) ;  // treeFilter.h
    //--------------------------------------------------------------------

    //if (qaMaker) masterTaskSet->add(qaMaker); //commented by Robert...

    HMdcTrackDSet::setTrackParam(beamtimeTrack);
    //HMdcTrackDSet::setFindOffVertTrkFlag(useOffVertex);
    if(!doMetaMatch)HMdcTrackDSet::setMetaMatchFlag(kFALSE,kFALSE);  //do not user meta match in clusterfinder
    if(doMetaMatchScale)HMetaMatchF2::setScaleCut(metaScale,metaScale,metaScale); // (tof,rpc,shower) increase matching window, but do not change normalization of MetaQA
    if(useWireStat) HMdcCalibrater1::setUseWireStat(kTRUE);

    HStart2HitFSim* starthitf = HStart2HitFSim::getHitFinder() ;
    if(starthitf) starthitf->setResolution(0.1);    // 100 ps start res

    //--------------------------------------------------------------------
    // find best initial params for segment fit (takes long!)
    if(doExtendedFit) {
 	HMdcTrackDSet::setCalcInitialValue(1);  // -  1 : for all clusters 2 : for not fitted clusters
    }
    //--------------------------------------------------------------------

    if (!gHades->init()){
	Error("init()","Hades did not initialize ... once again");
	exit(1);
    }


    //--------------------------------------------------------------------
    //----------------- Set not persistent categories --------------------
    HEvent *event = gHades->getCurrentEvent();

    /*for(UInt_t i=0;i<sizeof(notPersistentCat)/sizeof(Cat_t);i++){
	HCategory *cat = ((HCategory *)event->getCategory(notPersistentCat[i]));
	if(cat)cat->setPersistency(kFALSE);
    }*/

    //some Robert magic here
    if(evtBuild == 1) {
        for(UInt_t i=0;i<sizeof(notPersistentCatAll)/sizeof(Cat_t);i++){
            cat = ((HCategory *)event->getCategory(notPersistentCatAll[i]));
            if(cat)cat->setPersistency(kFALSE);
        }
    } 
    else {
        for(UInt_t i=0;i<sizeof(notPersistentCat)/sizeof(Cat_t);i++){
            cat = ((HCategory *)event->getCategory(notPersistentCat[i]));
            if(cat)cat->setPersistency(kFALSE);
        }
    }
    //--------------------------------------------------------------------

    //And most tricky Robert part... this thing (kFALSE for doTree):

    if(doTree){
	// output file
	gHades->setOutputFile((Text_t*)outFile.Data(),"RECREATE","Test",2);
	gHades->makeTree();
    }
    else{
        cout << "Event loop starts here!" << endl; //to check if Robert part kicked in
    }


    Int_t nProcessed = gHades->eventLoop(nEvents,startEvt);
    printf("Events processed: %i\n",nProcessed);

    cout<<"--Input file      : "<<inFile  <<endl;
    cout<<"--QA directory is : "<<outDirQA<<endl;
    cout<<"--Output file is  : "<<outFile <<endl;

    printf("Real time: %f\n",timer.RealTime());
    printf("Cpu time: %f\n",timer.CpuTime());
    if (nProcessed) printf("Performance: %f s/ev %f ev/s\n",timer.CpuTime()/nProcessed,nProcessed/timer.CpuTime());

    if(kParamFile) rtdb->saveOutput();

    delete gHades;
    timer.Stop();

    return 0;

}

#ifndef __CINT__
int main(int argc, char **argv)
{
    TROOT AnalysisDST("AnalysisDST","compiled analysisDST macros");


    TString nevents,startevent;
    switch (argc)
    {
    case 3:
	return analysisDST(TString(argv[1]),TString(argv[2])); // inputfile + outdir
	break;
    case 4:  // inputfile + outdir + nevents
	nevents=argv[3];

	return analysisDST(TString(argv[1]),TString(argv[2]),nevents.Atoi());
	break;
	// inputfile + nevents + startevent
    case 5: // inputfile + outdir + nevents + startevent
	nevents   =argv[3];
	startevent=argv[4];
	return analysisDST(TString(argv[1]),TString(argv[2]),nevents.Atoi(),startevent.Atoi());
	break;
    default:
	cout<<"usage : "<<argv[0]<<" inputfile outputdir [nevents] [startevent]"<<endl;
	return 0;
    }
}
#endif
