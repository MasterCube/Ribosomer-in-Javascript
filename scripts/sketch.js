// canvas size:
var canvas_width = 4900;
var canvas_height = 2000;
// execution time of draw loop:
var cycle_time_draw_loop = 0;
var rib_FOCUS = -1;
// buttons variables:
var buttonRPF;
var buttonRPFcolor;
var ButtonDefaultColor = 'rgb(255, 0, 0)'; //red
var ButtonColorGreen = 'rgb(0, 218, 0)'; //green

// buttons variable status:
var prolineStatus;
var EFPdepleteStatus;
var tunnelStatus;
var nonUniformmRNAabundanceStatus;
var nonUniformInitiationRatesStatus;

// species variable status:
var coliSpeciesStatus;
var yeastSpeciesStatus;
var chosenSpecies; // yeast or 'coli' or....

// dictionary from JSON files with hypoexponetial parameters for all 61 sense codons:
var yeastHypoExpParams, coliHypoExpParams;
// dictionary from JSON file with ADAT A to I editing enzyme silencing step 1 rate factors:
var adatSilencedDict;
// dictionary from JSON file with ribosome exit tunnel electrostatics axial forces
// for a window of 50 amino acids (position 1 is most carboxy-terminal, position 50 is N-terminal)
var tunnelElectrostaticsAF;

//let tRNA1pic;
// uncomment this to get a nice .png to draw a real tRNA
function preload() {
  //tRNA1pic = loadImage("http://localhost:8080/images/tRNAnormal.png");
  //tRNA1pic = loadImage("http://localhost:3000/tRNAnormal.png");
  //tRNA1pic = loadImage('./images/tRNAnormal.png');
  tRNA1pic = loadImage("http://localhost:8080/tRNAnormal.png");
  yeastHypoExpParams = loadJSON("http://localhost:8080/Data/dataJSONyeast.json");
  coliHypoExpParams = loadJSON("http://localhost:8080/Data/dataJSONcoli.json");
  tunnelElectrostaticsAF = loadJSON("http://localhost:8080/Data/tunnelElectrostaticsAF.json");
  adatSilencedDict = loadJSON("http://localhost:8080/Data/dataJSONadatFactors.json");
  uridine34SilencedDict = loadJSON("http://localhost:8080/Data/dataJSONu34Factors.json");
}
var picWidth = 68;  // 1700 x 790 px
var picHeight = 32;
var keyI;
var StoredGeneID = [];
var StoredmRNA = [];
var RiboSeqSPtruth = []; // this is an array that will eventually have 2D.
var RiboSeqRRTtruth = [];// idem
var RiboSeqFCsampled = []; // idem
var read_added_count = 0;
var readCount = [];
var readInitRate = [];
var riboRecruitScore = []; // is the readInitRate for each transcript but multiplied by 10 for 10% relative difference
var eventFlag = [];
var eventSnaps = [];
var iMax;
// store 36 event flags for 36 time intervals of 5 minutes:
for (i=0; i<36; i++){
  eventFlag[i] = 'false';
}
// initiate # event Snaps based on total simulation time and time time_between_snapshots:
var time_between_snapshots = 10; // 10 seconds as measured by time_sec.
var time_sec;
var total_sim_time = 108 // min
var snapshots_nb;
var footprintedFragmentsCount = 0;

var StoredTimes = [];
var listOfProlinePositions = [];
var listOfAAplusPositions = [];
var listOfAAminusPositions = [];

// Objects instance defined:
// lattice instances:
let lattice = [];
var transcriptCardinality;  // number of different transcript geneID
var transcriptFullCardinality; // total number of transcripts reads (add each geneID times each copy number)
var transcriptomeLength_nts = 0;
var transcriptomeLength_codons = 0;
var tr_displayed_set = [];
var riboPoolRatio;
var indicesListToSampleFrom = [];

// tRNAs instances:
let tRNAiso1_K;     // isoacceptor 1 tRNA for lysine
let tRNAiso2_K;     // isoacceptor 2 tRNA for Lysine
let tRNAsample = [];
var tRNAmolNumber = 10;
var tRNAlinearMaxAmplitude = 50; // max pixels amplitude per frame
var HALF_PI = 1.57;
var tRNAangularSpeed = HALF_PI/9.0; // 10 degrees per frame

// ribosomes variables/instantiations:
// free ribosomes instances:
let ribosome = []; // list of ribosomes objects
var freeRibosomeInitialNumber;  // initial number (at start) of free ribosomes (e.g., 10).
var freetype; // number of ribosomes of given type.
var initiatedtype; // number of ribosomes of given type.
var visibletype;
// Psite offset:
const riboFootprintLength = 33; // 33 nucleotides = 33 pixels.
const PsiteOffset = 12; // 4 codons (= 12 px) between 5' end footprint and P site of the ribosome (total footprint = 10 codons).
// this last constant PsiteOffset is never used...
// Translation rates for ribosomes:
// Ribosome Elongation Rate:
var RiboPosX = 415;
var x = RiboPosX;
var RiboElongationRate; // use 6:   6 codons/s = 18px/s for yeast

// Read the file with the ORFs:
var dropzone;
var tRNA;
var txt;
//var output;
//---------------------------------------------------------------------------
// genetic code 3D-array: 4 by 4 by 4 3D array with 64 elements:
var GenetiCode = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]];
//------------------------------------------------------------------------
var protCountTimePointTranscript = [];
var RPFCountTimePointTranscript = [];
var sequence = "ATGAAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCTAA";
var transcriptLength;
var peptide = [];
var codonParams = [[], [], []]; // 2D-array: 61 by 3
//------------------------------------------------------------------------
// Translation rates for ribosomes:
// Ribosome Elongation Rate:
var RiboPosX = 415;
var x = RiboPosX;
//var RiboElongationRate = 6; // use 6:   6 codons/s = 18px/s
// Initiation rate:
// average time between two initiations: 250 ms (4 initions per second):
var InitALPHA = 1;
// one initiation per time span of 60 s.
var InitBETA = 60000*3; // minimum 250 ms between two initiations on average
//*********************
// functions (callof)
function displayRPFOnOff(){
  // color the button in green/red
  // test button.style value:
  if (buttonRPF.style('background-color') == 'rgb(255, 0, 0)')
  {
    buttonRPF.style('background-color', ButtonColorGreen);
  }
  else
  {
    buttonRPF.style("background-color", ButtonDefaultColor);
  }
}

//*****************************************************************
// function setup
//*****************************************************************
function setup() {
  // put setup code here
  var canvasOne =  createCanvas(canvas_width, canvas_height); // 1920, 938
  canvasOne.parent('#canvasHere');
  //------------------------------------------------------------------------
  frameRate(60); // use 18 frames/s or 12/s // 6 frame/s = 167 ms between two frames.//12
  //------------------------------------------------------------------------
  //**** In this page, there will be only one button functionnality:
  // creation of buttons to activate the display of 31 nts RPF (protetced fragment reads):
  // DISPLAY RPF button
  buttonRPF = createButton("Display RPF (31 nts)");
  buttonRPF.parent("#canvasHere");
  buttonRPF.size(180, 40);// 200, 27
  buttonRPF.style('font-size', '16px');
  buttonRPF.position(1300, 140, 'relative');
  buttonRPF.style("background-color", ButtonDefaultColor);
  //****
  for (i=0; i<36; i++){
    //eventFlag[i] = 'false';
    // clear the local storage for timePoint# :
    if (getItem('timePoint#'+str(i)) != null){
      removeItem('timePoint#'+str(i));
    }
    if (getItem('protCountTimePoint#' + str(i)) != null){
      removeItem('protCountTimePoint#' + str(i));
    }
    if (getItem('rpfCountTimePoint#' + str(i)) != null){
      removeItem('rpfCountTimePoint#' + str(i));
    }
  }

  if (getItem('riboSeqSnapShotsNb') != null){
    removeItem('riboSeqSnapShotsNb');
  }
  if (getItem('footprintFragNb') != null){
    removeItem('footprintFragNb');
  }
  if (getItem('riboSeqFC') != null){
    removeItem('riboSeqFC');
  }
  var currentTime;

  snapshots_nb = floor(total_sim_time * 60 /time_between_snapshots);
  for (i=0; i<snapshots_nb; i++){
    eventSnaps[i] = 'false';
  }

  //console.log('yeastHypoParams["AAG"][0] = ', yeastHypoExpParams['AAG'][0]);
  //console.log('yeastHypoParams.length=', Object.keys(yeastHypoExpParams).length);

  // Retrieve the value of the average initiation rate for the ribosomes to initiate on the transcripts:
  if (getItem("initiationRateTime") != null){
    InitBETA = float(getItem("initiationRateTime"));
  }

  // Retrieve the ribosome pool ratio (# ribosomes / # transcripts):
  if (localStorage.getItem("riboPoolRatio") != null){
    riboPoolRatio = localStorage.getItem("riboPoolRatio");
  }

  // status of the used active buttons in the settings page:
  prolineStatus = getItem("prolineStatus");
  //console.log('proline button = ', prolineStatus);
  tunnelStatus = getItem("exitTunnelStatus");
  //console.log('tunnel electrostatic interaction button = ', tunnelStatus);
  EFPdepleteStatus = getItem("EFPdepleteStatus");
  //console.log('EFPdepleteStatus button = ', EFPdepleteStatus);
  nonUniformmRNAabundanceStatus = getItem("nonUniformmRNAabundanceStatus");
  // Retrieve the nonUniformmRNAabundanceStatus
  nonUniformInitiationRatesStatus = getItem("nonUniformInitiationRatesStatus");
  // Retrieve the nonUniformInitiationRatesStatus
  adatSilencedStatus = getItem("tRNAmultipleModificationStatus")
  // Retrieve the ADAT silencing status
  tRNAmodificationStatus = getItem("tRNAmodificationStatus")
  // status of the species used as expression vector in the protein synthesis similation:
  // use getItem and stuff... later
  yeastSpeciesStatus = getItem("yeastSpeciesStatus");
  coliSpeciesStatus = getItem("coliSpeciesStatus");
  if (coliSpeciesStatus=='true' && yeastSpeciesStatus=='true'){
    console.log('There is confusion in the species to use as expression vector. Two species were selected simulatneously.Please select only one species in the settings tab.');
    //text("Please select only one species in the settings tab.");
  }
  if (coliSpeciesStatus=='true'){
    chosenSpecies = 'coli';
  }
  if (yeastSpeciesStatus=='true'){
    chosenSpecies = 'yeast';
  }
  // tRNAs instantiations:
  tRNAiso1_K = new Trna();
  tRNAiso2_K = new Trna();
  for (let nb=0; nb <tRNAmolNumber; nb++){
    tRNAsample.push(new Trna());
  }
  // this generate a gamma distributed random value:
  alpha=1.0; //1.6438
  beta=250.0; //108.28
  rate_of_beta = 1.0/beta; // in ms minus 1.
  //console.log('random gamma 1 : ', jStat.gamma.sample(alpha, beta))
  //console.log('random gamma 2 : ', jStat.gamma.sample(alpha, beta))
  //console.log('random gamma 3 : ', jStat.gamma.sample(alpha, beta));
  //console.log('gamma distribution mean : ', jStat.gamma.mean(alpha, beta));
  //console.log('gamma distribution variance : ', jStat.gamma.variance(alpha, beta));
  //console.log('pdf for gamma of x=100 :', jStat.gamma.pdf(100, alpha, beta));
  //console.log('random exponential dist time 1 : ', jStat.exponential.sample(rate_of_beta));
  //console.log('random exponential dist time 2 : ', jStat.exponential.sample(rate_of_beta));
  //console.log('random exponential dist time 3 : ', jStat.exponential.sample(rate_of_beta));
  //console.log('mean of exponential dist times : ', jStat.exponential.mean(rate_of_beta));

  iMax = getItem("numbLines");
  for (i=0; i<iMax;i++){
    const keyOne = getItem(str(i));
    StoredGeneID[i] = keyOne
    StoredmRNA[i] = getItem(keyOne);
    if(nonUniformmRNAabundanceStatus=='true'){
      readCount[i] = getItem('readCopiesOf_'+str(i));
    }
    else{
      readCount[i] = 1;
    }
    console.log("true/false for initRate", nonUniformInitiationRatesStatus);
    if(nonUniformInitiationRatesStatus=='true'){
      readInitRate[i] = localStorage.getItem('initRateOf_'+str(i));
      console.log('readInitRate[i]', i, readInitRate[i]);
    }
    else{
      readInitRate[i] = 1;
    }
    console.log('readInitRate[i]', i, readInitRate[i]);
  }

  // initiate 2D array of the counts of terminated translation of each transcript at each time point:
  // initiate 2D array of the count of RPF number on each transcript at each time point
  for (iTimePoint=0; iTimePoint<36; iTimePoint++){
    protCountTimePointTranscript[iTimePoint] = [];
    RPFCountTimePointTranscript[iTimePoint] = [];
    for (iGeneID=0;iGeneID<iMax; iGeneID++){
      protCountTimePointTranscript[iTimePoint][iGeneID]=0;
      RPFCountTimePointTranscript[iTimePoint][iGeneID]=0;
    }
  }

    // Standard genetic code defined once and for all as a 3D array:
  // A, C, G, U maps 0, 1, 2, 3 indexes: GenetiCode[0][0][0] = K Lysine decoded by AAA...
  GenetiCode[0][0][0] = 'K'; // lysine Lys (AAA)
  GenetiCode[0][0][1] = 'N'; // asparagine Asn (AAC)
  GenetiCode[0][0][2] = 'K'; // lysine Lys (AAG)
  GenetiCode[0][0][3] = 'N'; // asparagine Asn (AAU)
  GenetiCode[0][1][0] = 'T'; // threonine Thr (ACA)
  GenetiCode[0][1][1] = 'T'; // threonine Thr (ACC)
  GenetiCode[0][1][2] = 'T'; // threonine Thr (ACG)
  GenetiCode[0][1][3] = 'T'; // threonine Thr (ACU)
  GenetiCode[0][2][0] = 'R'; // arginine Arg (AGA)
  GenetiCode[0][2][1] = 'S'; // serine Ser (AGC)
  GenetiCode[0][2][2] = 'R'; // arginine Arg (AGG)
  GenetiCode[0][2][3] = 'S'; // serine Ser (AGU)
  GenetiCode[0][3][0] = 'I'; // isoleucine Ile (AUA)
  GenetiCode[0][3][1] = 'I'; // isoleucine Ile (AUC)
  GenetiCode[0][3][2] = 'M'; // methionine Met (AUG) --> START
  GenetiCode[0][3][3] = 'I'; // isoleucine Ile (AUU)
  GenetiCode[1][0][0] = 'Q'; // glutamine Gln (CAA)
  GenetiCode[1][0][1] = 'H'; // histidine His (CAC)
  GenetiCode[1][0][2] = 'Q'; // glutamine Gln (CAG)
  GenetiCode[1][0][3] = 'H'; // histidine His (CAU)
  GenetiCode[1][1][0] = 'P'; // proline Pro (CCA)
  GenetiCode[1][1][1] = 'P'; // proline Pro (CCC)
  GenetiCode[1][1][2] = 'P'; // proline Pro (CCG)
  GenetiCode[1][1][3] = 'P'; // proline Pro (CCU)
  GenetiCode[1][2][0] = 'R'; // arginine Arg (CGA)
  GenetiCode[1][2][1] = 'R'; // arginine Arg (CGC)
  GenetiCode[1][2][2] = 'R'; // arginine Arg (CGG)
  GenetiCode[1][2][3] = 'R'; // arginine Arg (CGU)
  GenetiCode[1][3][0] = 'L'; // leucine Leu (CUA)
  GenetiCode[1][3][1] = 'L'; // leucine Leu (CUC)
  GenetiCode[1][3][2] = 'L'; // leucine Leu (CUG)
  GenetiCode[1][3][3] = 'L'; // leucine Leu (CUU)
  GenetiCode[2][0][0] = 'E'; // glutamate Glu (GAA)
  GenetiCode[2][0][1] = 'D'; // aspartate Asp (GAC)
  GenetiCode[2][0][2] = 'E'; // glutamate Glu (GAG)
  GenetiCode[2][0][3] = 'D'; // aspartate Asp (GAU)
  GenetiCode[2][1][0] = 'A'; // alanine Ala (GCA)
  GenetiCode[2][1][1] = 'A'; // alanine Ala (GCC)
  GenetiCode[2][1][2] = 'A'; // alanine Ala (GCG)
  GenetiCode[2][1][3] = 'A'; // alanine Ala (GCU)
  GenetiCode[2][2][0] = 'G'; // glycine Gly (GGA)
  GenetiCode[2][2][1] = 'G'; // glycine Gly (GGC)
  GenetiCode[2][2][2] = 'G'; // glycine Gly (GGG)
  GenetiCode[2][2][3] = 'G'; // glycine Gly (GGU)
  GenetiCode[2][3][0] = 'V'; // valine Val (GUA)
  GenetiCode[2][3][1] = 'V'; // valine Val (GUC)
  GenetiCode[2][3][2] = 'V'; // valine Val (GUG)
  GenetiCode[2][3][3] = 'V'; // valine Val (GUU)
  GenetiCode[3][0][0] = 'stop'; // STOP (UAA) --> STOP
  GenetiCode[3][0][1] = 'Y'; // thyrosine Tyr (UAC)
  GenetiCode[3][0][2] = 'stop'; // STOP (UAG) --> STOP
  GenetiCode[3][0][3] = 'Y'; // thyrosine Tyr (UAU)
  GenetiCode[3][1][0] = 'S'; // serine Ser (UCA)
  GenetiCode[3][1][1] = 'S'; // serine Ser (UCC)
  GenetiCode[3][1][2] = 'S'; // serine Ser (UCG)
  GenetiCode[3][1][3] = 'S'; // serine Ser (UCU)
  GenetiCode[3][2][0] = 'stop'; // STOP (UGA) --> STOP
  GenetiCode[3][2][1] = 'C'; // cysteine Cys (UGC)
  GenetiCode[3][2][2] = 'W'; // Tryptophan Trp (UGG)
  GenetiCode[3][2][3] = 'C'; // cysteine Cys (UGU)
  GenetiCode[3][3][0] = 'L'; // leucine Leu (UUA)
  GenetiCode[3][3][1] = 'F'; // phenylalanine Phe (UUC)
  GenetiCode[3][3][2] = 'L'; // leucine Leu (UUG)
  GenetiCode[3][3][3] = 'F'; // phenylalanine Phe (UUU)

  // lattices instantiations (create lattices objects):
  // The proline coding codons (4 of them) will be highlighted with a differnt color.
  // To do this, you need to spot all these P coding codons in the transcript to be instantiated and displayed.
  // The + and - charged aa coding codons will also be highlighted in different colors.
  // instantiate the transcripts:
  var interspersedDist = 65; //90
  // create and initiate the transcripts by reading them from the local storage retrieved list 'StoredmRNA'.
  // We take all the transcripts that were dropped in the settings html page.
  transcriptCardinality = StoredmRNA.length;

  // initiate 'riboRecruitScore' to generate a multinomial random sampling of
  // the transcripts for recruitment of the ribosomes for transcript initiation
  // relative probabilities (not uniform initiation probability):
  indicesListToSampleFrom = [];
  // is the list of transcript id number to sample from. In this list, the same
  // transcript id number will appear proportionnaly to the number of reads and
  // proportionnaly to the number of recruitment score for initiation.
  for (i_tr=0; i_tr<transcriptCardinality; i_tr++){
    console.log('+++++++++ -- ++ readInitRate[i]**************  :', i_tr, readInitRate[i_tr]);
    riboRecruitScore[i_tr] = int(10*readInitRate[i_tr]);
    // it is scaled up by a factor 10 to take into account the first
    // decimal digit of the initiation rates.
    for (i_read=0; i_read<readCount[i_tr]; i_read++){
      var read_str = str(StoredGeneID[i_tr])+'#'+str(i_read);
      for (i_score = 0; i_score < riboRecruitScore[i_tr]; i_score++){
        indicesListToSampleFrom.push(indexOfRead(StoredGeneID, readCount, read_str));
        // build indicesSample
      }
    }
  }// end of outer for loop.
  console.log('riboRecruitScore', riboRecruitScore);
  console.log('indicesListToSampleFrom.length', indicesListToSampleFrom.length);
  console.log('indicesListToSampleFrom[0]', indicesListToSampleFrom[0]);
  console.log('indicesListToSampleFrom[24]', indicesListToSampleFrom[24]);
  console.log('indicesListToSampleFrom[25]', indicesListToSampleFrom[25]);


  read_added_count = 0;
  for (i_tr=0; i_tr<transcriptCardinality; i_tr++){
    //arguments of Lattice: (tempStartX, tempStartY, tempLength, tempSequence, tempCopy, tempName, tempNumID, tempProList, tempProStatus, tempAAplusList, tempAAminusList, tempTunnelStatus)
    //lattice.push(new Lattice(130, 25, StoredmRNA[i].length/3, StoredmRNA[i], 1, StoredGeneID[i], 0));
    // Determine the positions of proline codons:
    listOfProlinePositions = [];
    for (j=0; j<listProlinePosition(StoredmRNA[i_tr]).length; j++){
      listOfProlinePositions.push(listProlinePosition(StoredmRNA[i_tr])[j]);
    }
    // Determine AAplus codons positions:
    listOfAAplusPositions = [];
    for (j=0; j<listAAplusPosition(StoredmRNA[i_tr]).length; j++){
      listOfAAplusPositions.push(listAAplusPosition(StoredmRNA[i_tr])[j]);
    }
    // Determine AAminus codons positions:
    listOfAAminusPositions = [];
    for (j=0; j<listAAminusPosition(StoredmRNA[i_tr]).length; j++){
      listOfAAminusPositions.push(listAAminusPosition(StoredmRNA[i_tr])[j]);
    }
    lattice.push(new Lattice(130, 50 + interspersedDist * i_tr, StoredmRNA[i_tr].length/3, StoredmRNA[i_tr], readCount[i_tr], StoredGeneID[i_tr], 0, listOfProlinePositions, prolineStatus, listOfAAplusPositions, listOfAAminusPositions, tunnelStatus));
    // Below we add the transcripts copies (if there are copies).
    //The y position of the copies starts below all the different transcripts names.
    // We do not add the copy #0 that was already pushed just before.
    if (readCount[i_tr]>1){
      for (i_copy=0; i_copy <readCount[i_tr]-1; i_copy++){
        lattice.push(new Lattice(130, 50 + interspersedDist * (transcriptCardinality + read_added_count + i_copy), StoredmRNA[i_tr].length/3, StoredmRNA[i_tr], readCount[i_tr], StoredGeneID[i_tr], i_copy+1, listOfProlinePositions, prolineStatus, listOfAAplusPositions, listOfAAminusPositions, tunnelStatus));
      }
      read_added_count += (readCount[i_tr]-1);
    }
  }
  // compute the full transcripts cardinality:
  transcriptFullCardinality = lattice.length;

  // compute the total length of the instantiated transcriptome:
  for (i=0; i<transcriptFullCardinality; i++){
    transcriptomeLength_nts += lattice[i].mRNAlength;
  }
  transcriptomeLength_codons = transcriptomeLength_nts/3;
  // check:
  console.log('total number of transcripts reads copies=', lattice.length);
  //for (i=0; i<lattice.length; i++){
  //  console.log('#', i+1, lattice[i].name, 'copy#=',lattice[i].numID, 'over', lattice[i].copyIndex);
  //}
  console.log('total length of sampling size=', indicesListToSampleFrom.length);
  console.log('counts of transcript id number references in sampling list:');
  for (i=0; i<lattice.length; i++){
    var ref_count = 0;
    for (i_sample=0; i_sample<indicesListToSampleFrom.length; i_sample++){
      if(indicesListToSampleFrom[i_sample] == i){
        ref_count +=1;
      }
    }
    console.log("::::::::::@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@::::::::::::::::::");
    console.log('#', i+1, lattice[i].name, 'copy#=',lattice[i].numID, 'over', lattice[i].copyIndex, '#samples=',ref_count);
  }
  // initialize a 'pseudomatrix' that will record the cumulated translocation
  // timer set points for each codon of each particular transcript read copy
  // that were exactly retained during the protein elongation simulation.
  // Initialize a second 'pseudomatrix'that will record the real ribosome dwell time
  // that were observed on each codon of each transcript read copy.
  // And again with a third pseudo matrix with just the ribosome
  // occurence footprint count upon sampling time point.
  for (tr=0; tr<transcriptFullCardinality; tr++){
    RiboSeqSPtruth[tr] = []; // initialize second dimension of array.
    RiboSeqRRTtruth[tr] = [];
    RiboSeqFCsampled[tr] = []
    for (i_codon=0; i_codon<lattice[tr].mRNAlength/3; i_codon++){
      // note that each row may have a different number of columns!
      // This beacuse each transcript read may have a different length.
      RiboSeqSPtruth[tr].push(0);
      RiboSeqRRTtruth[tr].push(0);
      RiboSeqFCsampled[tr].push(0);
    }
  }
  //check
  console.log("*********** R I B O - S E Q  check ****************");
  //console.log('index of (RPL4)#0=', indexOfRead(StoredGeneID, readCount, "(RPL4)#0"));
  var i_RPL4 = indexOfRead(StoredGeneID, readCount, "(RPL4)#0");
  //console.log("SP truth RPL4", RiboSeqSPtruth[i_RPL4]);
  //console.log("RRT truth RPL4", RiboSeqRRTtruth[i_RPL4]);
  //console.log("FC sampled RPL4", RiboSeqFCsampled[i_RPL4]);


  // Ribosomes instantiations (create ribosomes objects, see constructor and methods in 'ObjectsClasses.js').
  // (free) ribosomes instantiations:
  const riboInterY = 80;
  const riboInterX = 60;
  const riboPerLine = 20;
  freeRibosomeInitialNumber = floor(riboPoolRatio * transcriptFullCardinality);
  var numbRiboLines = int(freeRibosomeInitialNumber/riboPerLine) + 1;
  for (nbFree=0; nbFree < freeRibosomeInitialNumber; nbFree++){
    i_line = int(nbFree/riboPerLine);
    ribosome.push(new Ribo(80+(nbFree % riboPerLine)*riboInterX, 1200 + i_line * riboInterY, chosenSpecies));
  }
  //console.log('number of ribosomes instance=', ribosome.length);
}  //end setup function
//*****************************************************************
// end of function setup
//*****************************************************************
function indexOfRead(geneIDlist, readCopyList, readIDstring){
  //This function returns the index (integer) of the transcript read list (tr index)
  // corresponding to the geneIDlist given in first argument, the readCopyList given
  // in second argument and the unique read ID string given as the third argument (uniqueID to query for).
  var itemList = [];
  for (i=0; i<geneIDlist.length; i++){
    for (j=0; j<readCopyList[i]; j++){
      itemList.push(geneIDlist[i]+"#"+str(j));
    }
  }
  //console.log('list built in indexOfread function:', itemList);
  return itemList.indexOf(readIDstring);
}

function decode(codon){
  // This function receives a codon in argument as a 3 letter string.
  // This function returns the corresponding decoded amino acid
  // (called 'residue') using the standard genetic code.
  var i, j, k;
  var FirstLetter, SecondLetter, ThirdLetter;
  var alphabet = "ACGT"; // "ACGU"
  var residue;
  FirstLetter = codon.substring(0, 1);
  SecondLetter = codon.substring(1, 2);
  ThirdLetter = codon.substring(2, 3);
  i = alphabet.indexOf(FirstLetter);
  j = alphabet.indexOf(SecondLetter);
  k = alphabet.indexOf(ThirdLetter);
  residue = GenetiCode[i][j][k];
  return residue;
}

function CDStoRNA(tCodon){
  // this function converts a codon in DNA nucleotides to a codon in mRNA nucleotides.
  // It replaces the Ts by Us.
  var firstLetter, secondLetter, thirdLetter;
  var mRNACodon;
  var CDScodon = str(tCodon);
  firstLetter = CDScodon.substring(0, 1);
  secondLetter = CDScodon.substring(1, 2);
  thirdLetter = CDScodon.substring(2, 3);
  if (firstLetter == 'T'){
    firstLetter = 'U';
  }
  if (secondLetter == 'T'){
    secondLetter = 'U';
  }
  if (thirdLetter == 'T'){
    thirdLetter = 'U';
  }
  mRNACodon = str(firstLetter)+str(secondLetter)+str(thirdLetter);
  return mRNACodon;
}

function highlight(){
  dropzone.style('background-color', '#C12C6');//#787878 or rgb(205, 160, 135, 0.5)
}
function unhighlight(){
  dropzone.style('background-color', '#BC4F4F');
}
function listProlinePosition(seq){
  // this function returns a list of the positions of codons coding for proline
  // from its sequence mRNA string given in argument
  var listProline = [];
  var this_codon;
  var seqLength = seq.length;
  var codNumb = seqLength/3;
  for (i=0; i < codNumb; i++){
    this_codon = seq.substring(0+3*i, 3+3*i);
    if (decode(this_codon) == 'P'){
      listProline.push(i);
    }
  }
  return listProline;
} // end of the function
function listAAplusPosition(seq){
  // this function returns a list of the positions of codons coding for positively charged aa (R and K)
  // from its sequence mRNA string given in argument
  var listAAplus = [];
  var this_codon;
  var seqLength = seq.length;
  var codNumb = seqLength/3;
  for (i=0; i < codNumb; i++){
    this_codon = seq.substring(0+3*i, 3+3*i);
    if (decode(this_codon) == 'K' || decode(this_codon) == 'R'){
      listAAplus.push(i);
    }
  }
  return listAAplus;
} // end the function
function listAAminusPosition(seq){
  // this function returns a list of the positions of codons coding for positively charged aa (R and K)
  // from its sequence mRNA string given in argument
  var listAAminus = [];
  var this_codon;
  var seqLength = seq.length;
  var codNumb = seqLength/3;
  for (i=0; i < codNumb; i++){
    this_codon = seq.substring(0+3*i, 3+3*i);
    if (decode(this_codon) == 'D' || decode(this_codon) == 'E'){
      listAAminus.push(i);
    }
  }
  return listAAminus;
} // end the function

function MaxwellBoltzmannTunnelFactor(upstreamSeq){
  // this function uses the upstream window sequence of a codon of a given transcript to compute the Maxwell Boltzmann factor contributed by the tunnel
  // electrostatic interaction with the mobile window of the last max 50 aa upstream of P-site:
  // The unique input to compute this factor is the codon.usptreamWindow itself and the model published
  // in Physical Review E (Joiret et al.)
  // 1) We convert the upstreamSeq of nucleotides into an ordered sequence of 0, +1 and -1 for neutral, + and -
  // amino acids decoded from their triplets in the direction 3' to 5' !!!
  // Reverse direction, starting from carboxy-terminal P site aa to N-terminal amino acid at the distal part of the tunnel
  // When you pop the list, you first get the C-terminal end of the nascent chain (ListSign is a LIFO stack list).
  var listSign = [];
  var seqLength = upstreamSeq.length;
  var codNumb = seqLength/3;
  for (i=0; i < codNumb; i++){
    this_codon = upstreamSeq.substring(i*3, i*3+3);
    if (decode(this_codon) == 'D' || decode(this_codon) == 'E'){
      listSign.push(-1.0); // negatively charged aa
      }
    else {
      if (decode(this_codon) == 'K' || decode(this_codon) == 'R'){
        listSign.push(+1.0); // positively charged aa
        }
      else {
        listSign.push(0.0); // neutral aa
        }
      }
    }
  console.log('listCharge=', listSign);
  // 2) The list of axial forces mapping the positions of all aa in window of max 50 codons was
  // retrieved from JSON file (preloaded) to 'tunnelElectrostaticsAF' dictionary.
  // 2bis) If K or R (resp D or E) are at position 1 (Psite), the AF is due PTC electrostatics and
  // is much stronger (= -21.2 pN) than the contribution of tunnel (reference: Joiret et al. in CSBJ, 2023):
  // 3) The elementwise dot product of the two lists gives the elementwise list of axial forces
  // acting on the peptidyl-tRNA at P-site.
  // 4) The sum of these elementwise forces is the total force. You multiply with axial displacement
  // to get the mechanical work and for a given temperature, you get the MB factor affecting
  // the rate of step 2:
  var windowLength = listSign.length;
  var cumulatedAF = 0.0;
  var aaSign = 0.0;
  var afPTC = -21.2;
  for (i=0; i < windowLength; i++){
    aaSign = listSign.pop();
    if (aaSign != 0.0 && i==0){
      cumulatedAF += aaSign * afPTC;
    }
    cumulatedAF += aaSign * tunnelElectrostaticsAF[str(i+1)] // convert number integer value to string
  }
  var mechaWork;
  mechaWork = -1.0 * cumulatedAF * 0.25; // pN.nm
  // displacement parallel to axial forces : positive work!
  var kBT = 4.282; // pN.nm (at 310.15 K or 37 Celsius degree)
  var maxwBoltzFactor = exp(mechaWork/kBT); // exponential of Work/kBT
  return maxwBoltzFactor;
}// end maxwellBoltzmannFactor function

function HypoExpRandomTime(Pcodon, Acodon, species, upSeq){
  // This function returns a stochastic queueuing time for the elongation cycle
  // to be imposed on a ribosome when located at the P-site codon position.
  // This stochastic time is hypo-exponentially distributed. The sampling in the hypo-exponential
  // is equivalently built from 3 independent sampling from 3 exponential distributions with
  // pairwise differents rates describing the 3 elongation substeps.
  // The relevant rates are extracted from the species dictionaries that were preloaded in the
  // species specific JSON files.
  // Furthermore, there might be a Maxwell-Boltzmann correcting factor on a substep rate to
  // take into account a context factor like 'electrostatic interaction in the tunnel' or
  // or 'A to I tRNA editing ADAT enzyme silencing' effect on TAPSLIVR targeted A34 to I34 tRNAs or 'proline'
  // effect or 'secondary structure' effect. This is the reason of all the conditions (if then else).
  // The argument Pcodon is the string of the codon at P-site. Acodon is the string of the codon at A-site.
  // The argument species is the string name of the expression vector species.
  // The argument upSeq is the string sequence of nucleotides upstream of the P-codon.
  var expTime1, expTime2, expTime3; // iid exponential sampled times
  var rate1, rate2, rate3; // inverse of previous times
  var hypoExpTime; // sum of expTime1, expTime2, expTime3
  if (Acodon == 'TAA' || Acodon == 'TAG' || Acodon == 'TGA'){
    Acodon = Pcodon;
  }
  if (species=='yeast'){
    // rate1 is determined by substep one (accommodation) of the codon at A-site:
    rate1 = yeastHypoExpParams[Acodon][0];
    // rate2 is determined by substep two (peptide bond formation) of the tRNA cognate codon at P-site.
    rate2 = yeastHypoExpParams[Pcodon][1];
    // rate 3 is determined by substep three (eviction translocation) of the tRNA cognate codon at P-site.
    rate3 = yeastHypoExpParams[Pcodon][2];
  } // end if
  if (species=='coli'){
      // rate1 is determined by substep one (accommodation) of the codon at A-site:
      rate1 = coliHypoExpParams[Acodon][0];
      // rate2 is determined by substep two (peptide bond formation) of the tRNA cognate codon at P-site.
      rate2 = coliHypoExpParams[Pcodon][1];
      // rate 3 is determined by substep three (eviction translocation) of the tRNA cognate codon at P-site.
      rate3 = coliHypoExpParams[Pcodon][2];
  } // end if
  // Proline factor and EF-P or eIF5A elongation factor depletion (affects rate2)
  if (prolineStatus == 'true' && decode(Pcodon) == 'P'){
    // The peptide bond formation rate is at least t_1/2 = 5.102 ms (E.coli, Rodnina).
    if (species == 'coli'){
      rate2 = 0.6931/5102; // ln2 / t_1/2
    }
    if (species == 'yeast'){
      rate2 = 0.5 * 0.6931 / 5102; // ln2 / t_1/2
    }
  }//end if proline status
  if (EFPdepleteStatus == 'true' && decode(Acodon) == 'P'){
    // there is a fold change of 90 (90 times slower) in the rate when proline is at A site and EFP is depleted
    rate2 = rate2 / 90.0;
  }//end if EFP elongation factor depletion status
  // Exit tunnel electrostatic interaction factor (affects rate2):
  if (tunnelStatus == 'true'){
    rate2 = rate2 * MaxwellBoltzmannTunnelFactor(upSeq);
    //console.log('Maxwell Boltzmann =', MaxwellBoltzmannTunnelFactor(upSeq));
  }// end if Exit tunnel electrostatics interaction with the nascent chain (max 50 aa embedded in tunnel)
  // U34 tRNA modification ELP3- URM1- status:
  if (tRNAmodificationStatus == 'true' && ['K', 'Q', 'E'].indexOf(decode(Acodon))!==-1){
    console.log('ELP3- URM1- U34 tRNA modification enzyme is depleted', tRNAmodificationStatus);
    rate1 = rate1 * uridine34SilencedDict[Acodon]; // the uridine34SilencedDict has keys in dna (cds) codons.
  }
  // A to I tRNA editing ADAT enzyme silencing status
  // test if decode(RNAtoCDS(Acodon)) in ['T', 'A', 'P', 'S', 'L', 'V', 'R']:
  if (adatSilencedStatus == 'true' && ['T', 'A', 'P', 'S', 'L', 'I', 'V', 'R'].indexOf(decode(Acodon))!==-1){
    rate1 = rate1 * adatSilencedDict[CDStoRNA(Acodon)]; // recall that the adat dict. has keys in mRNA codons.
    //console.log('ADAT enzyme of A to I tRNA editing is silenced', adatSilencedStatus);
    //console.log('mRNA codon ADAT sensitive:', CDStoRNA(Acodon), ' rate 1 fold change=', adatSilencedDict[CDStoRNA(Acodon)]);
  }
  console.log('ADAT enzyme of A to I tRNA editing is silenced', adatSilencedStatus);
  console.log('ELP3- URM1- U34 tRNA modification enzyme is depleted', tRNAmodificationStatus);
  // mRNA secondatry structure downstream context (affects rate3):
  //console.log('tau 2 =', 1.0/rate2, 'decode P codon is', decode(Pcodon) == 'P');
  expTime1 = jStat.exponential.sample(rate1);
  expTime2 = jStat.exponential.sample(rate2);
  expTime3 = jStat.exponential.sample(rate3);
  hypoExpTime = expTime1 + expTime2 + expTime3;
  return hypoExpTime;
} // end the function
function riboTypeCounts(riboList){
  // counts the ribosome types of ribosomes:
  // riboList is the list of the ribosome objects that are currently instantiated
  // and either of type 'free' or 'initiated' (= translocating).
  // This function returns the tally of the free ribosomes and the tally of the initiated Ribosomes
  var freeRibCount = 0;
  var initiatedRibCount = 0;
  for (i=0; i<riboList.length; i++){
    //console.log('type of ribosome=', riboList[i].type);
    if (riboList[i].type == 'free'){
      freeRibCount += 1;
    }
    if (riboList[i].type == 'initiated' || riboList[i].type == 'translating'){
      initiatedRibCount += 1;
    }
  }// end for loop
  //console.log(freeRibCount, initiatedRibCount);
  return [freeRibCount, initiatedRibCount];
} // end function riboTypeCounts
function riboCountOnTranscript(transcriptID){
  // This function returns the number of ribosomes that are
  // currently elongating on the particular read of the transcript copy given in argument.
  var riboCountOnTranscript = 0;
  for (i=0; i<ribosome.length; i++){
    if (ribosome[i].uniqueReadID == transcriptID.uniqueID){
      riboCountOnTranscript+=1;
      //ribosome[i].print();
    }
  }
  return riboCountOnTranscript;
}

function riboCountOn16(riboList){
  var riboCount16 = 0;
  for (i=0; i<riboList.length; i++){
    if (riboList[i].type != 'free' && riboList[i].y < 1050){
      riboCount16 +=1;
    }
  }
  return riboCount16;
}

/*function fileready(lines){
  txt = join(lines, '\n');
}*/
/*
function loadFile(){
  loadStrings("Data/CDSfasta.txt", fileLoaded);
}
function fileLoaded(data){
  txt = data;
}
*/
function fileSelected(file){
  createP(file.name + " " + file.size + " " + file.type);
  //createP(file.data);
  //console.log(file.data);
}

function draw() {
  // enable the activation/desactivation of buttons for RPF display:
  buttonRPF.mousePressed(displayRPFOnOff);
  // callback function (see displayRPF() fnction)

  // test the status of the buttons:
  displayRPFStatus = str(buttonRPF.style('background-color')=='rgb(0, 218, 0)')
  // this displayRPFStatus is either == 'true' (green) or 'false' (red)

  array_cycle_draw_time = [];
  array_cycle_draw_time.push(millis());
  // put drawing code here
  //console.log(txt);
  //output.html(txt);
  background(235); //224

  // tRNAs jittering around:
  tRNAiso1_K.jitter();
  tRNAiso1_K.show();
  tRNAiso2_K.jitter();
  tRNAiso2_K.show();
  for (let i=0; i<tRNAmolNumber ; i++){
    tRNAsample[i].jitter();
    tRNAsample[i].show();
  }


  // mRNA transcripts:
  // The length of the transcript is proportional to number of codons
  // 1 codon = 3 pixels
  strokeWeight(4);
  strokeCap(ROUND);
  stroke(255, 0, 0);
  line(130, 1090, 1030, 1090); // 300 codons //y=1430
  //stroke(255, 0, 0);
  //line(400, 512, 1150, 512); // 250 codons
  //stroke(0, 0, 255);
  //line(400, 768, 1840, 768); // 480 codons

  // display only the first 16 lattices having a different nameID(< transcriptCardinality)
  var i_line = 0;
  var i_tr = 0;
  while (i_line<16){
    //console.log(StoredGeneID[i_line], 'compared to ', lattice[i_tr].name);
    if (lattice[i_tr].name == StoredGeneID[i_line] && lattice[i_tr].numID == 0){
        //console.log(StoredGeneID[i_line], '=', lattice[i_tr].name);
        tr_displayed_set[i_line] = i_tr;
        lattice[i_tr].displays();
        i_line +=1;
      }
    i_tr +=1;
  }
  //console.log('tr_displayed_set=', tr_displayed_set);

  // ---- UNIFORM ELONGATION MOTION ALONG Transcript with 300 codons: ------------
  // motion condition:
  if (yeastSpeciesStatus=='true'){
    RiboElongationRate = 1.4;
  }
  if (coliSpeciesStatus=='true'){
    RiboElongationRate = 2.8;
  }
  if (x <= 1030){
    x += RiboElongationRate;
  }
  else {
    x = 130; //30
  }
  // ---- RIBOSOME OBJECT DESIGN: ------------
  // define line center as ribosome center
  RiboX = x;
  RiboY = 1090; //1430
  strokeWeight(5); // mRNA tunnel in ribosome
  stroke(24, 120, 54, 180);
  line(RiboX-15, RiboY, RiboX+15, RiboY);
  strokeWeight(1);
  stroke(24);
  line(RiboX-15, RiboY, RiboX+15, RiboY);
  strokeWeight(4);
  stroke(24, 120, 54, 180);
  noStroke();
  // ---------- Large 60S subunit: -----------
  fill(24, 120, 54, 180);
  ellipse(RiboX, RiboY-22, 50, 40);
  // ---------- Small 40S subunit: -----------
  rect(RiboX-7, RiboY+3, 14, 16);
  arc(RiboX-7, RiboY+11, 16, 16, HALF_PI, PI+HALF_PI);
  arc(RiboX+7, RiboY+11, 16, 16, PI+HALF_PI, HALF_PI);
  // -----------protein exit tunnel:----------
  stroke(205, 25, 25, 180);
  line(RiboX-20, RiboY-31, RiboX, RiboY-15);
  // tests:
  // text type parameters:
  textFont("Helvetica");
  textSize(12);
  textAlign(LEFT);
  strokeWeight(0);
  stroke(120);
  //fill(0);
  // show the (free) ribosomes:
  for (let rfree=0; rfree<freeRibosomeInitialNumber; rfree++){
    ribosome[rfree].show();
  }

  // A full Monte Carlo Step recapitulates like this:
  // Scan all ribosomes and scan all transcripts :
  for (let rib = 0; rib < freeRibosomeInitialNumber; rib +=1){
    var i_multinomial_random = int(random(indicesListToSampleFrom.length));
    var pickedTranscript = int(indicesListToSampleFrom[i_multinomial_random]);

    for (let tr = 0; tr < transcriptFullCardinality; tr += 1){//transcriptCardinality was changed to Full...
      //console.log('tr', tr, ' -> current footprint of this transcript:', lattice[tr].currentFootprint);
      currentTime = millis();
      //attempt new initiations on any randomly chosen transcript:
      // pick randomly among all transcripts
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //var pickedTranscript = int(random(transcriptFullCardinality));//transcriptCardinality was changed to Full...
      // pick randomly according to the provided initiation probabilities file across transcripts:

      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ribosome[rib].initiates(lattice[pickedTranscript]);
      // record the successive timer set points and dwell times of a given unique read:
      //if (lattice[tr].uniqueID == '(RPL4)#0'){
        //var tr_RPL4read0 = tr;
        //if (ribosome[rib].uniqueReadID == lattice[tr_RPL4read0].uniqueID){
          //rib_FOCUS = rib;
          //ribosome[rib_FOCUS].show();
        //}
      //}
      // show the ribosomes that are on the displayed transcripts only:
      // tr_displayed_set is the list of these transcripts:
      if (tr_displayed_set.includes(tr)){
        //ribosome[rib].show();
      }

      //attempt new translocations:
      currentTime = millis();
      //var CumulatedTimeOnThisCodon = 0;
      ribosome[rib].translocates(lattice[tr]);
      //ribosome[rib_FOCUS].showSpecifics();
      // show the ribosomes that are on the displayed transcripts only:
      if (tr_displayed_set.includes(tr)){
        //ribosome[rib].show();
      }
    }// end inner nested loop

  }// end outer loop

  // update for each transcript name the number of produced translated proteins
  // i.e. aggregate the read copies Translated for each transcript name:
  var translated_count_list = [];
  for (tr=0; tr<iMax; tr++){
    var cumul_translated = 0;
    for (i_read=0; i_read<lattice.length; i_read++){
      if (lattice[i_read].name == StoredGeneID[tr]){
        cumul_translated += lattice[i_read].readTranslated;
      }
    }
    translated_count_list.push(cumul_translated);
  }
  //
  for (tr=0; tr<iMax; tr++){
    for (i_read=0; i_read<lattice.length; i_read++){
      if (lattice[i_read].name == StoredGeneID[tr]){
        lattice[i_read].countTranslated = translated_count_list[tr];
      }
    }
  }
  // update the ribosome protected fragment number on each transcript:
  for (i_read=0; i_read<lattice.length; i_read++){
    lattice[i_read].rpfCount = riboCountOnTranscript(lattice[i_read]);
  }

  // display the current counts of 'free' and 'initiated' ribsosomes and total initial number of ribosomes:
  freetype = riboTypeCounts(ribosome)[0];
  initiatedtype = riboTypeCounts(ribosome)[1];
  //console.log('free ribosomes=', freetype, 'initiated ribosomes=', initiatedtype, 'tot=', freeRibosomeInitialNumber);

  // text type parameters:
  textFont("Helvetica");
  textSize(22);
  textAlign(LEFT);
  push();
  strokeWeight(0);
  stroke(120);//stroke(120)
  fill(200, 0, 0);
  if (coliSpeciesStatus=='true' && yeastSpeciesStatus=='true'){
    //console.log('There is confusion in the species to use as expression vector. Two species were selected simulatneously.Please select only one species in the settings tab.');
    text("Error: please select only one species in the settings tab.", 30, RiboY + 40);
  }
  text(chosenSpecies, 18, RiboY+5);
  //
  // display the current value of the clock:
  fill(0);
  stroke(235);
  textSize(16);//24
  textAlign(LEFT);
  textFont("Helvetica"); //Helvetica, Georgia, Courrier new ?
  push();
  var time_running = millis();
  time_sec = time_running/1000;
  var time_minutes = floor(time_running/(1000*60));
  var time_seconds = int((time_running/1000) % 60);
  var timer_text = str(nfs(time_minutes, 3, 0))+" min:"+str(nfs(time_seconds, 2, 0))+" sec."
  text("Time running:"+nfs(time_minutes, 3, 0)+" min:"+nfs(time_seconds, 2, 0)+" sec.", 7, 17);
  pop();
  // store the last Ribo-SeqRRT and ribo-SeqSP:
  storeItem('RiboSeqRRT', RiboSeqRRTtruth);// all actual ribosome residence time on codons
  storeItem('RiboSeqSP', RiboSeqSPtruth);// all translocation timer set points

  // store the elapsed time every 5 minutes for 3 hours (180 minutes):
  // eventFlag
  // store the # of translated transcript for all transcripts:
  var timeInterval = 3;
  for (i=0; i<36; i++){
    if (time_minutes >= (i+1)*timeInterval && eventFlag[i] == 'false'){
      eventFlag[i] = 'true';
      // store on local storage (not session storage for the stored inf be available in another browser's page)
      StoredTimes[i] = timer_text;
      storeItem('timePoint#'+str(i), StoredTimes[i]);
      //for (){// for all transcripts ids, store the protein counts at this time point
      for (i_tr=0; i_tr<iMax; i_tr++){
        var cumul_translated_count = 0;
        var cumul_RPF_count = 0;
        for (i_read=0; i_read<lattice.length; i_read++){
          if (lattice[i_read].name == StoredGeneID[i_tr]){
            cumul_translated_count += lattice[i_read].readTranslated;
            cumul_RPF_count += lattice[i_read].rpfCount;
          }
        }
        protCountTimePointTranscript[i][i_tr] = cumul_translated_count;
        RPFCountTimePointTranscript[i][i_tr] = cumul_RPF_count;
      }
      storeItem('protCountTimePoint#'+str(i), protCountTimePointTranscript[i]);
      storeItem('rpfCountTimePoint#'+str(i), RPFCountTimePointTranscript[i]);
    }
  }
  //
  // Accumulating ribosome footprint count at codon resolution from systematic Ribo-Seq snapshots
  // Each and every 10 seconds, update and sum the ribosome count and footprint count on the transcripts reads.
  // Do this until the end of simulation. At the end of simulation, save the results in local storage
  // for display in Ribo-Seq output tab.The aggregation per transcript id over the read copies
  // will be carried out in the ribo_seq_sketch.
  //
  for (i_snap= 0; i_snap < snapshots_nb; i_snap++){
    if (time_sec >= (i_snap+1)*time_between_snapshots && eventSnaps[i_snap] == 'false'){
      eventSnaps[i_snap] = 'true';
      console.log('time_sec', time_sec);
      for (i_tr=0; i_tr<iMax; i_tr++){
        if (lattice[i_tr].rpfCount > 0){
          //console.log('i', i_tr, lattice[i_tr].uniqueID, lattice[i_tr].rpfCount, lattice[i_tr].currentFootprint);
          footprintedFragmentsCount += lattice[i_tr].rpfCount;
          // update the footprint count on this read after this snapshot event:
          for (i_rpf=0; i_rpf<lattice[i_tr].rpfCount; i_rpf++){
            // cover all footprinted codons of this rpf:
            var lowerLeftCodon = int(lattice[i_tr].currentFootprint[i_rpf]/3 - 4);
            var upperRightCodon = int(lattice[i_tr].currentFootprint[i_rpf]/3 + 6);
            // constrained by length of transcript in codon number:
            var lowerBound = int(constrain(lowerLeftCodon, 0, lattice[i_tr].length));
            var upperBound = int(constrain(upperRightCodon, 0, lattice[i_tr].length));
            console.log('left', lowerLeftCodon, 'right', upperRightCodon);
            console.log('lowerBound', lowerBound, 'upperBound', upperBound);
            for (i_codon = lowerBound; i_codon < upperBound; i_codon ++){
              RiboSeqFCsampled[i_tr][i_codon] += 1;
            }
          }
        }
      }// end inner for
    }// end if
  }// end outer for

  // store the cumulated sampled footprint counts in local storage (+associated rpf statistics)
  if (eventSnaps[snapshots_nb - 1] == 'true'){
    storeItem('riboSeqSnapShotsNb', snapshots_nb);
    storeItem('footprintFragNb', footprintedFragmentsCount); // total number of RPF sampled at the end of simulation
    storeItem('riboSeqFC', RiboSeqFCsampled);
  }

  //riboFreeCount, riboTranslatingCount = riboTypeCounts(ribosome);
  //freetype, initiatedtype = riboTypeCounts(ribosome);
  push();
  freetype = riboTypeCounts(ribosome)[0];
  initiatedtype = riboTypeCounts(ribosome)[1];
  visibletype = riboCountOn16(ribosome);

  textFont("Helvetica"); //Courier or Helvetica or Calibri or Garamond or Verdana or Lato
  textSize(16);
  fill(0, 90, 90); //(75, 0, 130)
  text("Ribosome pool "+": "+str(nfs(freeRibosomeInitialNumber, 3, 0)), 15, RiboY+38);
  text("Transcripts count : "+str(transcriptFullCardinality), 15, RiboY+58);
  text("Current count of translating ribosomes "+": "+str(initiatedtype), 195, RiboY+38);
  //text("Current count of ribosomes on the 16 displayed transcripts "+": "+str(visibletype), 550, RiboY+58);
  text("Current count of free ribosomes "+": "+str(freetype), 195, RiboY+58);
  text("Average distance between two ribosomes on this transcriptome "+str(floor(transcriptomeLength_nts/initiatedtype))+
  " [nts] or "+str(floor(transcriptomeLength_nts/(3*initiatedtype)))+" [codons].", 550, RiboY+38);
  text("Beta = 1 / Initiation rate = "+str(nfs(float(InitBETA/1000), 3, 1)) +" [s].", 550, RiboY+58);
  text("Half-life of a free ribosome before initiation = "+str(nfs(0.6931472*InitBETA/1000, 3, 1))+" [s].", 825, RiboY+58);
  //var queriedName = "(COL1A1)";
  //text("nb ribo on "+queriedName+":"+riboCountOnTranscript(queriedName), 25, RiboY+22);
  //queriedName = "(RPL4)";
  //text("nb ribo on "+queriedName+":"+riboCountOnTranscript(queriedName), 25, RiboY+42);
  //queriedName = "(RPL22)";
  //text("nb ribo on "+queriedName+":"+riboCountOnTranscript(queriedName), 25, RiboY+62);
  //queriedName = "(EIF5A)";
  //text("nb ribo on "+queriedName+":"+riboCountOnTranscript(queriedName), 25, RiboY+82);
  pop();

  //var new_codonP = new Codon(lattice[0], 1461);
  //new_codonPstring = new_codonP.type;
  //console.log('codon 1461 on COL1A1=', new_codonPstring);

  transcriptLength = sequence.length;
  var triplet;
  for (var n = 0; n < transcriptLength/3; n += 1){
    triplet = sequence.substring(3*n, 3*n+3);
    peptide[n] = decode(triplet);
  }

  array_cycle_draw_time.push(millis());
  //console.log('cycle time (ms) = ', array_cycle_draw_time[array_cycle_draw_time.length -1] - array_cycle_draw_time[array_cycle_draw_time.length -2]);
  cycle_time_draw_loop = array_cycle_draw_time[array_cycle_draw_time.length -1] - array_cycle_draw_time[array_cycle_draw_time.length -2];
  //cycle_time_draw_loop = 0;
} // end function draw

// objects classes:
// class Trna:
