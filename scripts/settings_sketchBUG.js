//--javascript sketch for html page: settings
// Marc Joiret - last update: Dec 31st, 2019.
// Coding history:
// Settings page design and transcripts' file drop zone + data processing of
// codons - Dec 29th, 2019.
// Save relevant information in 4 files ~ table format - Dec 31st, 2019.
// 'myTranscriptsLog.tsv', 'myRNAids.tsv', 'myRNAseq.tsv', 'myRNAmkTable.tsv'
// These 4 files are stored in the /scripts library next to the sketch.js to be
// used by other .js sketches from other html pages.
//------------------------------------------------------------------------------
// variable section:
//------------------------------------------------------------------------------
// genetic code 3D-array: 4 by 4 by 4 3D array with 64 elements:
//var GenetiCode = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]];
var dropzone;
var dropzone4_transcripts_copyNumber;
var riboPoolZone;
var riboSlider;
var sliderSetting = 10.0;
var rateInput;
var timeRateSet;
var rateBox; // for half-time life of a free ribosome
var rateColor = 'rgb(255, 255, 255, 255)';
var rateSetting; // input rate box
var rateSettingValue = 180; //---InitBeta = beta = 1 / initiation rate in seconds (time kinetics constant for initiation rate).
var textHalfTime;
var riboText;
var riboPoolText;
var riboPoolRatioValue;
var riboColor = 'white';
var fastaparag;
var fastatxt ="empty string";
var readtxt = 'empty string';
var readTag = [];
var readCount = [];
var readInitRate = [];
var sanity_check = 'true';
var sanity_check_length = 'true';
var canvas_settings;
var textProcessed = "STATUS: waiting for transcripts...";
var mRNAs = [];
var myRNAs = [];
var lastCodon = [];
var ORFflag = [];
var mRNAseq = [];
var s;
var fastaSplit =[];
var myarray = [];
var output1="1", output2= "2",  output = "0";
var output11 = [];// = ["11"];
var output22 =[]; // = ["22"];
var numbLines;
var gID = []
var uniqArray= [];
var uniqString;
var nucleotideSeqArray = [];
var uniqReadArray = [];
var uniqReadString;
//var uniqStringBIS = ""; // empty string
var CDSstartIndex=0, CDSstopIndex=10;
var CDSmRNA;
var upperLim;
var idTagsMaxNumb = 16;  // max number of transcripts that will be displayed on the settngs page by Ribosomer.
var outFile;
let myTranscriptsTable;
var keyString11;
var keyString22;
var keyI;
var prolineButtonColor;
var buttonProline;
var buttonEFPdeplete;
var buttonExitTunnel;
var buttonSecondaryStructure;
var buttonNonUniformmRNAabundance;
var buttonNonUniformInitiationRates;
var buttontRNAmodification;
var buttontRNAmultipleModification;
var buttonLimitedPoolRibosomes;
var buttonColi;
var buttonYeast;
var buttonWorm;
var buttonZebra;
var buttonMouse;
var buttonHuman;
var ButtonDefaultColor = 'rgb(255, 0, 0)'; //red
var ButtonColorGreen = 'rgb(0, 218, 0)'; //green
var prolineStatus;
var EFPdepleteStatus;
var exitTunnelStatus;
var secondaryStructureStatus;
var nonUniformmRNAabundanceStatus;
var nonUniformInitiationRatesStatus;
var tRNAmodificationStatus;
var tRNAmultipleModificationStatus;
var limitedPoolRibosomeStatus;
var coliSpeciesStatus;
var yeastSpeciesStatus;
var wormSpeciesStatus;
var zebraSpeciesStatus;
var mouseSpeciesStatus;
var humanSpeciesStatus;
var input_tRNA_anticodon;
var input_tRNA_reciprocal_lambda1;
var input_tRNA_reciprocal_lambda2;
var input_tRNA_reciprocal_lambda2;
//var ComplementaryTriplet;
var textCodon;
var legitTriplet = [];
//var aa_residue;
//var aa_res;
var textResidue;
var tRNAmod_tau = [45.0, 1.0, 20.0]; // [tau1, tau2, tau3]
var t_ms;
var tauList = [91, 1.09, 6.82];
// compute the points (2 hypoexponential functions) to display
var tSup = 300; //ms //200 ??
var tInf = 0; //ms
var yInf = 0; //pdf
var ySup = 0.012;
var pxInf = 10;
var pxSup = 510;
var pyInf = 410;
var pySup = 740;
var pxInf2 = 545; // 610
var pxSup2 = 1045; //1110
var pyInf2 = 410;
var pySup2 = 740;
var x_tick = [];
var px_tick = [];
var px_tick2 = [];
var numbPointstoDraw; // integer
var dt_ms = tSup/numbPointstoDraw;
var t_2draw = [];
var pdf_mod = []; // modified pdf distribution for elongation cycle;
var pdf_ref = []; // reference pdf for the same tRNA anticodon;
var pdf_mod_max;
var pdf_ref_max;
var px_ref = [];
var py_ref = [];
var px_mod = [];
var py_mod = [];
var px_syn = [];
var py_syn = [];
var py_pro = [];
var pdf_pro = [];
//---------------------------------------------------------------------------
//var codon;
// genetic code 3D-array: 4 by 4 by 4 3D array with 64 elements:
var GenetiCode = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]];
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


//------------------------------------------------------------------------------
// functions section:
//------------------------------------------------------------------------------
// dictionary from JSON files with hypoexponetial parameters for all 61 sense codons:
var yeastHypoExpParams, coliHypoExpParams;
var chosenSpecies; // yeast or 'coli' or....
function preload() {
  //tRNA1pic = loadImage("http://localhost:8080/tRNAnormal.png");
  yeastHypoExpParams = loadJSON("http://localhost:8080/Data/dataJSONyeast.json");
  coliHypoExpParams = loadJSON("http://localhost:8080/Data/dataJSONcoli.json");
}
function highlight(){
  dropzone.style("background-color", "grey");
}
function unhighlight(){
  dropzone.style('background-color', '#C8C4');//#C8C4
}
function highlightBis(){
  dropzone4_transcripts_copyNumber.style("background-color", "grey");
}
function unhighlightBis(){
  dropzone4_transcripts_copyNumber.style('background-color', '#C8C4');//#C8C4
}
function highlightRiboPool(){
  riboPoolZone.style("background-color", 'rgb(50, 215, 50)');
  riboPoolText.style("color", "black");
  riboTranscriptRatio.style("color", "black");

}
function unhighlightRiboPool(){
  riboPoolZone.style('background-color', '#CD7D6B');//#C8C4
  riboPoolText.style("color", "white");
  riboTranscriptRatio.style("color", "white");
  riboSlider.value(10);
}

function prolineOnOff(){
  // color the button in green/red
  // test button.style value:
  if (buttonProline.style('background-color') == 'rgb(255, 0, 0)')
  {
    buttonProline.style('background-color', ButtonColorGreen);
  }
  else
  {
    buttonProline.style("background-color", ButtonDefaultColor);
  }
}
function efpDepleteOnOff(){
  // color the button in green/red
  // test button.style value:
  if (buttonEFPdeplete.style('background-color') == 'rgb(255, 0, 0)')
  {
    buttonEFPdeplete.style('background-color', ButtonColorGreen);
  }
  else
  {
    buttonEFPdeplete.style("background-color", ButtonDefaultColor);
  }
}
function ExitTunnelOnOff(){
  // color the button in green/red
  // test button.style value:
  if (buttonExitTunnel.style('background-color') == 'rgb(255, 0, 0)')
  {
    buttonExitTunnel.style('background-color', ButtonColorGreen);
  }
  else
  {
    buttonExitTunnel.style("background-color", ButtonDefaultColor);
  }
}
function SecondaryStructureOnOff(){
  // color the button in green/red
  // test button.style value:
  if (buttonSecondaryStructure.style('background-color') == 'rgb(255, 0, 0)')
  {
    buttonSecondaryStructure.style('background-color', ButtonColorGreen);
  }
  else
  {
    buttonSecondaryStructure.style("background-color", ButtonDefaultColor);
  }
}
function NonUniformmRNAabundanceOnOff(){
  // color the button in green/red
  // test button.style value:
  if (buttonNonUniformmRNAabundance.style('background-color') == 'rgb(255, 0, 0)')
  {
    buttonNonUniformmRNAabundance.style('background-color', ButtonColorGreen);
  }
  else
  {
    buttonNonUniformmRNAabundance.style("background-color", ButtonDefaultColor);
  }
}
function NonUniformInitiationRatesOnOff(){
  // color the button in green/red
  // test button.style value:
  if (buttonNonUniformInitiationRates.style('background-color') == 'rgb(255, 0, 0)')
  {
    buttonNonUniformInitiationRates.style('background-color', ButtonColorGreen);
  }
  else
  {
    buttonNonUniformInitiationRates.style("background-color", ButtonDefaultColor);
  }
}
function tRNAmodificationOnOff(){
  // color the button in green/red
  // test button.style value:
  if (buttontRNAmodification.style('background-color') == 'rgb(255, 0, 0)')
  {
    buttontRNAmodification.style('background-color', ButtonColorGreen);
  }
  else
  {
    buttontRNAmodification.style("background-color", ButtonDefaultColor);
  }
}
function tRNAmultipleModificationOnOff(){
  // color the button in green/red
  // test button.style value:
  if (buttontRNAmultipleModification.style('background-color') == 'rgb(255, 0, 0)')
  {
    buttontRNAmultipleModification.style('background-color', ButtonColorGreen);
  }
  else
  {
    buttontRNAmultipleModification.style("background-color", ButtonDefaultColor);
  }
}
function limitedPoolRibosomesOnOff(){
  // color the button in green/red
  // test button.style value:
  if (buttonLimitedPoolRibosomes.style('background-color') == 'rgb(255, 0, 0)')
  {
    buttonLimitedPoolRibosomes.style('background-color', ButtonColorGreen);
    highlightRiboPool();
  }
  else
  {
    buttonLimitedPoolRibosomes.style("background-color", ButtonDefaultColor);
    unhighlightRiboPool();
  }
}
function coliSpeciesOnOff(){
  // color the button in green/red
  // test button.style value:
  if (buttonColi.style('background-color') == 'rgb(255, 0, 0)')
  {
    buttonColi.style('background-color', ButtonColorGreen);
    buttonYeast.style("background-color", ButtonDefaultColor);
  }
  else
  {
    buttonColi.style("background-color", ButtonDefaultColor);
    buttonYeast.style("background-color", ButtonColorGreen);

  }
}
function yeastSpeciesOnOff(){
  // color the button in green/red
  // test button.style value:
  if (buttonYeast.style('background-color') == 'rgb(255, 0, 0)')
  {
    buttonYeast.style('background-color', ButtonColorGreen);
    buttonColi.style('background-color', ButtonDefaultColor);
  }
  else
  {
    buttonYeast.style("background-color", ButtonDefaultColor);
    buttonColi.style('background-color', ButtonColorGreen);
  }
}

function getFile(file){
  // handle image and text differently:
  if (file.type != 'text'){
    fastaparag = createP("This file is wrong: a text file is required in FASTA format. ");
    fastaparag.style("fontSize", "20pt");
    fastaparag.style("family-font", "Lucida Console");
    fastaparag.style("textAlign", "center");
    fastaparag.style("color", "white");
    fastaparag.style("background-color", "red");
    fastaparag.position(10, 10, 'relative');
    fastaparag.parent('#canvas_settings_container'); //'#canvas_settings_container' or '#fastaHere'

  }
  else if (file.type == 'text') {
    // fastaparag = createP(file.data);
    // fastaparag.parent('#fastaHere');
    // if the file is a text, do something specific with its content:
    fastatxt = str(file.data);
    // Regex to have a unique single string without space, tab, carriage retun or end of line character:
    r_uniqString = /[^\n]/g; // regex to match everything that is not an end of line.
    uniqArray = fastatxt.match(r_uniqString);
    uniqString = join(uniqArray,""); // WARNING: this join function introduce a white space delimiter between each member of the array!
    // use string concatenation instead or a push function and loop through the array:
    // for (i=0; i<uniqArray.length; i++){
    //   uniqStringBIS += uniqArray[i];  // or equivalently uniqStringBIS.push(uniqArray[i])
    //   //There are no hidden white spaces here between the elements of the arrays being concatenated !!!
    // }

    // Regex to find the description field:
    r_descr =/\>.*\s?/g;
    // Regular expression to find the gene id tags:
    r_id=/\(\w+\)/g; // regular expression any one or more word or number between litteral parentheses(Flag is global).

    //output2 = r.test(fastatxt); // should be true as well
    output11 = fastatxt.match(r_descr); // there are arrays...
    output22 = fastatxt.match(r_id);
    //print('output22 length:', output22.length);
    //print('output22:', output22);
    //print('output11:', output11);
    //console.log('output22 length:', output22.length);

    result = createP('The number of retrieved gene id tags is ' + str(output22.length));
    // *******: try to display this § in canvas instead of #fastaHere...
    result.parent('#canvas_settings_container'); // '#canvas_settings_container' or '#fastaHere'
    // the following paragraph does not fit well in the html page. Check for css design !
    //createP('The number of retrieved gene id tags is ' + str(output22.length));

    //for (i=0; i<output22.length;i++){
    //  createP("Gene id tag retrieved:" + output22[i] + "- description: " + output11[i]);
    //}

    // chop the cds from the single uniqString:  possibly polluted with hidden white spaces CASE***
    CDSstartIndex = output11[0].length// length of first description field
    CDSstopIndex = uniqString.indexOf(output11[1]);
    CDSmRNA = uniqString.substring(CDSstartIndex, CDSstopIndex);
    CDSmRNAlength = CDSstopIndex-CDSstartIndex;

    // regex with a capturing group to match all mRNAs sequences:
    // reg_mRNA = /\b(ATG[ATCG\s]*TAG|TGA|TAA|[ATCG]{3})\s*\>?/g;
    reg_mRNA = /\b(ATG[ATCG\s]*(TAG|TGA|TAA|[ATCG]{3}))\s*\>?/g;
    // captured groups: group 0: the full match, group 1: mRNA seq, group 2: last codon.
    r_ATCG = /[ATCG]/g; // will be used to remove polluting white spaces!
    // we have a global g flag here and will require EXEC to manage the captured groups...
    // The group intended to be captured as $1 is a seq starting with ATG and ending with a stop codon (or not)
    // before the next '>' or the last codon.
    nucleotideSeqArray = uniqString.match(reg_mRNA);
    // each single execution of exec provides the next occurence of the match until you get null:
    capturedGroups = reg_mRNA.exec(uniqString);
    while (capturedGroups != null){
      //console.log('captured groups with exec: ', capturedGroups);
      capturedmRNA= capturedGroups[1]; // possibly contains whate spece...
      lastCodon.push(capturedGroups[2]);
      //console.log('captured mRNA=', capturedmRNA);
      numNucleotides = capturedmRNA.length; // possibly contains white spaces...
      //console.log('length of mRNA:', numNucleotides);
      // remove the white spaces in each particular captured group 1:
      capturedmRNAcorrect = capturedmRNA.match(r_ATCG);
      //console.log('capturedmRNAcorrect is ', capturedmRNAcorrect);
      ORFmRNA = "";
      for (i=0; i<capturedmRNAcorrect.length; i++){
        ORFmRNA += capturedmRNAcorrect[i];
      }
      //console.log(ORFmRNA);
      //console.log('correct length:', ORFmRNA.length);
      myRNAs.push(ORFmRNA);
      // check ORF status: 'yes' if last codon = STOP and seq length multiple of 3 nt.
      // 'frameshift' if seq length not a multiple of 3 nt and last codon = STOP.
      // '?' if last codon not a codon STOP.
      r_STOP = /TAG|TGA|TAA/; // list of stop codons
      if (r_STOP.test(capturedGroups[2])) {
        // the last codon is a stop codon:
        if (ORFmRNA.length % 3 == 0){
          // divisible by 3 was true
          flag = 'yes';
        }
        else {
          flag = 'shifted';
        }
      }

      if (!r_STOP.test(capturedGroups[2])) {
        // the last codon is not a stop codon
        flag = '?';
      }
      //console.log('ORFflag to be pushed:', flag);
      ORFflag.push(flag);

      capturedGroups = reg_mRNA.exec(uniqString);
    }

    //console.log('nucleotide seq array: ', nucleotideSeqArray);
    //console.log('all correct mRNAS:', myRNAs);

    // ----------- save the relevant setting results:   ---------------------------
    // // - the transcripts log:
    myTranscriptsTable = new p5.Table();
    myTranscriptsTable.addColumn('geneID');
    myTranscriptsTable.addColumn('description');
    myTranscriptsTable.addColumn('orf');
    myTranscriptsTable.addColumn('lastCodon');
    myTranscriptsTable.addColumn('mRNAlength');
    // we limited the number of rows to 16 here, but for all processing purposes,
    // the full table should be created...
    upperLim = max(idTagsMaxNumb, output22.length);
    // you should first clear all local storage data (BEFORE the for loop of course...)
    clearStorage();
    console.log('The local storage was cleared.');
    console.log(localStorage);
    console.log('upperLim=', upperLim);
    storeItem("numbLines", upperLim);
    for (i=0; i<upperLim; i++){
      let newRow = myTranscriptsTable.addRow();
      newRow.setString('geneID', output22[i].substring(1, output22[i].length-1).toUpperCase());
      newRow.setString('description', output11[i].substring(1, output11[i].length));
      newRow.setString('orf', str(ORFflag[i]));
      newRow.setString('lastCodon', lastCodon[i]);
      newRow.setNum('mRNAlength', myRNAs[i].length);
      // save the data in local storage:
      keyString22 = str(output22[i].substring(0, output22[i].length).toUpperCase());
      //console.log('keyString22', keyString22);
      storeItem(str(i), keyString22);
      storeItem(keyString22, myRNAs[i]);
      // store for later displaying/retrieving from local storage:
      storeItem('description' + str(i),  output11[i]);
      //console.log('ORFflag[i] to store', i, str(ORFflag[i]));
      storeItem('#orf#' + str(i), str(ORFflag[i]));
      storeItem('lastcodon' + str(i), lastCodon[i]);
    }
    // display in <id=fastaHere>:
  	//var i = 5;
  	//for (var i = 0; i < 14; i++) {
    // some statements
  	//}
    //const key1 = LocalStorage.getItem('1');
  	//document.getElementById("fastaHere").innerHTML = key1;
    //document.getElementById("fastaHere").innerHTML = i;

    saveTable(myTranscriptsTable, "myTranscriptsLog.tsv", "tsv");
    // This is saved in a file that will find itself in the downloads directory of the client,
    // not in datadase as service (server side), like FIREBASE (google).

    // Think to dump a file somewhere or in firebase via the session storage or local storage ??? or not

    // This alternative way of saving the table or a file enables the browser to prompt you for authorization...
    // // C:\Users\marcjoiret\Desktop\ONGOING\j5sketch RIBOSOMER\output
    //outFile = createWriter(C:\Users\marcjoiret\Dropbox\PC\Desktop\ONGOING\j5sketch RIBOSOMER\output\"myTranscripts.tsv", "tsv");
    outFile = createWriter("myTranscripts.tsv", "tsv");
    // // write the header line:
    outFile.print("geneID\tDescription\tORF\tlastCodon\tmRNAlength");
    upperLim = max(idTagsMaxNumb, output22.length);
    for (i=0; i<upperLim; i++){
      outFile.print(str(output22[i].substring(1, output22[i].length-1).toUpperCase()) + "\t" + str(output11[i].substring(1, output11[i].length))  + "\t" + str(ORFflag[i]) +"\t" + str(lastCodon[i]) +"\t" + str(myRNAs[i].length));
    }
    outFile.close();
    //-----------------------------------------------------------------------------

    // this number is not correct: there are other polluting characters in this string (possibly equals to the nb of lines...):
    // try and find out how many A, T, C, G there are and what are the hidden characters ?
    r_ATCG = /[ATCG]/g;
    nucleotidesSeq = CDSmRNA.match(r_ATCG);
    CDSmRNAlengthCorrect = nucleotidesSeq.length;
    //print('correct mRNA length', CDSmRNAlengthCorrect);
    //print('char at 69:', CDSmRNA.charAt(69),'char at 70:', CDSmRNA.charAt(70), 'char at 71:', CDSmRNA.charAt(71), 'char at 72:', CDSmRNA.charAt(72));
    //print('char at 139:', CDSmRNA.charAt(139),'char at 140:', CDSmRNA.charAt(140), 'char at 141:', CDSmRNA.charAt(141), 'char at 142:', CDSmRNA.charAt(142));
    // --> hidden characters are white spaces indeed !

    // chop the cds from the uniqStringBIS that is UNpolluted with white spaces   CASE*** CORRECT ***
    // not done here...

    }

    fastaSplit = splitTokens(fastatxt, ">");
    for (i=0; i<fastaSplit.length; i++){
      if (fastaSplit[i].length >= 50){
        mRNAs.push(str(fastaSplit[i]));
      }
    }

    count = 0;

    textProcessed = "STATUS: Your first "+ str(output22.length) + " transcripts have been processed and saved in a locally stored table to be used in the Ribosomer bench:";
}
//
function getRateFromInput(){
  rateSettingValue = rateSetting.value();
  timeRateSet.html(str(nfs(rateSettingValue, 3, 1)));
  timeRateSet.style("color", "white");
  timeRateSet.style('background-color', '#CD7D6B');
}

// function getRiboRatio(): returns the riboRatioText and save the ratio value in the appropriate variable
function getRiboRatioFromSlider(){
  // get the value from the ribo slider:
  riboPoolRatioValue = riboSlider.value();
  sliderSetting = riboPoolRatioValue;
  riboColor = 'blue';
  // store the value in local storage:
  //localStorage.storeItem("riboPoolRatio", riboPoolRatioValue);
  // return the value as text to be displayed by html method later on...
  var riboText = str(nfs(riboPoolRatioValue, 2, 1));
  //console.log('damned riboPoolRatioValue & riboText in the get function', riboPoolRatioValue, riboText);
  return riboPoolRatioValue;
}
// function getReads: get from dropped file the read counts of gene ID tags and optionnaly the initiation rate fold change
// This is where RNA-Seq data gets in:
function getReads(file){
  // handle image and text differently:
  if (file.type != 'text'){
    fastaparag = createP("This file is wrong: a text file is required in tab delimited format. ");
    fastaparag.style("fontSize", "20pt");
    fastaparag.style("family-font", "Lucida Console");
    fastaparag.style("textAlign", "center");
    fastaparag.style("color", "white");
    fastaparag.style("background-color", "red");
    fastaparag.parent('#canvas_settings_container'); //'#canvas_settings_container' or '#fastaHere'
    fastaparag.position(200, -30, 'relative');
  }
  else if (file.type == 'text') {
    // if the file is a text, do something specific with its content:
    readtxt = str(file.data);
    r_uniqString = /[^\n]/g; // regex to match everything that is not an end of line.
    uniqReadArray = readtxt.match(r_uniqString);
    uniqReadString = join(uniqReadArray,"");
    console.log('uniqReadString', uniqReadString);

    // Regular expression to find the gene id tags:
    read_id=/\(\w+\)/g; // regular expression any one or more word or number between litteral parentheses(Flag is global).

    readTag = readtxt.match(read_id);
    // there might be elements with lowercase. We want to convert in uppercase:
    for (i=0; i<readTag.length; i++){
      var upperString = readTag[i].toUpperCase();
      readTag[i] = upperString;
    }
    // Regex to find the read counts per transcript list and (optionnaly the initiation rates relative fold change):
    //(COL1A1)\t4\t1\r(RPL4)\t3\t1\r(RPL22)\t3\t1\r(EIF5A)\t1\t1\r(Bgn)\t1\t1\r(KIF4A)\t1\t1\r(KIF14)\t3\t1\r
    count_regex = /.+\t(\d+)\t(\d+)\s/g;
    //any character whatsoever followed by a frist group of any number of digits, followed by a tab followed by a
    // second group of any number of digits followed by an optional litteral dot followed by an optional digit,
    // followed by a space or a tab.

    // each single execution of exec provides the next occurence of the match until you get null:
    capturedReadGroups = count_regex.exec(uniqReadString);
    console.log('capturedReadGroups', capturedReadGroups);
    readCount = [];
    readInitRate = [];
    while (capturedReadGroups != null){
      readCount.push(float(capturedReadGroups[1]));
      console.log('read', capturedGroups[1]);
      readInitRate.push(float(capturedReadGroups[2]));
      console.log('init', capturedReadGroups[2]);
      capturedReadGroups = count_regex.exec(uniqReadString);
      console.log('capturedReadGroups', capturedReadGroups);
    }//end while loop

    if (output22.length == readTag.length){
      sanity_check_length = 'true';
      for (i=0; i<output22.length;i++){
        if (output22[i] != readTag[i]){
          sanity_check = 'false';
        }
      }
    }
    else{
      sanity_check_length = 'false';
    }

    if (sanity_check=='false' || sanity_check_length=='false'){
      // send error message
      fastaparag = createP("This file is wrong: the count reads gene tags do NOT match the transcripts tags. ");
      fastaparag.style("fontSize", "20pt");
      fastaparag.style("family-font", "Lucida Console");
      fastaparag.style("textAlign", "center");
      fastaparag.style("color", "white");
      fastaparag.style("background-color", "red");
      fastaparag.parent('#canvas_settings_container'); //'#canvas_settings_container' or '#fastaHere'
      fastaparag.position(200, -30, 'relative');
    }

    // save in local storage the transcripts copy number for all geneID:
    // idem for initiation rates...
    console.log('beeeeeennnnnnn here...@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@');
    for (i=0; i<readTag.length; i++){
      storeItem('readCopiesOf_'+str(i), readCount[i]);
      storeItem('initRateOf_'+str(i), readInitRate[i]);
    }// end for
  } // end else

}//end getReads function

// function threeLetterCode(): gives the amino acid in three letters code from single letter code:
function threeLetterCode(X){
  switch(X){
    case 'K': return 'Lys';
    case 'R': return 'Arg';
    case 'P': return 'Pro';
    case 'A': return 'Ala';
    case 'F': return 'Phe';
    case 'L': return 'Leu';
    case 'W': return 'Trp';
    case 'I': return 'Ile';
    case 'G': return 'Gly';
    case 'D': return 'Asp';
    case 'E': return 'Glu';
    case 'N': return 'Asn';
    case 'T': return 'Thr';
    case 'S': return 'Ser';
    case 'M': return 'Met';
    case 'C': return 'Cys';
    case 'Q': return 'Gln';
    case 'H': return 'His';
    case 'V': return 'Val';
    case 'Y': return 'Tyr';
    case 'stop': return 'STOP';
    default: return 'X';
  }
}
// function complementaryBase: finds the base complementary to the argument in RNA syntax.
function complementaryBase(N){
  switch(N){
    case 'A': return 'U';
    case 'U': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';
    default: return 'N';
  }
}
// display codon deduced from changed input anticodon:
function writeCodon(){
  var anticodonTriplet = input_tRNA_anticodon.value();
  //console.log('anticodontriplet', anticodonTriplet);
  var N1 = anticodonTriplet.substring(0, 1);
  //console.log('N1', N1);
  var N2 = anticodonTriplet.substring(1, 2);
  //console.log('N2', N2);
  var N3 = anticodonTriplet.substring(2, 3);
  //console.log('N3', N3);
  var tripletCodon = complementaryBase(N1) + complementaryBase(N2) + complementaryBase(N3);
  //console.log('triplet', tripletCodon);
  return tripletCodon;
}
function decode(){
  // This function receives an anticodon in argument as a 3 letter string.
  // This function returns the corresponding d amino acid
  // (called 'residue') using the standard genetic code.
  var anticodonTriplet = input_tRNA_anticodon.value();
  //console.log('anticodontriplet', anticodonTriplet);
  var N1 = anticodonTriplet.substring(0, 1);
  //console.log('N1', N1);
  var N2 = anticodonTriplet.substring(1, 2);
  //console.log('N2', N2);
  var N3 = anticodonTriplet.substring(2, 3);
  //console.log('N3', N3);
  var tripletCodon = complementaryBase(N1) + complementaryBase(N2) + complementaryBase(N3);
  var i, j, k;
  var FirstLetter, SecondLetter, ThirdLetter;
  var alphabet = "ACGU";
  var residue;
  FirstLetter = tripletCodon.substring(0, 1);
  SecondLetter = tripletCodon.substring(1, 2);
  ThirdLetter = tripletCodon.substring(2, 3);
  if (FirstLetter == null || FirstLetter == 'N' || SecondLetter == null || SecondLetter == 'N' || ThirdLetter == null || ThirdLetter == 'N'){
    residue = 'X';
  }
  else{
    i = alphabet.indexOf(FirstLetter);
    j = alphabet.indexOf(SecondLetter);
    k = alphabet.indexOf(ThirdLetter);
    residue = GenetiCode[i][j][k];
  }
  return residue;
}
function decodeCDS(triCDS){
  // This function receives a codon in argument as a 3 letter string.
  // This function returns the corresponding d amino acid
  // (called 'residue') using the standard genetic code.
  var tripletCodon = str(triCDS);
  var i, j, k;
  var FirstLetter, SecondLetter, ThirdLetter;
  var alphabet = "ACGT";
  var residue;
  FirstLetter = tripletCodon.substring(0, 1);
  SecondLetter = tripletCodon.substring(1, 2);
  ThirdLetter = tripletCodon.substring(2, 3);
  if (FirstLetter == null || FirstLetter == 'N' || SecondLetter == null || SecondLetter == 'N' || ThirdLetter == null || ThirdLetter == 'N'){
    residue = 'X';
  }
  else{
    i = alphabet.indexOf(FirstLetter);
    j = alphabet.indexOf(SecondLetter);
    k = alphabet.indexOf(ThirdLetter);
    residue = GenetiCode[i][j][k];
  }
  return residue;
}
function mRNAtoCDS(uCodon){
  // this function converts a codon in mRNA nucleotides to a codon in DNA nucleotides.
  // It replaces the Us by Ts.
  var firstLetter, secondLetter, thirdLetter;
  var cdsCodon;
  var mRNAcodon = str(uCodon);
  firstLetter = mRNAcodon.substring(0, 1);
  secondLetter = mRNAcodon.substring(1, 2);
  thirdLetter = mRNAcodon.substring(2, 3);
  if (firstLetter == 'U'){
    firstLetter = 'T';
  }
  if (secondLetter == 'U'){
    secondLetter = 'T';
  }
  if (thirdLetter == 'U'){
    thirdLetter = 'T';
  }
  cdsCodon = str(firstLetter)+str(secondLetter)+str(thirdLetter);
  return cdsCodon;
}
function findSynonymousList(resi){
  // this function returns the list of all synonymous codons for the
  // amino acid residue given in argument:
  // The function systematically scan the dictionary of the rates loaded from the JSON file:
  var synonymousList = [];
  //var trigram;
  var trigramList = Object.keys(yeastHypoExpParams);
  for (i=0; i<Object.keys(yeastHypoExpParams).length; i++){
    if (decodeCDS(trigramList[i]) == str(resi)){
      synonymousList.push(trigramList[i]);
    }
  }
  return synonymousList;
} // end function
function storeModifiedTaus(){
  if (input_tRNA_reciprocal_lambda1.value() != null){
    tRNAmod_tau[0] = float(input_tRNA_reciprocal_lambda1.value());
  }
  else{
    tRNAmod_tau[0] = tauList[0];
  }
  if (input_tRNA_reciprocal_lambda2.value() != null){
    tRNAmod_tau[1] = float(input_tRNA_reciprocal_lambda2.value());
  }
  else{
    tRNAmod_tau[1] = tauList[1];
  }
  if (input_tRNA_reciprocal_lambda3.value() != null){
    tRNAmod_tau[2] = float(input_tRNA_reciprocal_lambda3.value());
  }
  else{
    tRNAmod_tau[2] = tauList[2];
  }
}
function hypoexponential(t_ms, tauList){
  la1 = 1.0/tauList[0];
  la2 = 1.0/tauList[1];
  la3 = 1.0/tauList[2];
  var term1 = exp(-1.0*t_ms/tauList[0])/((la3 - la1)*(la2 - la1));
  var term2 = exp(-1.0*t_ms/tauList[1])/((la3 - la2)*(la1 - la2));
  var term3 = exp(-1.0*t_ms/tauList[2])/((la2 - la3)*(la1 - la3));
  return la1*la2*la3*(term1 + term2 + term3);
}

//------------------------------------------------------------------------------
// setup
//------------------------------------------------------------------------------
function setup(){
// noCanvas();
canvas_settings = createCanvas(1560, 2050); // 1320, 750
canvas_settings.parent("#canvas_settings_container");
frameRate(60);

dropzone = select('#dropzone4file');
dropzone.dragOver(highlight);
dropzone.dragLeave(unhighlight);
dropzone.drop(getFile, unhighlight);
//
// drop zone for transcripts \t copy number count file:
dropzone4_transcripts_copyNumber = select('#fastaHere');
dropzone4_transcripts_copyNumber.dragOver(highlightBis);
dropzone4_transcripts_copyNumber.dragLeave(unhighlightBis);
dropzone4_transcripts_copyNumber.drop(getReads, unhighlightBis);
// check if the geneID tags are matching exactly with the transcripts earlier privided:

// Create paragraph (for input box context) to set the average initiation rate of free ribosomes on transcripts
rateInput = createP("Initiation time rate in (s) = ");//rateSettingValue is the default value (converted in milliseconds);
rateInput.parent("myButtonContainer");
rateInput.position(1150, 50, 'relative');
rateInput.style('fontSize', '14px');
rateInput.style('textFont', 'Helvetica');
rateInput.style('color', rateColor); //rateColor

if (getItem("initiationRateTime") == null){
  rateSettingValue = 180;
}
else{
  rateSettingValue = getItem("initiationRateTime")/1000;
}

timeRateSet = createP(str(rateSettingValue));//rateSettingValue is the default value (converted in milliseconds);
timeRateSet.parent("myButtonContainer");
timeRateSet.position(1380, 50, 'relative');
timeRateSet.style('fontSize', '14px');
timeRateSet.style('textFont', 'Helvetica');
timeRateSet.style('color', 'white'); //rateColor
timeRateSet.style('background-color', '#B4B4B4');

// create the half-time rate displayed from above info
rateBox = createP("Initiation half-time in (s) = ");
rateBox.parent("myButtonContainer");
rateBox.position(1150, 68, 'relative');
rateBox.style('fontSize', '14px');
rateBox.style('textFont', 'Helvetica');
rateBox.style('color', rateColor); //rateColor

textHalfTime = createP("XXX.X");
textHalfTime.parent("#myButtonContainer");
textHalfTime.size(30);
textHalfTime.style('fontSize', "14px");
textHalfTime.style('textFont', "Helvetica");
textHalfTime.position(1325, 68, 'relative');

// Create Input box for the setting of the initiation rate:
rateSetting = createInput(str(rateSettingValue));
rateSetting.parent("#myButtonContainer");
//rateSetting.size(30);
rateSetting.position(1324, 63, 'relative');
rateSetting.style('width', '30px');
rateSetting.style('background-color', '#B4B4B4');
rateSetting.style('color', 'black');
rateSetting.input(getRateFromInput);

//riboPoolZone:
riboPoolZone = select("#RibosomePool");
//
if (getItem("riboPoolRatio") != null){
  //sliderSetting = getItem("riboPoolRatio");
  removeItem("riboPoolRatio");
  //console.log('riboPoolRatio as retrived from setup=', sliderSetting);
  //riboSlider.value(sliderSetting);
}
//if (limitedPoolRibosomeStatus == 'true'){
//  console.log('limited pool of ribosomes button is GREEN and status is ACTIVE.');
  // you need to change the default status at setup here:
  //limitedPoolRibosomesOnOff();
  //limitedPoolRibosomeStatus = 'false';
  //buttonLimitedPoolRibosomes.style("background-color", ButtonDefaultColor);
  //unhighlightRiboPool();
  //mousePressed(buttonLimitedPoolRibosomes);
//}

// Create slider for ratio of ribosomes cardinality over transcripts cardinality (included their copy number):
riboSlider = createSlider(0.1, 10.0, sliderSetting, 0.1);
riboSlider.parent("#RibosomePool");
riboSlider.position(870, 12);
riboSlider.style('width', '200px');
//
//getRiboRatioFromSlider();
// Create input to set ribosome over transcript ratio a given value
//riboInput = createInput(sliderSetting);
//riboInput.parent("#RibosomePool");
//riboInput.size(30);
//riboInput.position(950, 32);
//riboInput.changed(getRiboRatioFromInput, highlightRiboPool);
//
riboTranscriptRatio = createP("Ratio #Ribosomes / #Transcripts");
riboTranscriptRatio.parent('#RibosomePool');
riboTranscriptRatio.style('fontSize', "16px");
riboTranscriptRatio.style('textFont', "Helvetica");
riboTranscriptRatio.style('color', riboColor);
riboTranscriptRatio.position(1105, -5, 'relative');
//
riboPoolText = createP(sliderSetting);
riboPoolText.parent("#RibosomePool");
riboPoolText.size(30);
riboPoolText.style('fontSize', "18px");
riboPoolText.style('textFont', "Helvetica");
riboPoolText.style('color', riboColor);
riboPoolText.position(1105, 15, 'relative');
//
// create an input + cursors zone to set the tRNA modified data:
let text_tRNA = createP("Type the anticodon (3'-->5') of the tRNA for which you want to modify the rates :");
text_tRNA.parent("#myCanvas_tRNAmodification_container");
text_tRNA.style('fontSize', "16px");
text_tRNA.style('textFont', "Helvetica");
text_tRNA.position(10, 32, 'relative');

let codon_NNN = createP("5'-NNN-3' - codon");
codon_NNN.parent("#myCanvas_tRNAmodification_container");
codon_NNN.style('fontSize', "14px");
codon_NNN.style('textFont', "Helvetica");
codon_NNN.position(597, -5, 'relative');

let aa_residue = createP("amino acid residue");
aa_residue.parent("#myCanvas_tRNAmodification_container");
aa_residue.style('fontSize', "14px");
aa_residue.style('textFont', "Helvetica");
aa_residue.position(807, -5, 'relative');

textCodon = createP("NNN");
textCodon.parent("#myCanvas_tRNAmodification_container");
textCodon.size(30);
textCodon.style('fontSize', "14px");
textCodon.style('textFont', "Helvetica");
textCodon.position(612, 13, 'relative');

textResidue = createP("X - aaR");
textResidue.parent("#myCanvas_tRNAmodification_container");
textResidue.size(80);
textResidue.style('fontSize', "14px");
textResidue.style('textFont', "Helvetica");
textResidue.position(807, 13, 'relative');

let text_NNN = createP("3'-NNN-5' - anticodon (tRNA)");
text_NNN.parent("#myCanvas_tRNAmodification_container");
text_NNN.style('fontSize', "14px");
text_NNN.style('textFont', "Helvetica");
text_NNN.position(597, 57, 'relative');

input_tRNA_anticodon = createInput('NNN');
input_tRNA_anticodon.parent("#myCanvas_tRNAmodification_container");
input_tRNA_anticodon.size(30);
input_tRNA_anticodon.position(610, 50, 'relative');

input_tRNA_anticodon.changed(writeCodon, decode);
//input_tRNA_anticodon.mousePressed(decode);

let tRNA_lambda1 = createP("Step 1 exponential mean time (ms) :");
let tRNA_lambda2 = createP("Step 2 exponential mean time (ms) :");
let tRNA_lambda3 = createP("Step 3 exponential mean time (ms) :");
tRNA_lambda1.parent("myCanvas_tRNAmodification_container");
tRNA_lambda2.parent("myCanvas_tRNAmodification_container");
tRNA_lambda3.parent("myCanvas_tRNAmodification_container");
tRNA_lambda1.style('fontSize', "14px");
tRNA_lambda1.style('textFont', "Helvetica");
tRNA_lambda2.style('fontSize', "14px");
tRNA_lambda2.style('textFont', "Helvetica");
tRNA_lambda3.style('fontSize', "14px");
tRNA_lambda3.style('textFont', "Helvetica");

tRNA_lambda1.position(20, 85, 'relative');
tRNA_lambda2.position(20, 125, 'relative');
tRNA_lambda3.position(20, 165, 'relative');

input_tRNA_reciprocal_lambda1 = createInput();
input_tRNA_reciprocal_lambda2 = createInput();
input_tRNA_reciprocal_lambda3 = createInput();
input_tRNA_reciprocal_lambda1.parent("#myCanvas_tRNAmodification_container");
input_tRNA_reciprocal_lambda2.parent("#myCanvas_tRNAmodification_container");
input_tRNA_reciprocal_lambda3.parent("#myCanvas_tRNAmodification_container");
input_tRNA_reciprocal_lambda1.size(50);
input_tRNA_reciprocal_lambda2.size(50);
input_tRNA_reciprocal_lambda3.size(50);

input_tRNA_reciprocal_lambda1.position(260, 99, 'relative');
input_tRNA_reciprocal_lambda2.position(260, 139, 'relative');
input_tRNA_reciprocal_lambda3.position(260, 179, 'relative');

input_tRNA_reciprocal_lambda1.changed(storeModifiedTaus);
input_tRNA_reciprocal_lambda2.changed(storeModifiedTaus);
input_tRNA_reciprocal_lambda3.changed(storeModifiedTaus);

// descriptions of the 3 steps:
let tRNA_lambda1_desc = createP("accommodation & proofreading (&tau;_1)");
let tRNA_lambda2_desc = createP("peptide bond formation (&tau;_2)");
let tRNA_lambda3_desc = createP("translocation and E-site eviction (&tau;_3)");
tRNA_lambda1_desc.parent("myCanvas_tRNAmodification_container");
tRNA_lambda2_desc.parent("myCanvas_tRNAmodification_container");
tRNA_lambda3_desc.parent("myCanvas_tRNAmodification_container");
tRNA_lambda1_desc.style('fontSize', "14px");
tRNA_lambda1_desc.style('textFont', "Helvetica");
tRNA_lambda2_desc.style('fontSize', "14px");
tRNA_lambda2_desc.style('textFont', "Helvetica");
tRNA_lambda3_desc.style('fontSize', "14px");
tRNA_lambda3_desc.style('textFont', "Helvetica");

tRNA_lambda1_desc.position(340, 85, 'relative');
tRNA_lambda2_desc.position(340, 125, 'relative');
tRNA_lambda3_desc.position(340, 165, 'relative');

// draw the x-axis ticks (use map function):
for (i=0; i<6; i++){
  x_tick[i] = 50 + i * 50;
  px_tick[i] = map(x_tick[i], 0, 300, pxInf + 25, pxSup - 20);
  px_tick2[i] = map(x_tick[i], 0, 300, pxInf2 + 25, pxSup2 - 20);
}

// create a button to submit new values from above input:



// creation of buttons to activate the different factors separately
// PROLINE SLOW DOWN
buttonProline = createButton("Proline slowdown");
buttonProline.parent("myButtonContainer");
buttonProline.size(160,27);
buttonProline.style('font-size', '14px');
buttonProline.position(610, 35, 'relative');
// check button status in local storage: if not empty, recall the stored status:
if (getItem("prolineStatus") == null){
  buttonProline.style('background-color', ButtonDefaultColor);
}
if (getItem("prolineStatus") == 'true')
{
  buttonProline.style('background-color', 'rgb(0, 218, 0)');
}
if (getItem("prolineStatus") == 'false')
{
  buttonProline.style('background-color', 'rgb(255, 0, 0)');
}

// EF-P or eIF5A depletion
buttonEFPdeplete = createButton("EF-P or eIF5A depleted");
buttonEFPdeplete.size(170,27);
buttonEFPdeplete.position(610-6, 70, 'relative');
buttonEFPdeplete.parent("myButtonContainer");
buttonEFPdeplete.style('font-size', '14px');
// check button status in local storage: if not empty, recall the stored status:
if (getItem("EFPdepleteStatus") == null){
  buttonEFPdeplete.style('background-color', ButtonDefaultColor);
}
if (getItem("EFPdepleteStatus") == 'true')
{
  buttonEFPdeplete.style('background-color', 'rgb(0, 218, 0)');
}
if (getItem("EFPdepleteStatus") == 'false')
{
  buttonEFPdeplete.style('background-color', 'rgb(255, 0, 0)');
}

// Exit tunnel Electrostatic Interaction
buttonExitTunnel = createButton("Exit tunnel electrostatic interaction");
buttonExitTunnel.size(270, 27);
buttonExitTunnel.position(300, 70, 'relative');
buttonExitTunnel.parent("myButtonContainer");
buttonExitTunnel.style('font-size', '14px');
// check button status in local storage: if not empty, recall the stored status:
if (getItem("exitTunnelStatus") == null){
  buttonExitTunnel.style('background-color', ButtonDefaultColor);
}
if (getItem("exitTunnelStatus") == 'true')
{
  buttonExitTunnel.style('background-color', 'rgb(0, 218, 0)');
}
if (getItem("exitTunnelStatus") == 'false')
{
  buttonExitTunnel.style('background-color', 'rgb(255, 0, 0)');
}

// mRNA secondary Structure button
buttonSecondaryStructure = createButton("mRNA secondary structure");
buttonSecondaryStructure.size(210, 27);
buttonSecondaryStructure.position(330, 35, 'relative');
buttonSecondaryStructure.parent("myButtonContainer");
buttonSecondaryStructure.style('font-size', '14px');
// check button status in local storage: if not empty, recall the stored status:
if (getItem("secondaryStructureStatus") == null){
  buttonSecondaryStructure.style('background-color', ButtonDefaultColor);
}
if (getItem("secondaryStructureStatus") == 'true')
{
  buttonSecondaryStructure.style('background-color', 'rgb(0, 218, 0)');
}
if (getItem("secondaryStructureStatus") == 'false')
{
  buttonSecondaryStructure.style('background-color', 'rgb(255, 0, 0)');
}

// Non Uniform mRNA relative abundance button
buttonNonUniformmRNAabundance = createButton("Non-uniform mRNA relative abundance");
buttonNonUniformmRNAabundance.size(280, 27);
buttonNonUniformmRNAabundance.position(860, 35, 'relative');
buttonNonUniformmRNAabundance.parent("myButtonContainer");
buttonNonUniformmRNAabundance.style('font-size', '14px');
// check button status in local storage: if not empty, recall the stored status:
if (getItem("nonUniformmRNAabundanceStatus") == null){
  buttonNonUniformmRNAabundance.style('background-color', ButtonDefaultColor);
}
if (getItem("nonUniformmRNAabundanceStatus") == 'true')
{
  buttonNonUniformmRNAabundance.style('background-color', 'rgb(0, 218, 0)');
}
if (getItem("nonUniformmRNAabundanceStatus") == 'false')
{
  buttonNonUniformmRNAabundance.style('background-color', 'rgb(255, 0, 0)');
}

// Non uniform Initiation Rates button
buttonNonUniformInitiationRates = createButton("Non-uniform initiation rates");
buttonNonUniformInitiationRates.size(210, 27);
buttonNonUniformInitiationRates.position(900, 70, 'relative');
buttonNonUniformInitiationRates.parent("myButtonContainer");
buttonNonUniformInitiationRates.style('font-size', '14px');
// check button status in local storage: if not empty, recall the stored status:
if (getItem("nonUniformInitiationRatesStatus") == null){
  buttonNonUniformInitiationRates.style('background-color', ButtonDefaultColor);
}
if (getItem("nonUniformInitiationRatesStatus") == 'true')
{
  buttonNonUniformInitiationRates.style('background-color', 'rgb(0, 218, 0)');
}
if (getItem("nonUniformInitiationRatesStatus") == 'false')
{
  buttonNonUniformInitiationRates.style('background-color', 'rgb(255, 0, 0)');
}

// single tRNA modification button
buttontRNAmodification = createButton("modification on single tRNA");
buttontRNAmodification.size(235, 27);
buttontRNAmodification.position(32, 35, 'relative');
buttontRNAmodification.parent("myButtonContainer");
buttontRNAmodification.style('font-size', '14px');
// check button status in local storage: if not empty, recall the stored status:
if (getItem("tRNAmodificationStatus") == null){
  buttontRNAmodification.style('background-color', ButtonDefaultColor);
}
if (getItem("tRNAmodificationStatus") == 'true')
{
  buttontRNAmodification.style('background-color', 'rgb(0, 218, 0)');
}
if (getItem("tRNAmodificationStatus") == 'false')
{
  buttontRNAmodification.style('background-color', 'rgb(255, 0, 0)');
}
// multiple tRNA modifications (ADAT like) button
buttontRNAmultipleModification = createButton("modification on multiple tRNAs (ADAT2)");
buttontRNAmultipleModification.size(270, 27);
buttontRNAmultipleModification.position(16, 70, 'relative');
buttontRNAmultipleModification.parent("myButtonContainer");
buttontRNAmultipleModification.style('font-size', '14px');
// check button status in local storage: if not empty, recall the stored status:
if (getItem("tRNAmultipleModificationStatus") == null){
  buttontRNAmultipleModification.style('background-color', ButtonDefaultColor);
}
if (getItem("tRNAmultipleModificationStatus") == 'true')
{
  buttontRNAmultipleModification.style('background-color', 'rgb(0, 218, 0)');
}
if (getItem("tRNAmultipleModificationStatus") == 'false')
{
  buttontRNAmultipleModification.style('background-color', 'rgb(255, 0, 0)');
}
// limited Pool of ribosomes
buttonLimitedPoolRibosomes = createButton("Limited pool of ribosomes");
buttonLimitedPoolRibosomes.size(180, 27);
buttonLimitedPoolRibosomes.position(1200, 35, 'relative');
buttonLimitedPoolRibosomes.parent("myButtonContainer");
buttonLimitedPoolRibosomes.style('font-size', '14px');
// check button status in local storage: if not empty, recall the stored status:
if (getItem("limitedPoolRibosomeStatus") == null){
  buttonLimitedPoolRibosomes.style('background-color', ButtonDefaultColor);
}
if (getItem("limitedPoolRibosomeStatus") == 'true')
{
  buttonLimitedPoolRibosomes.style('background-color', 'rgb(0, 218, 0)');
  // retrieve the stored slider value: riboPool or sliderSetting
  //sliderSetting = float(localStorage.getItem("riboPoolRatio"));
  //riboSlider.value(sliderSetting);
  getRiboRatioFromSlider();
  // retrieve the black color for riboPoolText
  // retrive the background for riboPoolZone
  highlightRiboPool();
}
if (getItem("limitedPoolRibosomeStatus") == 'false')
{
  buttonLimitedPoolRibosomes.style('background-color', 'rgb(255, 0, 0)');
  unhighlightRiboPool();
}
// expression vector species selection button for coli:
buttonColi = createButton("E. coli");
buttonColi.size(80, 27);
buttonColi.position(760, 7, 'relative');
buttonColi.parent("tRNAHere");
buttonColi.style('font-size', '16px');
// check button status in local storage: if not empty, recall the stored status:
if (getItem("coliSpeciesStatus") == null){
  buttonColi.style('background-color', ButtonDefaultColor);
}
if (getItem("coliSpeciesStatus") == 'true')
{
  buttonColi.style('background-color', 'rgb(0, 218, 0)');
}
if (getItem("coliSpeciesStatus") == 'false')
{
  buttonColi.style('background-color', 'rgb(255, 0, 0)');
}
// expression vector species selection button for yeast:
buttonYeast = createButton("Yeast");
buttonYeast.size(80, 27);
buttonYeast.position(860, 7, 'relative');
buttonYeast.parent("tRNAHere");
buttonYeast.style('font-size', '16px');
// check button status in local storage: if not empty, recall the stored status:
if (getItem("yeastSpeciesStatus") == null){
  buttonYeast.style('background-color', ButtonDefaultColor);
}
if (getItem("yeastSpeciesStatus") == 'true')
{
  buttonYeast.style('background-color', 'rgb(0, 218, 0)');
}
if (getItem("yeastSpeciesStatus") == 'false')
{
  buttonYeast.style('background-color', 'rgb(255, 0, 0)');
}
// expression vector species selection button for worm:
buttonWorm = createButton("Worm");
buttonWorm.size(80, 27);
buttonWorm.position(960, 7, 'relative');
buttonWorm.parent("tRNAHere");
buttonWorm.style('font-size', '16px');

// expression vector species selection button for zebrafish:
buttonZebra = createButton("Zebrafish");
buttonZebra.size(100, 27);
buttonZebra.position(1060, 7, 'relative');
buttonZebra.parent("tRNAHere");
buttonZebra.style('font-size', '16px');

// expression vector species selection button for mouse:
buttonMouse = createButton("Mouse");
buttonMouse.size(80, 27);
buttonMouse.position(1180, 7, 'relative');
buttonMouse.parent("tRNAHere");
buttonMouse.style('font-size', '16px');

// expression vector species selection button for human:
buttonHuman = createButton("Human");
buttonHuman.size(80, 27);
buttonHuman.position(1280, 7, 'relative');
buttonHuman.parent("tRNAHere");
buttonHuman.style('font-size', '16px');

// columns header in latex for the table of the retrieved rates from rates dictionary if chosen species:
//let latex = createP();
//latex.size(8, 8);
//latex.style('font-size', '16px');
//latex.style('background-color', '#008080');
//latex.style('size', (60, 18));
//latex.position(1210, 445, 'relative');
//latex.parent("#canvas_settings_container");
//katex.render('\\tau', latex.elt);

// The function storeData in argument will Store permanently locally in your browser through LocalStorage.
// The magic is that this stored data can be retrieved from another html page!
//dropzone.changed(storeData);

// connect to FIREBASE as a database service:
// Import the functions you need from the SDKs you need
//import { initializeApp } from "firebase/app";
//import { getAnalytics } from "firebase/analytics";
//const { initializeApp } = require('/firebase-admin/app/database');
//import { initializeApp } from "https://www.gstatic.com/firebasejs/9.22.1/firebase-app.js";
//import { getAnalytics } from "https://www.gstatic.com/firebasejs/9.22.1/firebase-analytics.js";
// TODO: Add SDKs for Firebase products that you want to use
// https://firebase.google.com/docs/web/setup#available-libraries

// Your web app's Firebase configuration
// For Firebase JS SDK v7.20.0 and later, measurementId is optional
const firebaseConfig = {
  apiKey: "AIzaSyDXG3A9LZNt0YIDXpqCphwwWHfXsvcEgLU",
  authDomain: "ribosomedigitaltwindb.firebaseapp.com",
  //databaseURL: "https://ribosomedigitaltwindb-default-rtdb.europe-west1.firebasedatabase.app",
  databaseURL: "https://ribosomedigitaltwindb-default-rtdb.europe-west1.firebasedatabase.app.firebaseio.com",
  projectId: "ribosomedigitaltwindb",
  storageBucket: "ribosomedigitaltwindb.appspot.com",
  messagingSenderId: "433909121779",
  appId: "1:433909121779:web:476c181d4f4380ea9d6d86",
  measurementId: "G-4NF3167Q49"
};

// Initialize Firebase (this was done in the settings.html)
//const app = initializeApp(firebaseConfig);
//const analytics = getAnalytics(app);
//const app = initializeApp(firebaseConfig);
//import {getDatabase, set, get, update, remove, ref, child} from "https://www.gstatic.com/firebasejs/9.22.1/firebase-database.js"
//require(getDatabase, set, get, update, remove, ref, child)
//from "https://www.gstatic.com/firebasejs/9.22.1/firebase-database.js"

//firebase.initializeApp(firebaseConfig);
//console.log(firebase);
//console.log(firebaseConfig);

//const app = initializeApp(firebaseConfig);
//const db = getFirestore(app);
//firebase.InitializeApp(firebaseConfig);
//console.log(app)
//console.log(db)

// Get a list of cities from your database
//async function getCities(db) {
//  const citiesCol = collection(db, 'cities');
//  const citySnapshot = await getDocs(citiesCol);
//  const cityList = citySnapshot.docs.map(doc => doc.data());
//  return cityList;
//}

// provide data to the firebase database:
//var ref = database.ref('transcripts');

//var database = firebase.database();
//var ref = database.ref('transcripts');
//var data = {
//  geneID: "COL1",
//  mRNAseq: "UUCGAAAAUUU"
//}
//ref.push(data);
} // end setup()

//------------------------------------------------------------------------------
// draw
//------------------------------------------------------------------------------
function draw(){
  /*background(0, 90, 90); *///Gainsboro 220, 220, 220 or teal green
  /*textSize(18);*/
  /*fill(255);*/
  background(255, 255, 255); //white-background
  textSize(18);
  fill(90);  //black = 0, gray=10-100, white=255

  textFont("Helvetica");
  if (getItem("1") == null){
    text(textProcessed, 10, 25);
  }
  textSize(14);

  // enable the activation/desactivation of buttons for the factors:
  buttonProline.mousePressed(prolineOnOff);
  buttonEFPdeplete.mousePressed(efpDepleteOnOff);
  buttonExitTunnel.mousePressed(ExitTunnelOnOff);
  buttonSecondaryStructure.mousePressed(SecondaryStructureOnOff);
  buttonNonUniformmRNAabundance.mousePressed(NonUniformmRNAabundanceOnOff);
  buttonNonUniformInitiationRates.mousePressed(NonUniformInitiationRatesOnOff);
  buttontRNAmodification.mousePressed(tRNAmodificationOnOff);
  buttontRNAmultipleModification.mousePressed(tRNAmultipleModificationOnOff);
  buttonLimitedPoolRibosomes.mousePressed(limitedPoolRibosomesOnOff);
  buttonColi.mousePressed(coliSpeciesOnOff);
  buttonYeast.mousePressed(yeastSpeciesOnOff);

  // test the status of all the buttons:
  prolineStatus = str(buttonProline.style('background-color')=='rgb(0, 218, 0)')

  EFPdepleteStatus = str(buttonEFPdeplete.style('background-color')=='rgb(0, 218, 0)')
  if (EFPdepleteStatus == 'true')
  {
    //console.log('EF-P deplete button is GREEN and status is ACTIVE.');
  }
  else
  {
    //console.log('EF-P deplete button is red and status is not active.');
  }

  exitTunnelStatus = str(buttonExitTunnel.style('background-color')=='rgb(0, 218, 0)')
  if (exitTunnelStatus == 'true')
  {
    //console.log('Exit tunnel interaction button is GREEN and status is ACTIVE.');
  }
  else
  {
    //console.log('Exit tunnel interaction button is red and status is not active.');
  }

  secondaryStructureStatus = str(buttonSecondaryStructure.style('background-color')=='rgb(0, 218, 0)')
  if (secondaryStructureStatus == 'true')
  {
    //console.log('mRNA secondary structure button is GREEN and status is ACTIVE.');
  }
  else
  {
    //console.log('mRNA secondary structure button is red and status is not active.');
  }

  nonUniformmRNAabundanceStatus = str(buttonNonUniformmRNAabundance.style('background-color')=='rgb(0, 218, 0)')
  if (nonUniformmRNAabundanceStatus == 'true')
  {
    //console.log('non uniform mRNA relative abundance button is GREEN and status is ACTIVE.');
  }
  else
  {
    //console.log('non uniform mRNA relative abundance button is red and status is not active.');
  }

  nonUniformInitiationRatesStatus = str(buttonNonUniformInitiationRates.style('background-color')=='rgb(0, 218, 0)')
  if (nonUniformInitiationRatesStatus == 'true')
  {
    //console.log('non uniform initiation rates button is GREEN and status is ACTIVE.');
  }
  else
  {
    //console.log('non uniform initiation rates button is red and status is not active.');
  }

  tRNAmodificationStatus = str(buttontRNAmodification.style('background-color')=='rgb(0, 218, 0)')
  if (tRNAmodificationStatus == 'true')
  {
    //console.log('tRNA modification button is GREEN and status is ACTIVE.');
  }
  else
  {
    //console.log('tRNA modification button is red and status is not active.');
  }
  tRNAmultipleModificationStatus = str(buttontRNAmultipleModification.style('background-color')=='rgb(0, 218, 0)');

  limitedPoolRibosomeStatus = str(buttonLimitedPoolRibosomes.style('background-color')=='rgb(0, 218, 0)');
  if (limitedPoolRibosomeStatus == 'true')
  {
    //console.log('limited pool of ribosomes button is GREEN and status is ACTIVE.');
  }
  else
  {
    //console.log('limited pool of ribosomes button is red and status is not active.');
  }
  coliSpeciesStatus = str(buttonColi.style('background-color')=='rgb(0, 218, 0)');
  yeastSpeciesStatus = str(buttonYeast.style('background-color')=='rgb(0, 218, 0)');

  // store all the button status in local storage to be available in ribosomer.html:
  storeItem("prolineStatus", prolineStatus);
  //console.log("prolineStatus is currently ", getItem("prolineStatus"));
  storeItem("EFPdepleteStatus", EFPdepleteStatus);
  //console.log("EFPdepleteStatus is currently ", getItem("EFPdepleteStatus"));
  storeItem("exitTunnelStatus", exitTunnelStatus);
  //console.log("exitTunnelStatus is currently ", getItem("exitTunnelStatus"));
  storeItem("secondaryStructureStatus", secondaryStructureStatus);
  //console.log("secondaryStructureStatus is currently ", getItem("secondaryStructureStatus"));
  storeItem("nonUniformmRNAabundanceStatus", nonUniformmRNAabundanceStatus);
  //console.log("nonUniformmRNAabundanceStatus is currently ", getItem("nonUniformmRNAabundanceStatus"));
  storeItem("nonUniformInitiationRatesStatus", nonUniformInitiationRatesStatus);
  //console.log("nonUniformInitiationRatesStatus is currently ", getItem("nonUniformInitiationRatesStatus"));
  storeItem("tRNAmodificationStatus", tRNAmodificationStatus);
  //console.log("tRNAmodificationStatus is currently ", getItem("tRNAmodificationStatus"));
  storeItem("tRNAmultipleModificationStatus", tRNAmultipleModificationStatus);
  //console.log("tRNAmultipleModificationStatus is currently ", getItem("tRNAmultipleModificationStatus"));
  storeItem("limitedPoolRibosomeStatus", limitedPoolRibosomeStatus);
  //console.log("limitedPoolRibosomeStatus is currently ", getItem("limitedPoolRibosomeStatus"));
  storeItem("coliSpeciesStatus", coliSpeciesStatus);
  //console.log("coliSpeciesStatus is currently ", getItem("coliSpeciesStatus"));
  storeItem("yeastSpeciesStatus", yeastSpeciesStatus);
  //console.log("yeastSpeciesStatus is currently ", getItem("yeastSpeciesStatus"));

  // provide gene id tags, length of peptide, name of gene
  // ... while getItem not empty keep the previously stored data in the local storage.
  if (localStorage.getItem("1") == null)
  {
    console.log('The local storage is null at the moment.');
  }

  else
  { var kStorage = localStorage.getItem("1");
    console.log('There is something in local storage, e.g. : ', getItem("1"));
    numbLines = getItem('numbLines');
    for (i=0; i < numbLines; i++)
    {
      output22[i] = getItem(str(i));
      output11[i] = getItem('description' + str(i));
      ORFflag[i] = getItem('#orf#' + str(i));
      //print('ORFflag[i] = ', ORFflag[i]);
      lastCodon[i] = getItem('lastcodon' + str(i));
      myRNAs[i] = getItem(getItem(str(i)));
    }
    if (getItem("readCopiesOf_0") == null || nonUniformmRNAabundanceStatus == 'false'){
      console.log('no copy numbers of transcripts were provided at this moment. Or they are ignored by Button status.');
      for (i=0; i<numbLines; i++){
        readCount[i] = 1;
        readInitRate[i] = 1;
      }
    }
    else
    {
      for (i=0; i<numbLines; i++){
        readCount[i] = getItem('readCopiesOf_'+str(i));
        readInitRate[i] = getItem('initRateOf_'+str(i));
      }
    }// end inner else
  }// end outer else

  if (output22[0] != null){
    text('Here is the list of the first '+ str(idTagsMaxNumb)+ ' provided transcripts variants (capitalized letters gene ID tags): ', 30, 50);
    text('ORF', 927, 50);
    text('last codon', 989, 50);
    text('peptide length', 1064, 50);
    text('mRNA length', 1175, 50);
    text('read copy #', 1285, 50);
    text('init. rate fold change', 1385, 50);
    upperLim = max(idTagsMaxNumb, output22.length);
    for(i=0; i<idTagsMaxNumb; i++){
      text(output22[i].substring(1, output22[i].length-1).toUpperCase(), 30, 80+20*i);
      text(output11[i].substring(1, output11[i].length), 100, 80+20*i);
      text(ORFflag[i], 927, 80+20*i);
      text(lastCodon[i], 989, 80+20*i);
      if (myRNAs[i].length % 3 == 0){
        peplength = (myRNAs[i].length-3)/3;
      }
      else {
        peplength = '?'
      }
      text(str(peplength), 1064, 80+20*i);
      text(myRNAs[i].length, 1175, 80+20*i);
      if (readCount[i] != null){
        text(readCount[i], 1285, 80+20*i);
        text(nfs(readInitRate[i], 2, 1), 1385, 80+20*i);
      }
    }
  }

  // see function callback 'getRateFromInput()':
  //rateSetting.changed(getRateFromInput, highlightRateBox);
  var goodHalfTime = str(nfs(rateSettingValue * 0.6931472, 3, 1)); // * Ln2
  textHalfTime.html(goodHalfTime);

  //highlightRateBox();

  // store the initiation time rate in ms in local storage for later retrieval in sketch.js:
  var initBETA = float(rateSettingValue*1000); // converts in milliseconds before sending to localStorage
  storeItem("initiationRateTime", initBETA);

  // see function callback 'getRiboRatio'
  riboText = getRiboRatioFromSlider();
  riboPoolRatioValue = float(riboText);
  storeItem("riboPoolRatio", riboPoolRatioValue);
  //storeItem("riboPoolRatio", sliderSetting);

  riboPoolText.html(riboText);

   // see function callback 'writeCodon
  var goodText = writeCodon();
  textCodon.html(goodText);
  // deduce the amino acid residue associated to this codon:
  var aa_res = decode();
  //textResidue.html(aa_res);
  var XXX = threeLetterCode(aa_res);
  var stringResidue = aa_res + ' - ' + XXX;
  textResidue.html(stringResidue);
  //console.log('output22 length:', output22.length);

    // Build table of elongation rate with defaulted parameters from database for the chosen species
  // box:
  push();
  strokeWeight(2);
  stroke(150); // stroke(100, 0, 170)
  var x_bottom_left = pxSup2 + 30; //810
  var y_bottom_left = pySup; //680
  var y_line_length = 330;
  var x__line_length = 420;
  line(x_bottom_left, y_bottom_left, x_bottom_left, y_bottom_left-y_line_length);
  line(x_bottom_left, y_bottom_left, x_bottom_left + x__line_length, y_bottom_left);
  line(x_bottom_left + x__line_length, y_bottom_left, x_bottom_left + x__line_length, y_bottom_left-y_line_length);
  line(x_bottom_left, y_bottom_left - y_line_length, x_bottom_left + x__line_length, y_bottom_left-y_line_length);

  pop();
  // title:
  chosenSpecies ='unknown';
  var lineOne = "Synonymous codons elongation rates retrieved parameters of";
  if (getItem("coliSpeciesStatus")=='true'){
    chosenSpecies = 'coli';
  }
  if (getItem("yeastSpeciesStatus")=='true'){
    chosenSpecies = 'yeast';
  }
  // retrieve the codon complementary to the anti-codon provided in the tRNA modification box:
  var currentTriplet = str(writeCodon());
  //console.log('currentTriplet modified or not from input=', currentTriplet);
  // determine its residue and in 3 letters qs well:
  var currentAA = decode(currentTriplet);
  //console.log('that is our currentAA=', currentAA);
  var currentAA3letters = str(threeLetterCode(currentAA));
  var lineTwo = 'amino acid residue ' + str(currentAA) + '-' + str(currentAA3letters) + ' for ' + str(chosenSpecies) + ' species:';
  textFont("Helvetica"); //Courier or Helvetica
  textSize(14);
  fill(0, 90, 90); //(75, 0, 130)
  //text(lineOne, x_bottom_left + 30, 430);
  text(lineOne, x_bottom_left + 15, 430);
  text(lineTwo, x_bottom_left + 15, 446);
  //text('synonymous codons', x_bottom_left + 15, 470)
  var lineThree = 'synonymous codons ' + '\t' + '\t'+   '  tau'+'_1'  + '\t' + '   tau'+'_2' + '\t'+ '\t' + '  tau' + '_3' + '\t'+ '\t' + '\t' + '\t '+ '      tot';
  text(lineThree, x_bottom_left + 15, 470);
  var lineFour = '(ms)' + '\t' + '\t' + '  (ms)' + '\t' + '\t' + '    (ms)' + '\t' + '\t' + '\t' +'\t' + '    (ms)';
  text(lineFour, x_bottom_left + 180, 486);
  //console.log('aa=', currentAA, 'has synonymous codons:', findSynonymousList(currentAA));
  //var colorChoice = ["255, 0, 0", "0, 255, 0", "0, 0, 255", "0, 0, 0", "100, 100, 100", "255, 0, 255"]; // red, green, blue, black, grey, magenta
  var colorChoice;
  for (i=0; i<findSynonymousList(currentAA).length;i++){
    textFont("Helvetica"); //Courier or Helvetica
    textSize(14);
    if (i==0){
      colorChoice='red';
    }
    if (i==1){
      colorChoice = 'green';
    }
    if (i==2){
      colorChoice = 'blue';
    }
    if (i==3){
      colorChoice = 'black';
    }
    if (i==4){
      colorChoice = 'grey';
    }
    if (i==5){
      colorChoice = 'magenta';
    }
    fill(colorChoice);
    text(str(findSynonymousList(currentAA)[i]), x_bottom_left + 60, 524 + i*32);
    var paramRates = [];
    var totTime;
    if (yeastSpeciesStatus=='true'){
      var rateTrigram = findSynonymousList(currentAA)[i];
      paramRates = yeastHypoExpParams[rateTrigram];
      totTime = 1.0/paramRates[0] + 1.0/paramRates[1] + 1.0/paramRates[2];
    }//end if
    if (coliSpeciesStatus=='true'){
      var rateTrigram = findSynonymousList(currentAA)[i];
      paramRates = coliHypoExpParams[rateTrigram];
      totTime = 1.0/paramRates[0] + 1.0/paramRates[1] + 1.0/paramRates[2];
    }//end if
    var lineRates = nf(1.0/paramRates[0], 4, 1) +'\t'+'\t'+ nf(1.0/paramRates[1], 4, 1) +'\t'+'\t'+ nf(1.0/paramRates[2], 4, 1) +'\t'+'\t'+'\t'+'\t'+ nf(totTime, 4, 1);
    text(lineRates, x_bottom_left + 160, 524 + i*32);
  }// end for

  //******
  //  display the hypoexponential function in graphic box:
  // draw the border of the graphic box:
  push();
  strokeWeight(4);
  stroke(0, 90, 90);
  line(pxInf, pyInf, pxSup, pyInf);
  line(pxInf, pyInf, pxInf, pySup);

  line(pxInf, pySup, pxSup, pySup);
  line(pxSup, pyInf, pxSup, pySup);
  pop();

  push();
  strokeWeight(2);
  stroke(150);
  line(35, 710, 490, 710);
  line(35, 440, 35, 710);

  line(35, 1850, 490, 1850);//**************

  pop();

  // graphic box title and axes labels
  var lineOne = "Ribosome residence time distribution for modified tRNA";
  textFont("Helvetica"); //Courier or Helvetica
  textSize(14);
  fill(0, 90, 90); //(75, 0, 130)
  text(lineOne, 70, 430);

  var ord_label = "pdf";
  textFont("Helvetica"); //Courier or Helvetica
  textSize(14);
  fill(100, 105, 100); //(75, 0, 130)
  text(ord_label, 27, 434);

  var abs_label = "time (ms)";
  textFont("Helvetica"); //Courier or Helvetica
  textSize(14);
  fill(100, 105, 100); //(75, 0, 130)
  text(abs_label, 438, 725);

  // draw ticks :
  push();
  strokeWeight(2);
  stroke(150); // stroke(100, 0, 170)
  var ticklength = 7;
  line(px_tick[0], pySup - 30 - ticklength, px_tick[0], pySup - 30);
  line(px_tick[1], pySup - 30 - ticklength, px_tick[1], pySup - 30);
  line(px_tick[2], pySup - 30 - ticklength, px_tick[2], pySup - 30);
  line(px_tick[3], pySup - 30 - ticklength, px_tick[3], pySup - 30);
  line(px_tick[4], pySup - 30 - ticklength, px_tick[4], pySup - 30);
  line(px_tick[5], pySup - 30 - ticklength, px_tick[5], pySup - 30);
  pop();

  // scale for x-axis ticks:
  push();
  noStroke();
  fill(120, 125, 120);
  textSize(12);
  text(nfs(50, 2, 0), px_tick[0] - 10, pySup - 15);
  text(nfs(100, 3, 0), px_tick[1] - 12, pySup - 15);
  text(nfs(150, 3, 0), px_tick[2] - 12, pySup - 15);
  text(nfs(200, 3, 0), px_tick[3] - 12, pySup - 15);
  text(nfs(250, 3, 0), px_tick[4] - 12, pySup - 15);
  pop();

  // initialize 2000 points of hypoexponential functions (ref and mod):
  numbPointstoDraw = 2000;
  dt_ms = tSup/numbPointstoDraw;
  tauList = [];
  input_tRNA_anticodon.changed(writeCodon, decode);
  currentTriplet = str(writeCodon());
  //console.log('here is the current triplet in U mRNA', currentTriplet);
  currentCDStriplet = mRNAtoCDS(currentTriplet);
  //console.log('here is the current triplet in CDS', currentCDStriplet);
  legitTriplet = Object.keys(yeastHypoExpParams);
  if (yeastSpeciesStatus == 'true' && legitTriplet.includes(currentCDStriplet)){
    paramRates = yeastHypoExpParams[currentCDStriplet];
    tauList.push(1.0/paramRates[0]);
    tauList.push(1.0/paramRates[1]);
    tauList.push(1.0/paramRates[2]);
  }//end if
  if (coliSpeciesStatus == 'true' && legitTriplet.includes(currentCDStriplet)){
    paramRates = coliHypoExpParams[currentCDStriplet];
    tauList.push(1.0/paramRates[0]);
    tauList.push(1.0/paramRates[1]);
    tauList.push(1.0/paramRates[2]);
  }//end if

  for (i=0; i<numbPointstoDraw; i++){
    t_ms = tInf + i * dt_ms;
    t_2draw[i] = t_ms;
    pdf_mod[i] = hypoexponential(t_ms, tRNAmod_tau);
    pdf_ref[i] = hypoexponential(t_ms, tauList);
    }

  // find max of the two pdfs to be used to scale the plot as best as possible.
  pdf_mod_max = max(pdf_mod);
  pdf_ref_max = max(pdf_ref);
  var ySup_pdf = max(pdf_mod_max, pdf_ref_max);

  //scale the y-axis with map and the canvas positions:
  var py_tick_max;
  py_tick_max = map(ySup_pdf, 0, ySup_pdf*1.05, 710, 440);
  push();
  strokeWeight(2);
  stroke(150); // stroke(100, 0, 170)
  line(pxInf+25, py_tick_max, pxInf+25+ticklength, py_tick_max);
  pop();
  // legend for the tRNA being displayed:
  push();
  noStroke();
  fill(0, 90, 90); // 120, 125, 120
  textSize(18);
  tRNAcodon = str(currentTriplet);
  tRNAlegend = 'tRNA';
  tRNAaa = str(currentAA)+'-'+str(threeLetterCode(currentAA));
  addXoffset = 260;
  addYoffset = 15;
  text(tRNAlegend, pxInf + addXoffset + 130, py_tick_max+addYoffset+2);
  textSize(12);
  text(tRNAaa, pxInf + addXoffset + 172, py_tick_max+addYoffset-7);
  text(tRNAcodon, pxInf +addXoffset + 102, py_tick_max+addYoffset+7);
  textSize(14);
  fill(0, 0, 0);
  text('standard rates: - black line', pxInf +addXoffset + 45, py_tick_max+55);
  fill(0, 90, 90);
  text('modified rates: - teal line', pxInf +addXoffset + 45, py_tick_max+80);
  pop();

  // scale for y-axis ticks:
  push();
  noStroke();
  fill(120, 125, 120);
  textSize(12);
  text(nfs(ySup_pdf, 1, 4), pxInf + 30, py_tick_max+2);
  pop();

  // draw the points:
  push(); // for the reference distribution (gray)
  for (var i = 0; i < numbPointstoDraw; i++){
    strokeWeight(2);
    stroke(55, 55, 55); // grey
    px_ref[i] = map(t_2draw[i], 0, 300, 35, 490);
    py_ref[i] = map(pdf_ref[i], 0, ySup_pdf*1.05, 710, 440); // ???
    point (px_ref[i], py_ref[i]);
  }
  pop();

  push(); // for the tRNA modified distribution (teal)
  for (var i = 0; i < numbPointstoDraw; i++){
    strokeWeight(2);
    stroke(0, 90, 90);
    px_mod[i] = map(t_2draw[i], 0, 300, 35, 490);
    py_mod[i] = map(pdf_mod[i], 0, ySup_pdf*1.05, 710, 440); // ???
    point (px_mod[i], py_mod[i]);
  }
  pop();

  //******
  //  Synonymous codons elongation cycle time distribution for a given aa
  //  display the hypoexponential function in graphic box:
  //  draw the border of the graphic box:
  push();
  strokeWeight(4);
  stroke(0, 90, 90);
  line(pxInf2, pyInf, pxSup2, pyInf);
  line(pxInf2, pyInf, pxInf2, pySup);

  line(pxInf2, pySup, pxSup2, pySup);
  line(pxSup2, pyInf, pxSup2, pySup);
  pop();

  push();
  strokeWeight(2);
  stroke(150);
  line(575, 710, 1030, 710);
  line(575, 440, 575, 710);

  pop();

  // graphic box title and axes labels
  var lineOne = "Ribosome residence time distribution for synonymous codons";
  textFont("Helvetica"); //Courier or Helvetica
  textSize(14);
  fill(0, 90, 90); //(75, 0, 130)
  text(lineOne, 610, 430);

  var ord_label = "pdf";
  textFont("Helvetica"); //Courier or Helvetica
  textSize(14);
  fill(100, 105, 100); //(75, 0, 130)
  text(ord_label, 567, 434);

  var abs_label = "time (ms)";
  textFont("Helvetica"); //Courier or Helvetica
  textSize(14);
  fill(100, 105, 100); //(75, 0, 130)
  text(abs_label, 978, 725);

  // draw ticks :
  push();
  strokeWeight(2);
  stroke(150); // stroke(100, 0, 170)
  var ticklength = 7;
  line(px_tick2[0], pySup - 30 - ticklength, px_tick2[0], pySup - 30);
  line(px_tick2[1], pySup - 30 - ticklength, px_tick2[1], pySup - 30);
  line(px_tick2[2], pySup - 30 - ticklength, px_tick2[2], pySup - 30);
  line(px_tick2[3], pySup - 30 - ticklength, px_tick2[3], pySup - 30);
  line(px_tick2[4], pySup - 30 - ticklength, px_tick2[4], pySup - 30);
  line(px_tick2[5], pySup - 30 - ticklength, px_tick2[5], pySup - 30);
  pop();

  // scale for x-axis ticks:
  push();
  noStroke();
  fill(120, 125, 120);
  textSize(12);
  text(nfs(50, 2, 0), px_tick2[0] - 10, pySup - 15);
  text(nfs(100, 3, 0), px_tick2[1] - 12, pySup - 15);
  text(nfs(150, 3, 0), px_tick2[2] - 12, pySup - 15);
  text(nfs(200, 3, 0), px_tick2[3] - 12, pySup - 15);
  text(nfs(250, 3, 0), px_tick2[4] - 12, pySup - 15);
  pop();

  // array the pdf functions for all synonymous codons:
  var pdf_syn = [];
  paramRates = [];

  // do the next for loop only if you have a legit trigram or a known aa:
  if (currentAA != 'X'){
    for (j=0; j<findSynonymousList(currentAA).length;j++){// for each synonymous codon j:
      tauList = [];
      pdf_syn[j] = [];
      py_syn[j] = [];
      var rateTrigram = findSynonymousList(currentAA)[j];
      // select the species:
      switch (chosenSpecies){
        case 'coli':
          paramRates = coliHypoExpParams[rateTrigram];
          break;
        case 'yeast':
          paramRates = yeastHypoExpParams[rateTrigram];
          break;
      }
      tauList.push(1.0/paramRates[0]); // retrieve the rates for this codon j
      tauList.push(1.0/paramRates[1]);
      tauList.push(1.0/paramRates[2]);

      for (i=0; i<numbPointstoDraw; i++){//for each point i:
        t_ms = tInf + i * dt_ms;
        t_2draw[i] = t_ms;
        pdf_syn[j].push(hypoexponential(t_ms, tauList));
      }//end inner loop
      var yMax_pdf = max(pdf_syn[j]);
      if (yMax_pdf > ySup_pdf){
        ySup_pdf = yMax_pdf;
      }
      // draw the points:
      push(); // for the synonymous codon retrieved parameters
      switch(j){
        case 0:
          colorChoice='red';
          break;
        case 1:
          colorChoice='green';
          break;
        case 2:
          colorChoice='blue';
          break;
        case 3:
          colorChoice='black';
          break;
        case 4:
          colorChoice='grey';
          break;
        case 5:
          colorChoice='magenta';
          break;
      }
      for (var i = 0; i < numbPointstoDraw; i++){
        strokeWeight(2);
        stroke(colorChoice); // grey (55, 55, 55)
        px_syn[i] = map(t_2draw[i], 0, 300, 575, 1030);
        py_syn[j].push(map(pdf_syn[j][i], 0, ySup_pdf*1.05, 710, 440)); // ???
        point (px_syn[i], py_syn[j][i]);
      }
      pop();
    }// end outer loop

    // scale for y-axis ticks:
    push();
    noStroke();
    fill(120, 125, 120);
    textSize(12);
    text(nfs(ySup_pdf, 1, 4), pxInf2 + 35, py_tick_max+2);
    pop();
    //scale the y-axis with map and the canvas positions:
    push();
    strokeWeight(2);
    stroke(150); // stroke(100, 0, 170)
    line(pxInf2+30, py_tick_max, pxInf2+30+ticklength, py_tick_max);
    pop();

    // legend for the aa of the synonymous codons being displayed:
    push();
    noStroke();
    fill(0, 90, 90); // 120, 125, 120
    tRNAaa = str(currentAA)+'-'+str(threeLetterCode(currentAA));
    addXoffset = 260;
    addYoffset = 15;
    textSize(16);
    text(tRNAaa, pxInf2 + addXoffset + 146, py_tick_max+addYoffset-7);
    pop();
  }// end if

  // an extra distribution profile display for proline depending
  //on proline button status and the EF-P depletion button status
  var slowProlineRate = [];
  if (prolineStatus == 'true' && currentAA == 'P' && legitTriplet.includes(currentCDStriplet)){
    // select the species:
    switch (chosenSpecies){
      case 'coli':
        slowProlineRate[0] = Object.values(coliHypoExpParams[currentCDStriplet])[0];
        slowProlineRate[1] = 0.6931/5102; // ln2 / t_1/2 0.6931/5102
        slowProlineRate[2] = Object.values(coliHypoExpParams[currentCDStriplet])[2];
        break;
      case 'yeast':
        slowProlineRate[0] = Object.values(yeastHypoExpParams[currentCDStriplet])[0];
        slowProlineRate[1] = 0.6931*0.5/5102; // ln2 / t_1/2   //0.6931/120
        slowProlineRate[2] = Object.values(yeastHypoExpParams[currentCDStriplet])[2];
        break;
    } // end switch cases
  } // end if

  if (EFPdepleteStatus == 'true' && currentAA == 'P' && legitTriplet.includes(currentCDStriplet)){
    // select the species:
    switch (chosenSpecies){
      case 'coli':
        slowProlineRate[0] = Object.values(coliHypoExpParams[currentCDStriplet])[0];
        //slowProlineRate[1] = Object.values(coliHypoExpParams[currentCDStriplet])[1] / 90.0; // it is 90 times slower
        slowProlineRate[1] = 0.6931/(5102*90.0); // it is still 90 times slower than with EF_P
        slowProlineRate[2] = Object.values(coliHypoExpParams[currentCDStriplet])[2];
        break;
      case 'yeast':
        slowProlineRate[0] = Object.values(yeastHypoExpParams[currentCDStriplet])[0];
        //slowProlineRate[1] = Object.values(yeastHypoExpParams[currentCDStriplet])[1] / 90.0; // ln2 / t_1/2
        slowProlineRate[1] = 0.6931*0.5 / (5102*90.0);
        slowProlineRate[2] = Object.values(yeastHypoExpParams[currentCDStriplet])[2];
        break;
    } // end switch cases
  } // end if
  // calculate the pdf for slow down proline to draw:
  tauPro = [];
  tauPro.push(1.0/slowProlineRate[0]); // retrieve the rates for this proline codon with slow down
  tauPro.push(1.0/slowProlineRate[1]);
  tauPro.push(1.0/slowProlineRate[2]);

  // display the slowdown proline profile according to above proline button status:
  if ((prolineStatus=='true' || EFPdepleteStatus =='true') && currentAA == 'P' && legitTriplet.includes(currentCDStriplet)){
    push();
    pdf_pro = [];
    py_pro = [];
    for (var i = 0; i < numbPointstoDraw; i++){
      t_ms = tInf + i * dt_ms;
      t_2draw[i] = t_ms;
      pdf_pro.push(hypoexponential(t_ms, tauPro));
      strokeWeight(2);
      stroke('orange'); // grey (55, 55, 55) or colorChoice = 'orange'
      px_syn[i] = map(t_2draw[i], 0, 300, 575, 1030);
      py_pro.push(map(pdf_pro[i], 0, ySup_pdf*1.05, 710, 440)); // ???
      point(px_syn[i], py_pro[i]);
    }
    pop();
  }

} //end function draw
//-------END-------------------------------------------------------------------
