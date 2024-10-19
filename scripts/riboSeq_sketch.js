// riboSeq_sketch
var riboRRT_string;
var riboSP_string;
var riboFC_string = '';
var riboRRT = []; // to be a 2D array
var riboSP = []; // to be a 2D array
var riboFC = []; // to be a 2D array
var snapshots_nb;
var footprintedFragmentsCount;
var transcriptsLengthInCodons;
var geneTag;
var tagInput;
var tagNameBox;
var tagInputName = '?';
var submitButton;
var TagNameIDstring = '';
var geneTagString = '';
var geneName ='';
var geneIndex = 0;
var read_indices = [];
var readNumber;
var transcriptsFullNb;
var readTagString;
var iMax;
var keyOne;
var StoredGeneID = []; //without parenthèses
var geneTagID = []; // with the parenthèses
var StoredTransCount = [];
var myRNA = [];
var RNAseq;
var codonSeqList = [];
var codonOrderedNumb = [];
var aaSeqList = [];
var RRT = []; // total time spent by all ribosomes on each codon of a given transcript.
var RRT_type_1 = []; // time in ms spent by a single ribosome on average on a codon of a specified transcript.
var RRT_type_2 = []; // time fraction in 1 per thousand unit a single ribosome spent on a codon of a specified transcript.
var SP = []; // total timer set point attributed to all ribosomes on each codon of a given transcript.
var SP_type_1 = [];
var SP_type_2 = [];
var FC = []; // ribosome footprint count normalized by number of transcript read copies for a specified transcript
var text_codon = '';
var text_aa = '';
var buttonLSsave;
var defaultColor = 'rgb(255, 0, 0)'; //'rgb(150, 150, 150, 125)';
var fileName = "LSbu00";
var dropzone4LS;
var parag;
//var list_Of_keys = [];
var keysForLS = [];
var valuesForLS = [];
var list_Of_values = [];
var list_Of_pairsKeyValue = [];
var reg_key;
var reg_key_value;
var reg_value;
var splitByCommas;
var splitBySemiColon;
var ls_json;
var i_length = 36;
var StoredProtCount = []; // to be instantiated later as a 2D array
var StoredTimes = [];
var StoredFCcount = []; //to be instantiated later as a 2D array
var string_array = '';
var arrayOfStrings = [];
var ButtonDefaultColor = 'rgb(255, 0, 0)'; //red
var ButtonColorGreen = 'rgb(0, 218, 0)'; //green
// graph and map variables:
var pxInf = 90;
var pyInf = 500;
var pyInfMap = 240;
var pxSup = 1465;
var pySup = 50;
var dx_triplet = 36;
var px_occ;
var py_occ;
var lower_px = pxInf - 12;
var lower_codon = 0;
var start_codon = 0;
let rangeSlider;
var window_range = 130; //130
var y_mapped_max;
var rrtStatus;
var spStatus;
var fcStatus;

//---------------------------------------------------------------------------
// genetic code 3D-array: 4 by 4 by 4 3D array with 64 elements:
var GenetiCode = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]];
//------------------------------------------------------------------------
// function section:
//-------------------------------------------------------------------------
//function preload(){
//  LS_bu = loadJSON("http://localhost:8080/output/SimulationResults/LSbu00.txt");
//}
function highlight(){
  dropzone4LS.style("background-color", "grey");
}
function unhighlight(){
  dropzone4LS.style('background-color', '#C8C4');//#C8C4
}
function retrieveLS(file){
  // callback function upon dropping a file (json type with a previously stored local storage):
  //takes the following actions:
  // 1) clear the current local storage:
  clearStorage();
  console.log('The local storage was cleared.');
  //console.log(localStorage);

  // 2) retrieve the saved (back-up) LS json file and write the local storage keys and
  // values pairs to the active local storage:

  if (file.type != 'text'){
    parag = createP("This file is wrong: a text file is required in JSON format with all local storage of the key-value pairs. ");
    parag.parent('#canvasRiboSeq'); //'#canvas_settings_container' or '#fastaHere'
    parag.style("fontSize", "16pt");
    parag.style("family-font", "Lucida Console");
    parag.style("textAlign", "center");
    parag.style("color", "white");
    parag.style("background-color", "red");
    parag.position(25, 600, 'relative');
    highlight();
  }
  else if (file.type == 'text') {
    if (parag != null){
      parag.remove();
    }
    // if the file is a text, do something specific with its content:
    var ls_json_full = str(file.data);
    ls_json = ls_json_full.substring(2, ls_json_full.length-2);
    //console.log('ls_json=', ls_json);

    // we use regex (with capturing group or not) to retrieve all keys and all values sequentially:
    // all keys starts with either [literal { or litteral,] and always ends with ":".
    // all values always starts with ":" and ends with litteral" and [litteral, or litteral }].
    //reg_key = /[{|,]\"\#?\(?\w*#?\)?\w*\":/g;
    //reg_key = /[{|,]\"(\#?\(?\w*#?\)?\w*)\":/g; // capturing group
    //reg_key_value = /[{|,]\"\#?\(?\w*#?\)?\w*\":"([([]?\w*\d*\,*.*\]?)\"[,}]/g; // capturing group
    //reg_key_value = /[{|,]\"\#?\(?\w*#?\)?\w*\":"[([]?\w*\d*\,*.*\]?\"[,}]/g;
    //reg_value = /:\"([([]?\w*\d*\,*.*\]?\"?\,?\"?[,}]?/g;
    //reg_value = /:\"[([]?\w*\d*\,*.*\]?\"?\,?\"?[,}]?/g;
    splitByCommas = /\"\,\"/g;
    //splitBySemiColon = /\":\"/g;

    list_Of_pairsKeyValue = ls_json.split(splitByCommas);

    for (var i = 0; i<list_Of_pairsKeyValue.length; i++){
      var indexMaxKeys = list_Of_pairsKeyValue[i].match(/.*\":\"/)[0].length;
      keysForLS.push(list_Of_pairsKeyValue[i].match(/.*\":\"/)[0].substring(0, indexMaxKeys-3));

      var indexMaxValues = list_Of_pairsKeyValue[i].match(/\":\".*/)[0].length;
      valuesForLS.push(list_Of_pairsKeyValue[i].match(/\":\".*/)[0].substring(3, indexMaxValues));
    }
    for (var i = 0; i<list_Of_pairsKeyValue.length; i++){
      //console.log('key=', keysForLS[i], 'value=', valuesForLS[i])
      // generate the local storage from the retrieved file:
      storeItem(str(keysForLS[i]), valuesForLS[i]);
    }
    console.log('The local storage was updated from file drop.');
    //console.log(localStorage);
  } // end else
}// end function retrieveLS

function saveLS(){
  // call back function upon clicking on button "save local storage to file".
  // Take the following action:
  // fileName when LS is saved:
  fileName = "LSbu00";
  // change the background color of the button:
  if (buttonLSsave.style('background-color') == defaultColor){
    buttonLSsave.style('background-color', 'rgb(0, 255, 0)');
  }
  else{
    buttonLSsave.style('background-color', defaultColor);
  }
  /* dump local storage to string */

  var a = {};
  for (var i = 0; i < localStorage.length; i++) {
    var k = localStorage.key(i);
    var v = localStorage.getItem(k);
    a[k] = v;
  }

  /* save as blob */
  var textToSave = JSON.stringify(a)
  var textToSaveAsBlob = new Blob([textToSave], {
    type: "text/plain"
  });
  var textToSaveAsURL = window.URL.createObjectURL(textToSaveAsBlob);

  /* download without button hack */
  var downloadLink = document.createElement("a");
  downloadLink.download = fileName;
  downloadLink.innerHTML = "Download File";
  downloadLink.href = textToSaveAsURL;
  downloadLink.onclick = function () {
    document.body.removeChild(event.target);
  };
  downloadLink.style.display = "none";
  document.body.appendChild(downloadLink);
  downloadLink.click();
}// end saveLS() function

function getTagNameID(){
  geneName = tagInput.value();
  TagNameIDstring = str(geneName);
  geneTagString = "(" + TagNameIDstring + ")";
  geneIndex = int(geneTagID.indexOf(geneTagString));
  if (geneIndex < 0){
    TagNameIDstring = "WRONG TAG NAME..."+str(geneTagID[0]).substring(1, geneTagID[0].length - 1)+" was taken by default.";
    tagInputName = geneTagID[0].substring(1, geneTagID[0].length -1);
  }
  else{
    tagInputName = geneTagID[geneIndex].substring(1, geneTagID[geneIndex].length -1);
    TagNameIDstring = str(tagInputName);
    geneTagString = "(" + TagNameIDstring + ")";
  }
  tagNameBox.html(TagNameIDstring);
  tagNameBox.style("color", "rgb(255, 255, 255)");
  tagNameBox.style('background-color', '#CD7D6B');/* '#CD7D6B' */
} // end getTagNameID() function

function indexOfRead(geneIDlist, readCopyList, readIDstring){
  //This function returns the index (integer) of the transcript read list (tr index)
  // corresponding to the geneIDlist given in first argument, the readCopyList given
  // in second argument and the unique read ID string given as the third argument (uniqueID to query for).
  var itemList = [];
  for (i=0; i<geneIDlist.length; i++){
    for (j=0; j<readCopyList[i]; j++){
      itemList.push("("+geneIDlist[i]+")#"+str(j));
    }
  }
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

function buttonRRTOnOff(){
  // color the button in green/red
  // test button.style value:
  if (buttonRRT.style('background-color') == 'rgb(255, 0, 0)')
  {
    buttonRRT.style('background-color', ButtonColorGreen);
    buttonSP.style("background-color", ButtonDefaultColor);
    buttonFC.style("background-color", ButtonDefaultColor);
  }
  else
  {
    buttonRRT.style("background-color", ButtonDefaultColor);
  }
}

function buttonSPOnOff(){
  // color the button in green/red
  // test button.style value:
  if (buttonSP.style('background-color') == ButtonDefaultColor)
  {
    buttonSP.style('background-color', ButtonColorGreen);
    buttonRRT.style("background-color", ButtonDefaultColor);
    buttonFC.style("background-color", ButtonDefaultColor);
  }
  else
  {
    buttonSP.style("background-color", ButtonDefaultColor);
  }
}

function buttonFCOnOff(){
  // color the button in green/red
  // test button.style value:
  if (buttonFC.style('background-color') == 'rgb(255, 0, 0)')
  {
    buttonFC.style('background-color', ButtonColorGreen);
    buttonRRT.style("background-color", ButtonDefaultColor);
    buttonSP.style("background-color", ButtonDefaultColor);
  }
  else
  {
    buttonFC.style("background-color", ButtonDefaultColor);
  }
}

function occupancyMapShow(array_Y, pixMaxInf_y, pixMaxSup_y){
  // find max:
  var occupancyMax = max(array_Y);
  // iterate x-positions
  for (i=0; i<array_Y.length; i++){
    px_occ = pxInf - 12 + 2 + i*3;
    // scale on y-axis to find y coordinate:
    py_occ = map(array_Y[i], 0, occupancyMax * 1.15, 0, pixMaxInf_y - pixMaxSup_y); // pixMaxInf_y, pixMaxSup_y
    push();
    strokeWeight(1);
    stroke(200); // stroke(100, 0, 170)
    line(px_occ, pyInfMap - 40, px_occ, pyInfMap -40 - py_occ);
    stroke(0, 90, 90);// color of point (teal)
    strokeWeight(4);
    point(px_occ, pyInfMap -40 - py_occ);
    pop();
  }
  // add y-scale axis and axis label:
  y_mapped_max = map(occupancyMax, 0, occupancyMax * 1.15, 0, pixMaxInf_y - pixMaxSup_y);
  push();
  strokeWeight(2);
  stroke(0, 90, 90); // stroke(100, 0, 170) 50
  line(pxInf - 24, pyInfMap - 40, pxInf - 24, pyInfMap -40 - y_mapped_max);
  pop();
  if (rrtStatus == 'true'){
    push();
    textFont("Helvetica"); //Courier or Helvetica or Calibri
    textSize(14);
    fill(0, 90, 90); //(75, 0, 130) // coral color 255, 127, 80
    text("RRT (ms)", pxInf - 50, pyInfMap -40 - y_mapped_max - 16); // if RRT_type_1
    //text("RRT [-] (x 1,000)", pxInf - 32, pyInfMap -40 - y_mapped_max - 16); // if RRT_type_2
    textSize(12);
    text(str(nfs(occupancyMax, 4,1)), pxInf - 50, pyInfMap -40 - y_mapped_max + 5);
    pop();
  }
  if (spStatus == 'true'){
    push();
    textFont("Helvetica"); //Courier or Helvetica or Calibri
    textSize(14);
    fill(0, 90, 90); //(75, 0, 130) // coral color 255, 127, 80
    text("timer SP (ms)", pxInf - 40, pyInfMap -40 - y_mapped_max - 16); // if SP_type_1
    //text("RRT [-] (x 1,000)", pxInf - 32, pyInfMap -40 - y_mapped_max - 16); // if SP_type_2
    textSize(12);
    text(str(nfs(occupancyMax, 4,1)), pxInf - 60, pyInfMap -40 - y_mapped_max + 5);
    pop();
  }
  if (fcStatus == 'true'){
    push();
    textFont("Helvetica"); //Courier or Helvetica or Calibri
    textSize(14);
    fill(0, 90, 90); //(75, 0, 130) // coral color 255, 127, 80
    text("NFC", pxInf - 50, pyInfMap -40 - y_mapped_max - 16); // if RRT_type_1
    //text("RRT [-] (x 1,000)", pxInf - 32, pyInfMap -40 - y_mapped_max - 16); // if RRT_type_2
    textSize(12);
    text(str(nfs(occupancyMax, 4,1)), pxInf - 50, pyInfMap -40 - y_mapped_max + 5);
    pop();
  }

}// end function occupancyMapShow()

function zoomOccupancy(array_Y, start_codon, window_range, py_inf, py_sup){
  // find max of the range:
  var occupancyMaxRange = max(array_Y.slice(start_codon, start_codon + window_range));
  // iterate x-positions (1 data y on px = dx_triplet/2, px+dx_triplet , etc...) in the scale dx_triplet per codon:
  for (i=0; i<window_range; i++){
    px_occ = pxInf - 18 + dx_triplet*(i + 0.5);
    // scale on y-axis to find y coordinate:
    py_occ = map(array_Y[start_codon+i], 0, occupancyMaxRange * 1.15, 0, py_inf - py_sup);
    push();
    strokeWeight(4);
    stroke(200); // stroke(100, 0, 170)
    line(px_occ, pyInf-8, px_occ, pyInf-8 - py_occ);
    stroke(0, 90, 90);// color of point (teal)
    strokeWeight(10);
    point(px_occ, pyInf-8 - py_occ);
    pop();
  }
  // add y-scale axis and axis label:
  y_mapped_max = map(occupancyMaxRange, 0, occupancyMaxRange * 1.15, 0, py_inf - py_sup);
  push();
  strokeWeight(4);
  stroke(0, 90, 90); // stroke(100, 0, 170) 50
  line(pxInf - 24, pyInf - 8, pxInf - 24, pyInf-8 - y_mapped_max);
  pop();
  if (rrtStatus == 'true'){
    push();
    textFont("Helvetica"); //Courier or Helvetica or Calibri
    textSize(14);
    fill(0, 90, 90); //(75, 0, 130) // coral color 255, 127, 80
    text("RRT (ms)", pxInf - 80, pyInf - 8 - y_mapped_max - 10); // if RRT_type_1
    //text("RRT [-] (x 1,000)", pxInf - 50, pyInf -8 - y_mapped_max - 16); // if RRT_type_2
    textSize(12);
    text(str(nfs(occupancyMaxRange, 4, 1)), pxInf - 80, pyInf -8 - y_mapped_max + 5);
    pop();
  }
  if (spStatus == 'true'){
    push();
    textFont("Helvetica"); //Courier or Helvetica or Calibri
    textSize(14);
    fill(0, 90, 90); //(75, 0, 130) // coral color 255, 127, 80
    text("timer SP (ms)", pxInf - 89, pyInf - 8 - y_mapped_max - 10); // if RRT_type_1
    //text("RRT [-] (x 1,000)", pxInf - 50, pyInf -8 - y_mapped_max - 16); // if RRT_type_2
    textSize(12);
    text(str(nfs(occupancyMaxRange, 4, 1)), pxInf - 89, pyInf -8 - y_mapped_max + 5);
    pop();
  }
  if (fcStatus == 'true'){
    push();
    textFont("Helvetica"); //Courier or Helvetica or Calibri
    textSize(14);
    fill(0, 90, 90); //(75, 0, 130) // coral color 255, 127, 80
    text("NFC", pxInf - 80, pyInf - 8 - y_mapped_max - 10); // if RRT_type_1
    //text("RRT [-] (x 1,000)", pxInf - 50, pyInf -8 - y_mapped_max - 16); // if RRT_type_2
    textSize(12);
    text(str(nfs(occupancyMaxRange, 4, 1)), pxInf - 80, pyInf -8 - y_mapped_max + 5);
    pop();
  }
} // end function zoom occupancy

//--end function section---------------------------------------------------

function setup(){
  canvas_settings = createCanvas(4900, 2050); // 1320, 750  1560,2050 (52850) for collagen
  canvas_settings.parent("#canvasRiboSeq");
  frameRate(2);

  buttonLSsave = createButton("save local storage to file");
  buttonLSsave.parent("myLocalStorageContainer");
  buttonLSsave.size(180, 27);
  buttonLSsave.style('font-size', '14px');
  buttonLSsave.style('text-align', 'center');
  buttonLSsave.position(738, 8, 'relative');
  buttonLSsave.style('background-color', defaultColor);
  // call back save Local Storage if button is pressed. Take the filename by default.
  buttonLSsave.mouseClicked(saveLS);

  // Retrieve a back-up local storage file to update the local storage.
  // Will erase current active local storage first !
  // Events upon drop file zone are defined in a call back function called 'retrieveLS'.
  dropzone4LS = select('#dropzone2LSretrieve');
  dropzone4LS.dragOver(highlight);
  dropzone4LS.dragLeave(unhighlight);
  dropzone4LS.drop(retrieveLS, unhighlight);

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
  // retrieve the (27) stored transcript ids
  // retrieve the (27) stored transcript ids
  iMax = getItem("numbLines");
  for (i=0; i<iMax;i++){
    keyOne = getItem(str(i));
    // retrieve the readtagID:
    geneTagID[i] = keyOne;
    // remove the parentheses () around the geneID:
    StoredGeneID[i] = keyOne.substring(1, keyOne.length-1);
    // retrieve the copy number for each transcript:
    if (getItem("nonUniformmRNAabundanceStatus")=='true'){
      StoredTransCount[i] = getItem('readCopiesOf_'+str(i));
    }
    else{
      StoredTransCount[i] = 1;
    }
    //retrieve the mRNA sequence:
    myRNA[i] = getItem(keyOne);
  }
  // total number of transcripts (including copies):
  transcriptsFullNb = 0;
  transcriptsLengthInCodons = 0;
  for (i=0; i<iMax; i++){
    transcriptsFullNb += int(StoredTransCount[i]);
    for (i_read = 0; i_read<StoredTransCount[i]; i_read++){
      // compute the transcriptome length:
      transcriptsLengthInCodons += int(myRNA[i].length/3);
    }
  }

  // get proteins counts info:
  for (i=0; i<i_length; i++){// i is time index here
    StoredProtCount[i] = []; //remember this is a 2D array instantiation
    if (getItem('timePoint#'+str(i)) != null){
      StoredTimes[i] = getItem('timePoint#'+str(i));
      // for each of the 36 time points:
      // retrieve the protein counts from local storage:
      StoredProtCount[i] = getItem("protCountTimePoint#"+str(i));
      // if the retrieved field is recognized as a string, convert it into an array of integer numbers!
      // [178,163,187,56,40,52,149,78,196,94,60,232,65,45,60,47,53,61,573,51,48,52,55,78,47,57,291]
      // convert this string into an array of integers...
      if (typeof StoredProtCount[i] == "string"){
        // remove both brackets:
        string_array = StoredProtCount[i].substring(1, StoredProtCount[i].length - 1);
        StoredProtCount[i] = [];
        arrayOfStrings = string_array.split(/,/);
        for (ind=0; ind < arrayOfStrings.length; ind++){// ind is transcript number index here
          StoredProtCount[i][ind] = int(arrayOfStrings[ind]);
        }
      }
    }
  }

  // get RiboSeq info:
  riboRRT_string = getItem('RiboSeqRRT');
  riboSP_string = getItem('RiboSeqSP');
  if(getItem('riboSeqFC') != null){
    riboFC_string = getItem('riboSeqFC');
    snapshots_nb = getItem('riboSeqSnapShotsNb');
    footprintedFragmentsCount = getItem('footprintFragNb');
    console.log('riboFC_string[7]');
    console.log(riboFC_string[7]);
    console.log('riboFC_string[7][0]',riboFC_string[7][0]);
    console.log('snapshots_nb', snapshots_nb);
    console.log('footprintedFragmentsCount', footprintedFragmentsCount);
  }
  // (and check for format to avoid string/array of float problems...)
  //console.log('riboRRT[tr][codon] is a string?', typeof riboRRT_string == "string");
  if (typeof riboRRT_string == "string"){
    // [[3.0, 2.1, 1.1],[4.0, 5.5, 6.8],[1.2, 2.22, 3.3]]
    var stringArray1D = riboRRT_string.substring(2, riboRRT_string.length-2);
    var array1DofStrings = stringArray1D.split(/],\[/);
    // we have 58 comma separated values of # codons values separated by commas...
    //console.log('array1DofStrings.length RRT', array1DofStrings.length);
    for (tr=0; tr<array1DofStrings.length; tr++){
      //console.log('array1DofStrings[tr]', array1DofStrings[tr]);
      riboRRT[tr] = []; // second dimension of array initiated.
      var arrayOfRRTonCodon = array1DofStrings[tr].split(/,/)
      for (i_codon=0; i_codon<arrayOfRRTonCodon.length; i_codon++){
        var rrt_float = float(arrayOfRRTonCodon[i_codon]);
        riboRRT[tr].push(rrt_float);
      }
    }
  }

  //console.log('riboSP[tr][codon] is a string?', typeof riboSP_string == "string");
  if (typeof riboSP_string == "string"){
    // [[3.0, 2.1, 1.1],[4.0, 5.5, 6.8],[1.2, 2.22, 3.3]]
    var stringArray1D = riboSP_string.substring(2, riboSP_string.length-2);
    var array1DofStrings = stringArray1D.split(/],\[/);
    // we have 58 comma separated values of # codons values separated by commas...
    //console.log('array1DofStrings.length SP', array1DofStrings.length);
    for (tr=0; tr<array1DofStrings.length; tr++){
      //console.log('array1DofStrings[tr] SP', array1DofStrings[tr]);
      riboSP[tr] = []; // second dimension of array initiated.
      var arrayOfSPonCodon = array1DofStrings[tr].split(/,/)
      for (i_codon=0; i_codon<arrayOfSPonCodon.length; i_codon++){
        var sp_float = float(arrayOfSPonCodon[i_codon]);
        riboSP[tr].push(sp_float);
      }
      //console.log('floated riboSP[tr]', tr, riboSP[tr]);
    }
  }
  //console.log('riboFC[tr][codon] is a string?', typeof riboFC_string == "string");
  console.log('true / false ?', typeof riboFC_string == "string");
  if (typeof riboFC_string == "string" && riboFC_string != ''){
    // [[3.0, 2.1, 1.1],[4.0, 5.5, 6.8],[1.2, 2.22, 3.3]]
    var stringArray1D = riboFC_string.substring(2, riboFC_string.length-2);
    var array1DofStrings = stringArray1D.split(/],\[/);
    // we have 58 comma separated values of # codons values separated by commas...
    console.log('array1DofStrings.length FC', array1DofStrings.length);
    for (tr=0; tr<array1DofStrings.length; tr++){
      //console.log('array1DofStrings[tr]', array1DofStrings[tr]);
      riboFC[tr] = []; // second dimension of array initiated.
      var arrayOfFConCodon = array1DofStrings[tr].split(/,/)
      for (i_codon=0; i_codon<arrayOfFConCodon.length; i_codon++){
        var fc_float = float(arrayOfFConCodon[i_codon]);
        riboFC[tr].push(fc_float);
      }
    }
  }
  else{
    console.log('bro is on the jesus motorbike');
    console.log('transcriptsFullNb', transcriptsFullNb);
    for (tr=0; tr < transcriptsFullNb; tr++){
      riboFC[tr] = []; // second dimension of array initiated.
      console.log('riboFC_string[tr].length', riboFC_string[tr]);
      for (i_codon=0; i_codon < riboFC_string[tr].length; i_codon++){
        console.log('bro is on the jesus motorbike again and again');
        riboFC[tr][i_codon] = riboFC_string[tr][i_codon];
      }
    }
  }
  console.log('we collected riboFC[7][*]:', riboFC[7]);
  console.log('we collected riboFC[7][0]:', riboFC[7][0]);

  // Create Input box for the gene id tag name:
  tagInput = createInput(tagInputName); //tagInputName as argument
  tagInput.parent("#myButtonContainer");
  //tagInput.size(60);
  tagInput.position(825, 8, 'relative');
  tagInput.style('width', '90px');
  tagInput.style('height', '22px');
  tagInput.style('background-color', '#B4B4B4');
  tagInput.style('color', 'black');
  tagInput.style('fontSize', '20px');
  tagInput.style('text-align', 'center');
  //tagInput.input(getTagNameID);
  tagInput.changed(getTagNameID);

  // Create paragraph (for input box context) to set gene/transcript tag name (without parenthesis):
  //TagNameIDstring = StoredGeneID[0];
  tagNameBox = createP(TagNameIDstring);// is the default value;
  tagNameBox.parent("myButtonContainer");
  tagNameBox.position(985, -13, 'relative');
  tagNameBox.style('fontSize', '20px');
  tagNameBox.style('textFont', 'Helvetica');
  tagNameBox.style('color', 'white');

  // creation of buttons to select the type of RiboSeq profile to display
  // Ribo-Seq of type RRT (actual ribosome residence time during simulation)
  buttonRRT = createButton("actual RRT during simulation");
  buttonRRT.parent("myButtonContainer");
  buttonRRT.size(200, 27);
  buttonRRT.style('font-size', '14px');
  buttonRRT.style('text-align', 'center');
  buttonRRT.position(825, 52, 'relative');
  buttonRRT.style('background-color', 'rgb(0, 218, 0)');//'rgb(0, 255, 0, 255)'

  // Ribo-Seq of type RRT (actual ribosome residence time during simulation)
  buttonSP = createButton("queueing time set points");
  buttonSP.parent("myButtonContainer");
  buttonSP.size(200, 27);
  buttonSP.style('font-size', '14px');
  buttonSP.style('text-align', 'center');
  buttonSP.position(1045, 52, 'relative');
  buttonSP.style('background-color', ButtonDefaultColor);

  // Ribo-Seq of type RRT (actual ribosome residence time during simulation)
  buttonFC = createButton("RPF count sampling");
  buttonFC.parent("myButtonContainer");
  buttonFC.size(200, 27);
  buttonFC.style('font-size', '14px');
  buttonFC.style('text-align', 'center');
  buttonFC.position(1265, 52, 'relative');
  buttonFC.style('background-color', ButtonDefaultColor);

  // Create slider for ratio of ribosomes cardinality over transcripts cardinality (included their copy number):
  rangeSlider = createSlider(10, 130, 130, 10);
  rangeSlider.parent("myButtonContainer");
  rangeSlider.position(1055, 22);
  rangeSlider.style('width', '390px');
  // Create paragraph with codon window range displayed:
  rangeDisplay = createP("codon range window: ");
  rangeDisplay.parent("myButtonContainer");
  rangeDisplay.position(1125, -16);
  rangeDisplay.style('fontSize', '16px');
  rangeDisplay.style('textFont', 'Helvetica');
  rangeDisplay.style('color', 'blue');

}// end setup

function draw(){// draw starts here
  // display the current value of the clock:
  background(255);
  //fill(0, 90, 90);
  //strokeWeight(1);
  //stroke(235);//235
  //textSize(24);//24
  //textAlign(LEFT);
  //textFont("Helvetica"); //Helvetica, Georgia, Courrier new ?
  //push();
  var time_running = millis();
  var time_minutes = floor(time_running/(1000*60));
  var time_seconds = int((time_running/1000) % 60);
  var timer_text = str(nfs(time_minutes, 3, 0))+" min:"+str(nfs(time_seconds, 2, 0))+" sec.";
  //text("Time running:"+str(nfs(time_minutes, 3, 0))+" min:"+str(nfs(time_seconds, 2, 0))+" sec.", 7, 20);
  //text("Time running:"+ str(timer_text), 7, 17);
  //pop();

  // get status of buttons for Ribosome occupancy mode of display:
  buttonRRT.mousePressed(buttonRRTOnOff);
  // test the status of all the buttons:
  rrtStatus = str(buttonRRT.style('background-color') == 'rgb(0, 218, 0)');

  buttonSP.mousePressed(buttonSPOnOff);
  spStatus = str(buttonSP.style('background-color') == 'rgb(0, 218, 0)');

  buttonFC.mousePressed(buttonFCOnOff);
  fcStatus = str(buttonFC.style('background-color') == ButtonColorGreen);

  //console.log('statuses=', rrtStatus, spStatus, fcStatus);

  // get the geneID tag for which you want to display the Ribo-Seq profile:

  console.log('geneIndex', geneIndex);
  geneTagString = "(" + TagNameIDstring + ")";
  console.log("First draw gene tag=", geneTagString);
  TagNameIDstring = geneTagString.substring(1, geneTagString.length-1);
  console.log("TagNameIDstring=", TagNameIDstring);
  console.log('geneTagID list?', geneTagID);
  console.log('index in geneTagID of the geneTagString', geneTagID.indexOf(geneTagString))
  console.log('geneIndex in draw', geneIndex);

  // retrieve the number of reads for this geneTag:
  readNumber = StoredTransCount[geneIndex];

  // get the set of read indices for this geneTag:
  read_indices = [];
  for (var i = 0; i < readNumber; i++){
    readTagString = geneTagID[geneIndex] + "#" + str(i);
    read_indices.push(int(indexOfRead(StoredGeneID, StoredTransCount, readTagString))); // search in the readlist!
    }
  // get the mRNA sequence of this geneTag:
  mRNAseq = myRNA[geneIndex];
  // build a codon numbered list from 5' to 3' mRNA:
  codonSeqList = [];
  aaSeqList = [];
  codonOrderedNumb = [];
  for (let i = 0; i < mRNAseq.length/3; i++){
    codonSeqList.push(mRNAseq.substring(0+i*3, (i+1)*3));
    codonOrderedNumb.push(i);
    aaSeqList.push(decode(mRNAseq.substring(0+i*3, (i+1)*3)));
  }
  // get RiboSeq info:
  // sum all RRT, (resp. SP) over the read copies of this geneID tag for each codon belonging to this transcript:
  RRT = [];
  SP = [];
  FC = [];
  RRT_type_1 = [];
  RRT_type_2 = [];
  SP_type_1 = [];
  SP_type_2 = [];
  //FC_type_1 = [];
  //FC_type_2 = [];

  for (i_codon = 0; i_codon < mRNAseq.length/3; i_codon++){// number of codons in this transcript
    var RRT_cumul = 0;
    var SP_cumul = 0;
    var FC_cumul = 0;
    for (i = 0; i < readNumber; i++){
      RRT_cumul += riboRRT[read_indices[i]][i_codon];
      SP_cumul += riboSP[read_indices[i]][i_codon];
      FC_cumul += riboFC[read_indices[i]][i_codon];
    }
    // for FC_cumul we divide by the number of reads to get a normalized footprint count:
    FC.push(FC_cumul/readNumber);
    // for RRT or SP, do not divide by the number of reads because we will divide by the number of ribosomes!
    RRT.push(RRT_cumul);
    SP.push(SP_cumul);
    // How many ribosomes have been through this transcript?
    // Answer: the protein count for this transcript at last time point:
    // if the count exists at the final time (otherwise take one):
    if (getItem('timePoint#'+str(i_length - 1)) != null){
      var riboThroughCount = StoredProtCount[i_length - 1][geneIndex];
    }
    else {
      var riboThroughCount = 1;
    }
    // So, if RRT (resp. SP) is divided by the total number of ribosomes that have translated this transcripts by the end
    // of the simulation, you get, 'RRT_type_1', the total time spent by one ribosome on average on each codon
    // of this transcript geneID.
    // The unit is a time unit in ms.
    RRT_type_1.push(RRT_cumul/riboThroughCount);
    SP_type_1.push(SP_cumul/riboThroughCount);
  }// end for

  // optionnaly, if you divide the previous result RRT_type_1 by the total time spent
  // on this transcript by a ribosome on average, you have RRT_type_2 the relative fraction
  // of time a ribosome has spent on average on a codon of this transcript. The unit is adimensional.
  // Adimensinal but scaled to per thousand (pour mille).
  var totalTimeOnTranscript = 0;
  for (i_codon = 0; i_codon < RRT_type_1.length; i_codon++){
    totalTimeOnTranscript += RRT_type_1[i_codon];
  }
  for (i_codon = 0; i_codon < RRT_type_1.length; i_codon++){
    RRT_type_2.push(RRT_type_1[i_codon] * 1000/totalTimeOnTranscript);
  }
  // likewise for SP:
  var totalSPOnTranscript = 0;
  for (i_codon = 0; i_codon < SP_type_1.length; i_codon++){
    totalSPOnTranscript += SP_type_1[i_codon];
  }
  for (i_codon = 0; i_codon < SP_type_1.length; i_codon++){
    SP_type_2.push(SP_type_1[i_codon] * 1000/totalSPOnTranscript);
  }

  //******************************************************************************************************
  //  D I S P L A Y   R I B O - S E Q   profile with RRT per window of window_range = 130 codons (390 nts)
  //******************************************************************************************************

  background(255); //235

  // Create paragraph with selected slider value for codon window range displayed:
  rangeDisplayed = createP(str(nfs(window_range, 3, 0)));
  rangeDisplayed.parent("myButtonContainer");
  rangeDisplayed.position(1295, -16);
  rangeDisplayed.style('fontSize', '16px');
  rangeDisplayed.style('textFont', 'Helvetica');
  rangeDisplayed.style('color', 'blue');
  rangeDisplayed.style('background', '#CD7D6B');

  // graphic box title
  var lineOne = "Predicted ribosome occupancy map (P-site reference frame)";
  push();
  textFont("Helvetica"); //Courier or Helvetica or Calibri
  textSize(24);
  fill(0, 90, 90); //(75, 0, 130) // coral color 255, 127, 80 // teal 0, 90, 90
  textAlign(LEFT);
  text(lineOne, pxInf - 12, pySup - 20);
  pop();

  // Display line of transcript with scale identical to settings tab (1 codon = 3px, 1 nt=1px):
  push();
  strokeWeight(4);
  strokeCap(ROUND);
  stroke(248, 196, 211); // rgb code for fast codon
  line(pxInf-12, pyInfMap -40, pxInf-12 + mRNAseq.length, pyInfMap - 40); // (1 px per nt)
  pop();

  // add name of the transcript to the left as text:
  fill(0);
  stroke(235);
  textSize(18);//24
  textAlign(LEFT);
  push();
  text(geneName, pxInf - 89, pyInfMap - 35);
  pop();
  // add length of transcript in codons (and nts):
  fill(0);
  stroke(235);
  textSize(12);//24
  textAlign(CENTER);
  push();
  text(nfs(mRNAseq.length/3, 3, 0), pxInf - 60, pyInfMap - 22);
  text("("+str(nfs(mRNAseq.length, 4, 0))+" nts)", pxInf - 60, pyInfMap - 10);
  pop();

  // get the window range from the slider value:
  window_range = rangeSlider.value();

  // Depending on the status of the selected (pressed) button: display the chosen occupancy map:
  if (rrtStatus == 'true'){
    // display RRT as measured on each codon during simulation and aggregated over
    // RNA-seq copy numbers and avergaed by the number of ribosomes that have terminated this trnacript:
    // call the OccupancyMapShow function(array, pixMaxInf_y, pixMaxSup_y)
    occupancyMapShow(RRT_type_1, pyInfMap - 40, pySup);
  }
  if (spStatus == 'true'){
    // display timer set points as attributed on each codon during simulation and aggregated over
    // RNA-seq copy numbers and avergaed by the number of ribosomes that have terminated this transcript:
    // call the OccupancyMapShow function(array, pixMaxInf_y, pixMaxSup_y)
    occupancyMapShow(SP_type_1, pyInfMap - 40, pySup);
  }
  if (fcStatus == 'true'){
    // display footprint count as snapshot during simulation and aggregated over
    // RNA-seq copy numbers:
    // call the OccupancyMapShow function(array, pixMaxInf_y, pixMaxSup_y)
    occupancyMapShow(FC, pyInfMap - 40, pySup);
  }

  if (mouseX >= pxInf-12 && mouseX <= pxInf-12 + mRNAseq.length && mouseY <= pyInfMap -40 && mouseY >= pySup){
    console.log('mouseX', mouseX, 'mouseY', mouseY);
    if (mouseIsPressed == true && mouseButton == LEFT){
      lower_px = mouseX;
    }
    // draw a dashed coral vertical line at lower and upper bound of 130 codons wide window:
    push();
    strokeWeight(3);
    stroke(255, 127, 80); // coral  // stroke(100, 0, 170) stroke(150, 125);
    drawingContext.setLineDash([8, 16]);
    line(lower_px, pyInfMap - 25, lower_px, pySup);
    line(lower_px + window_range*3, pyInfMap - 25, lower_px + window_range*3, pySup);
    pop();

    lower_codon = int(map(lower_px, pxInf-12, pxInf-12 + mRNAseq.length, 0, mRNAseq.length/3));
  }

  //lower_codon = int(map(lower_px, pxInf-12, pxInf-12 + mRNAseq.length, 0, mRNAseq.length/3));

  // change lower or start codon when arrow is pressed...
  if (keyIsPressed){
    if (keyCode == 39){// 39 is keyCode for RIGHT_ARROW
      lower_codon += window_range;
    }
    if (keyCode == 37){// 37 is keyCode for LEFT_ARROW
      lower_codon -= window_range;
    }
    lower_codon = constrain(lower_codon, 0, mRNAseq.length/3 - window_range);
    lower_px = map(lower_codon, 0, mRNAseq.length/3, pxInf-12, pxInf-12 + mRNAseq.length);
    push();
    strokeWeight(3);
    stroke(255, 127, 80); // coral  // stroke(100, 0, 170) stroke(150, 125);
    drawingContext.setLineDash([8, 16]);
    line(lower_px, pyInfMap - 25, lower_px, pySup);
    line(lower_px + window_range*3, pyInfMap - 25, lower_px + window_range*3, pySup);
    pop();
  }

  // The code below will be called in a function whose argument will be the chosen codon starting position
  // The function will be a callback function of some action on the mouse...or on arrow keys pressed.
  // Display axis x of codons ordered by sequence; alternate color blue and gray by triplets:

  fill(0, 0, 255, 125);
  stroke(235);
  textSize(16);
  textAlign(CENTER);
  //push();
  // codon:
  text_codon = '';
  text_aa = '';

  start_codon = lower_codon;
  var end_codon = start_codon + window_range;
  if (end_codon > mRNAseq.length/3){
    end_codon = mRNAseq.length/3;
  }
  //for (var i_codon=0;i_codon<mRNAseq.length/3; i_codon++){
  for (var i_codon=start_codon;i_codon<end_codon; i_codon++){
    // codon text
    text_codon = '';
    text_codon = str(codonSeqList[i_codon]);
    text_aa = '';
    text_aa = str(aaSeqList[i_codon]);
    if (i_codon % 2 == 0){ // even color = blue:
      fill(0, 0, 255, 125);
    }
    else {// odd color
      fill(150, 50, 50, 250);
    }
    //push();
    var i_window = i_codon - start_codon;
    text(text_codon, pxInf + i_window * dx_triplet, pyInf + 20);
    // Display axis x by amino acid
    text(text_aa, pxInf + i_window * dx_triplet, pyInf + 40);
    // Display axis x by codon number:
    text(str(i_codon), pxInf + i_window * dx_triplet, pyInf + 60);
    //pop();
  }
  // nucleotide numbering:
  // Display axis x by nucleotide number by multiple of 3   (2)
  fill(10, 10, 10, 125);
  stroke(235);
  textSize(10);
  textAlign(RIGHT);
  push();
  //for (var i_codon=0;i_codon<mRNAseq.length/3; i_codon++){
  for (var i_codon=start_codon;i_codon<end_codon; i_codon++){
    // nt number by 3:
    var i_window = i_codon - start_codon;
    text(str(nfs((i_codon+1)*3, 4, 0)), pxInf + i_window * dx_triplet + 16, pyInf + 5);
  }
  pop();
  // Display gene tag name on the left of codons sequence:
  // add name of the transcript to the left as text:
  fill(0);
  stroke(235);
  textSize(18);//24
  textAlign(LEFT);
  push();
  text(geneName, pxInf - 89, pyInf + 20);
  pop();
  // Display scale and label of y-axis on the left of the zoomed map:
  // depends on the window local maximum:
  // depends on rrt or sp or fc Status:
  // call function zoomOccupancy(start_codon, window_range);
  // Depending on the status of the selected (pressed) button: display the chosen occupancy map:
  push();
  // display the x-axis line as the transcript window:
  strokeWeight(4);
  strokeCap(ROUND);
  stroke(248, 196, 211); // rgb code for fast codon
  line(pxInf-12, pyInf-8, pxInf-12 + window_range * dx_triplet - 12, pyInf-8); // (1 px per nt)
  pop();
  if (rrtStatus == 'true'){
    // display RRT as measured on each codon during simulation and aggregated over
    // RNA-seq copy numbers and avergaed by the number of ribosomes that have terminated this trnacript:
    // call the zoomOccupancy
    zoomOccupancy(RRT_type_1, start_codon, window_range, pyInf - 8, pyInfMap-10);
    }
  if (spStatus == 'true'){
    // display SP as attributed on each codon during simulation and aggregated over
    // RNA-seq copy numbers and avergaed by the number of ribosomes that have terminated this trnacript:
    // call the zoomOccupancy
    zoomOccupancy(SP_type_1, start_codon, window_range, pyInf - 8, pyInfMap-10);
    }
  if (fcStatus == 'true'){
    // display footprint count as snapshot during simulation and aggregated over
    // RNA-seq copy numbers:
    // call the zoomOccupancy
    zoomOccupancy(FC, start_codon, window_range, pyInf - 8, pyInfMap-10);
    // display Ribo-Seq coverage statistics:
    var lineStat_1 = 'Number of ribosome footprint snapshots = '+str(snapshots_nb);
    var lineStat_2 = 'Number of ribosome footprint fragments = '+str(footprintedFragmentsCount);
    var lineStat_3 = 'Averaged number of ribosomes on all transcripts per snapshot = '+str(nfs(footprintedFragmentsCount/snapshots_nb, 3, 1));
    var lineStat_4 = 'Exact Ribo-Seq coverage = '+str(nfs(footprintedFragmentsCount*10/transcriptsLengthInCodons, 2, 1));
    fill(0);
    stroke(235);
    textSize(18);//24
    textAlign(LEFT);
    push();
    text(lineStat_1, pxInf - 89, pyInf + 110);
    text(lineStat_2, pxInf - 89, pyInf + 135);
    text(lineStat_3, pxInf - 89, pyInf + 160);
    text(lineStat_4, pxInf - 89, pyInf + 185);

    pop();
    }
}// end draw
