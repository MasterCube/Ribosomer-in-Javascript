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
var dropzone;
var fastaparag;
var fastatxt ="empty string";
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
var uniqArray= [];
var uniqString;
var nucleotideSeqArray = [];
//var uniqStringBIS = ""; // empty string
var CDSstartIndex=0, CDSstopIndex=10;
var CDSmRNA;
var upperLim;
var idTagsMaxNumb = 16;  // max number of transcripts that will be processed by Ribosomer.
var outFile;
let myTranscriptsTable;
//------------------------------------------------------------------------------
// functions section:
//------------------------------------------------------------------------------

function highlight(){
  dropzone.style("background-color", "grey");
}
function unhighlight(){
  dropzone.style('background-color', '#C8C4');//#C8C4
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
    fastaparag.parent('#fastaHere');

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
    print('output22 length:', output22.length);
    print('output22:', output22);
    print('output11:', output11);

    result = createP('Results zone for all ORF: ');
    result.parent('#resultsZone');
    // the following paragraph does not fit well in the html page. Check for css design !
    //createP('The number of retrieved gene id tags is ' + str(output22.length));

    for (i=0; i<output22.length;i++){
      createP("Gene id tag retrieved:" + output22[i] + "- description: " + output11[i]);
    }

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
      console.log('captured groups with exec: ', capturedGroups);
      capturedmRNA= capturedGroups[1]; // possibly contains whate spece...
      lastCodon.push(capturedGroups[2]);
      console.log('captured mRNA=', capturedmRNA);
      numNucleotides = capturedmRNA.length; // possibly contains white spaces...
      console.log('length of mRNA:', numNucleotides);
      // remove the white spaces in each particular captured group 1:
      capturedmRNAcorrect = capturedmRNA.match(r_ATCG);
      console.log('capturedmRNAcorrect is ', capturedmRNAcorrect);
      ORFmRNA = "";
      for (i=0; i<capturedmRNAcorrect.length; i++){
        ORFmRNA += capturedmRNAcorrect[i];
      }
      console.log(ORFmRNA);
      console.log('correct length:', ORFmRNA.length);
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
      ORFflag.push(flag);

      capturedGroups = reg_mRNA.exec(uniqString);
    }

    console.log('nucleotide seq array: ', nucleotideSeqArray);
    console.log('all correct mRNAS:', myRNAs);

    // ----------- save the relevant setting results:   ---------------------------
    // // - the transcripts log:
    myTranscriptsTable = new p5.Table();
    myTranscriptsTable.addColumn('geneID');
    myTranscriptsTable.addColumn('description');
    myTranscriptsTable.addColumn('orf');
    myTranscriptsTable.addColumn('lastCodon');
    myTranscriptsTable.addColumn('mRNAlength');

    upperLim = min(idTagsMaxNumb, output22.length);
    for (i=0; i<upperLim; i++){
      let newRow = myTranscriptsTable.addRow();
      newRow.setString('geneID', output22[i].substring(1, output22[i].length-1).toUpperCase());
      newRow.setString('description', output11[i].substring(1, output11[i].length));
      newRow.setString('orf', str(ORFflag[i]));
      newRow.setString('lastCodon', lastCodon[i]);
      newRow.setNum('mRNAlength', myRNAs[i].length);
    }
    saveTable(myTranscriptsTable, "myTranscriptsLog.tsv", "tsv");
    // This is saved in a file that will find itself in the downloads directory of the client,
    // not in datadase as service (server side), like FIREBASE (google).

    // This alternative way of saving the table or a file enables the browser to prompt you for authorization...
    // // C:\Users\marcjoiret\Desktop\ONGOING\j5sketch RIBOSOMER\output
    //outFile = createWriter(C:\Users\marcjoiret\Dropbox\PC\Desktop\ONGOING\j5sketch RIBOSOMER\output\"myTranscripts.tsv", "tsv");
    outFile = createWriter("myTranscripts.tsv", "tsv");
    // // write the header line:
    outFile.print("geneID\tDescription\tORF\tlastCodon\tmRNAlength");
    upperLim = min(idTagsMaxNumb, output22.length);
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
    //for (k=0; k<mRNAs.length; k++){
    //  s=mRNAs[k].substring(mRNAs[k].indexOf('ATG'), mRNAs[k].length);
    //  myarray = splitTokens(s,"''' ''\s\s''\s''\t''\t\t''\n''\n\n'");
    //  ss = join(myarray,'');
    //  createP('transcript number ' + str(k+1) + ' nt length= '+ str(ss.length) + ' last codon ' + str(ss.charAt(ss.length-4))+str(ss.charAt(ss.length-3))+str(ss.charAt(ss.length-2)));
    //  createP(ss);
    //  }
    //  createP('length of single string: '+ str(uniqString.length) + ' and is computed as: ' + str(CDSmRNAlength) + ' or '+ str(CDSmRNAlengthCorrect));
    //  createP('length of first mRNA:'+ str(CDSmRNA.length));
    //  createP('first mRNA:'+'<br>'+CDSmRNA);
    //  createP('uniqString: ' + '<br>' + uniqString);

    //----- correct mRNAs: mRNAseq:
    //createP('CORRECT mRNA TRANSCRIPTS LIST:'+'<br>'+ '+++++++++++++++++++++++++++++++')

    // for (k = 0; k<mRNAseq.length; k++){
    //   createP('transcript number ' + str(k+1) + 'nt length= '+ str(mRNAseq[k].length)+' last codon: '+mRNAseq[k].substring(mRNAseq[k].length-3, mRNAseq[k].length));
    //   createP('transcript code of the mRNA:'+'<br>'+ mRNAseq[k]);
    // }

    textProcessed = "STATUS: Your first "+ str(idTagsMaxNumb) + " transcripts have been processed and saved in a temporary table to be used in the Ribosomer bench:"

  }

//------------------------------------------------------------------------------
// setup
//------------------------------------------------------------------------------
function setup(){
// noCanvas();
canvas_settings = createCanvas(1320, 400);
canvas_settings.parent("#canvas_settings_container");
dropzone = select('#dropzone4file');
dropzone.dragOver(highlight);
dropzone.dragLeave(unhighlight);
dropzone.drop(getFile, unhighlight);

// connect to FIREBASE as a database service:
// { initializeApp } = require('firebase-admin/app
// Your web app's Firebase configuration
// For Firebase JS SDK v7.20.0 and later, measurementId is optional
// const instead of var:
const firebaseConfig = {
  apiKey: "AIzaSyDXG3A9LZNt0YIDXpqCphwwWHfXsvcEgLU",
  authDomain: "ribosomedigitaltwindb.firebaseapp.com",
  projectId: "ribosomedigitaltwindb",
  storageBucket: "ribosomedigitaltwindb.appspot.com",
  messagingSenderId: "433909121779",
  appId: "1:433909121779:web:476c181d4f4380ea9d6d86",
  measurementId: "G-4NF3167Q49"
};

// Initialize Firebase
const app = initializeApp(firebaseConfig);
const analytics = getAnalytics(app);
//firebase.initializeApp(firebaseConfig);

// provide data to the firebase database:
var database = app.database();
var ref = database.ref('transcripts');
var data = {
  geneID: "COL1",
  mRNAseq: "UUCGAAAAUUU"
}
ref.push(data);
} // end setup()
//------------------------------------------------------------------------------
// draw
//------------------------------------------------------------------------------
function draw(){
  background(0, 90, 90); //Gainsboro 220, 220, 220 or teal green
  textSize(18);
  fill(255);
  textFont("Helvetica");
  text(textProcessed, 10, 25);
  textSize(14);


  // provide gene id tags, length of peptide, name of gene

  if (output22[0] != null){
    text('Here is the list of the first '+ str(idTagsMaxNumb)+ ' out of the '+ str(output22.length) + ' provided transcripts variants (capitalized letters gene ID tags): ', 30, 50);
    text('ORF', 927, 50);
    text('last codon', 989, 50);
    text('peptide length', 1064, 50);
    text('mRNA length', 1175, 50);
    upperLim = min(idTagsMaxNumb, output22.length);
    for(i=0; i<upperLim; i++){
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
    }
  }
}

//-------END-------------------------------------------------------------------
