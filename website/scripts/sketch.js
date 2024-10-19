// random modules:



//let tRNA1pic;
// uncomment this to get a nice .png to draw a real tRNA
function preload() {
  //tRNA1pic = loadImage("http://localhost:8080/images/tRNAnormal.png");
  tRNA1pic = loadImage("http://localhost:3000/tRNAnormal.png");
  //tRNA1pic = loadImage("http://localhost:8080/tRNAnormal.png");
  //tRNA1pic = loadImage('./images/tRNAnormal.png');
}

var picWidth = 34;  // 1700 x 790 px
var picHeight = 16;

// Objects instance defined:
// tRNAs instances:
let tRNAiso1_K;     // isoacceptor 1 tRNA for lysine
let tRNAiso2_K;     // isoacceptor 2 tRNA for Lysine
var tRNAlinearMaxAmplitude = 50; // max pixels amplitude per frame
var HALF_PI = 1.57;
var tRNAangularSpeed = HALF_PI/9.0; // 10 degrees per frame

// Read the file with the ORFs:
var dropzone;
var tRNA;
var txt;
//var output;
//---------------------------------------------------------------------------
// genetic code 3D-array: 4 by 4 by 4 3D array with 64 elements:
var GenetiCode = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]]
//------------------------------------------------------------------------
var sequence = "AUGAAAAACAAGAAUACAACCACGACUAGAAGCAGGAGUAUAAUCUAA";
var transcriptLength;
var peptide = [];
//------------------------------------------------------------------------
// Translation rates for ribosomes:
// Ribosome Elongation Rate:
var RiboPosX = 415;
var x = RiboPosX;
var RiboElongationRate = 6; // use 6:   6 codons/s = 18px/s
function setup() {
  // put setup code here
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
  // display tRNA picture as a dom element in the header1:
  // tRNA = createP();
  // tRNA.position(1100, 90);
  // tRNA.parent('tRNAHere');
  // var img = createImg('images/tRNAnormal.png');
  // img.size(272, 116);
  // img.parent(tRNA);

  //noCanvas();
  var canvasOne =  createCanvas(1920, 938);
  canvasOne.parent('#canvasHere');
  //------------------------------------------------------------------------
  frameRate(12); // use 18 frames/s
  //------------------------------------------------------------------------


  // Read the file:
  // from directory accessed by the server:
  //loadStrings('Data/CDSfasta.txt', fileready);// callback function: fileready
  // or from a drop zone in the html page index:
  /*
  var button = select("#loadfile");
  button.mousePressed(loadFile);
  */
  // The user will chose his own file by clicking on a button in the html page:
  //createFileInput(fileSelected);
  //output = select("#output");

  // select file in the dropzone
  //dropzone = select("filedropzone");
  // dropzone = select(<a href="settings.html#filedropzone">);
  // dropzone.mouseOver(highlight);
  // dropzone.mouseOut(unhighlight);


  // tRNAs instantiations:
  tRNAiso1_K = new Trna();
  tRNAiso2_K = new Trna();

  console.log('tRNAiso1_K : ', tRNAiso1_K);
  // this generate a gamma distributed random value:
  alpha=1.6438;
  beta=108.28;
  console.log('random gamma : ', jStat.gamma.sample(alpha, beta));
  console.log('gamma distribution mean : ', jStat.gamma.mean(alpha, beta));
  console.log('gamma distribution variance : ', jStat.gamma.variance(alpha, beta));
  console.log('pdf for gamma of x=100 :', jStat.gamma.pdf(100, alpha, beta));

  // Standard genetic code defined once and for all as a 3D array:
  // A, C, G, U maps 0, 1, 2, 3 indexes: GenetiCode[0][0][0] = K Lysine decoded by AAA...
  GenetiCode[0][0][0] = 'K'; // lysine Lys (AAA)
  GenetiCode[0][0][1] = 'N'; // asparagine Asn (AAC)
  GenetiCode[0][0][2] = 'K'; // lysine Lys (AAG)
  GenetiCode[0][0][3] = 'I'; // isoleucine Ile (AAU)
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
}
function decode(codon){
  // This function receives a codon in argument as a 3 letter string.
  // This function returns the corresponding decoded amino acid
  // (called 'residue') using the standard genetic code.
  var i, j, k;
  var FirstLetter, SecondLetter, ThirdLetter;
  var alphabet = "ACGU";
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
function highlight(){
  dropzone.style('background-color', '#C12C6');//#787878 or rgb(205, 160, 135, 0.5)
}
function unhighlight(){
  dropzone.style('background-color', '#BC4F4F');
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
  console.log(file.data);
}

function draw() {
  // put drawing code here
  //console.log(txt);
  //output.html(txt);
  background(235); //224

  // tRNAs jittering around:
  tRNAiso1_K.jitter();
  tRNAiso1_K.show();
  tRNAiso2_K.jitter();
  tRNAiso2_K.show();


  // mRNA transcripts:
  // The length of the transcript is proportional to number of codons
  // 1 codon = 3 pixels
  strokeWeight(4);
  strokeCap(ROUND);
  stroke(0, 255, 0);
  line(400, 256, 1300, 256); // 300 codons
  stroke(255, 0, 0);
  line(400, 512, 1150, 512); // 250 codons
  stroke(0, 0, 255);
  line(400, 768, 1840, 768); // 480 codons
  // Ribosome design variant 1:
  //stroke(24);
  //line(300, 256, 330, 256)
  //strokeWeight(2);
  //fill(24, 24, 24, 150);
  //ellipse(315, 231, 50, 50)
  //arc(315, 256, 30, 30, TWO_PI, PI)
  // ------ Ribosome design variant 2:----------------
  strokeWeight(5); // mRNA tunnel in ribosome
  stroke(24, 120, 54, 180);
  line(300, 297, 330, 297);
  strokeWeight(1);
  stroke(24);
  line(300, 297, 330, 297);
  strokeWeight(4);
  stroke(24, 120, 54, 180);
  noStroke();
  // ---------- Large 60S subunit: -----------
  fill(24, 120, 54, 180);
  ellipse(315, 275, 50, 40);
  // ---------- Small 40S subunit: -----------
  rect(308, 300, 14, 16);
  arc(308, 308, 16, 16, HALF_PI, PI+HALF_PI);
  arc(322, 308, 16, 16, PI+HALF_PI, HALF_PI);
  // -----------protein exit tunnel:----------
  stroke(205, 25, 25, 180);
  line(295, 266, 315, 282);
  // ---- ELONGATION MOTION ALONG Transcript: ------------
  // motion condition:
  if (x <= 1300){
    x += RiboElongationRate;
  }
  else {
    x = 200;
  }
  // ---- RIBOSOME OBJECT DESIGN: ------------
  // define line center as ribosome center
  RiboX = x;
  RiboY = 256;
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
  //text("Is this working like hell?", 30, 200);
  //text(GenetiCode[0][0][0], 40, 40);
  //text(GenetiCode[0][0][1], 40, 80);
  //text(GenetiCode[0][0][2], 40, 120);
  //text(GenetiCode[0][0][3], 40, 160);
  transcriptLength = sequence.length;
  var triplet;
  for (var n = 0; n < transcriptLength/3; n += 1){
    triplet = sequence.substring(3*n, 3*n+3);
    peptide[n] = decode(triplet);
  }
  console.log(sequence);
  console.log(peptide);
} // end function draw



// to lauch this sketch: type the 3 keys Ctrl + Alt + L

// objects classes:
// class Trna:
class Trna {
  constructor() {
    // has an image (attribute)

    // has a position
    this.x = 150;
    this.y = 300;
    // has an max angle of oscillation
    this.angle = 0.0;
    }  // end constructor

    // has functionalities (methods) - they don't need to be declared with syntax 'function'
    // because these "functions" are encapsulated in a class and are thus now "methods" of this object class.
    // jitters
    jitter() {// start jitter method
      this.x += tRNAlinearMaxAmplitude*random(-1, 1); // needs to be constrained in the canvas.
      this.y += tRNAlinearMaxAmplitude*random(-1, 1);
      this.x = constrain(this.x, 350, width);
      this.y = constrain(this.y, 0, height);
      this.angle += tRNAangularSpeed*random(-1.5, 1.5);
    } // end jitter method
    show(){// start show method
      strokeWeight(5);
      stroke(0, 25);
      fill(100, 250, 100, 35);
      push();
      translate(this.x, this.y);
      rotate(this.angle);
      // rect(0, 0, picWidth*5, picHeight);
      image(tRNA1pic, this.x, this.y, picWidth, picHeight);
      pop();

    }// end show method

} // end Trna class
