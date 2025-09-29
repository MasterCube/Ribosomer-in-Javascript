// Objects classes for Ribosomer (ribosome digital twins platform for translation simulation)
// class Trna:
var codonProcessed;
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
      this.x = constrain(this.x, 100, width-100);
      this.y = constrain(this.y, 300, height-400);
      this.angle += tRNAangularSpeed*random(-2, 2);
    } // end jitter method
    show(){// start show method
      strokeWeight(5);
      stroke(0, 25);
      fill(100, 250, 100, 35);
      push();
      translate(this.x, this.y);
      //rotate(this.angle);
      image(tRNA1pic, this.x, this.y, picWidth, picHeight);
      pop();

    }// end show method
  } // end Trna class

  //*************************************************************************
  // transcript class, properties and methods:
  //------------------------------------------

  // mRNA transcripts:
  // The length of the transcript is proportional to number of codons
  // 1 codon = 3 pixels
  class Lattice{
    constructor(tempStartX, tempStartY, tempLength, tempSequence, tempCopy, tempName, tempNumID, tempProList, tempProStatus, tempPlusList, tempMinusList, tempTunnelStatus){
      //properties:
        //property: position for display
      this.x=tempStartX;
      this.y=tempStartY;
      this.length=tempLength; // length is the number of codons in the transcript (stop codon included).
      this.sequence=tempSequence;
      this.mRNAlength= tempSequence.length; // mRNAlength is the number of nucleotides in the transcripts stop codon included
      this.copyIndex=tempCopy;
      this.name=tempName;  // should make sure it is an individual unique name (unique transcript ID)
      this.numID=tempNumID; // unique numeration of the different transcripts ID. (j ranges from 0 to J-1).
      this.uniqueID = str(tempName)+'#'+str(tempNumID);
      this.posProlines=tempProList; // this will content the list of the proline coding codons' positions for this transcript.
      this.proStatus = tempProStatus; // boolean is true if proline slowdown was activated in settings page
      this.plusList= tempPlusList; // list of the positions of + charged aa coding codons for this transcript.
      this.minusList=tempMinusList; // list of the positions of - charged aa coding codons for this transcript.
      this.tunnelStatus = tempTunnelStatus; //boolean is true if tunnel elec interaction was activated in settings page
      this.currentFootprint = [-1]; // array of ribosome positions when at P sites (nucleotide scale).
      this.translatingTime = []; // list of all translating times that each tanslating ribosomes have spent on this transcript
      this.readTranslated = 0; // only this copy of the transcript is iterated upon translation.
      this.countTranslated = 0; // CAUTION - aggregated over all the copies of this transcript name.
      this.rpfCount = 0; // number of RPF on the transcript.
  } // end constructor
      //methods:
      displays(){
        strokeWeight(4);
        strokeCap(ROUND);
        stroke(248, 196, 211); // rgb code for fast codon
        line(this.x, this.y, this.x + this.length*3, this.y); // 300 codons

        // highlight in purple the proline coding codons:
        if (this.posProlines.length > 0 && this.proStatus == "true"){
          var flagHere = 1;
          for (i=0; i<this.posProlines.length; i++){
            strokeWeight(2);
            stroke(153, 51, 102); // rgb for slow proline codon 129, 110, 208
            push();
            line(this.x + (this.posProlines[i])*3, this.y, this.x + (this.posProlines[i] + 1)*3, this.y);
            pop();
          }
        }// end if

        // paint in red the aa+ coding codons:
        if (this.plusList.length > 0 && this.tunnelStatus == "true"){
          var flagHere = 2;
          for (i=0; i<this.plusList.length; i++){
            strokeWeight(2);
            stroke(255, 0, 0); // rgb for positively charged aa coding codon 255, 0, 0
            push();
            line(this.x + (this.plusList[i])*3, this.y, this.x + (this.plusList[i] + 1)*3, this.y);
            pop();
          }
        }// end if

        // paint in blue the aa- coding codons:
        if (this.minusList.length > 0 && this.tunnelStatus == "true"){
          var flagHere = 3;
          for (i=0; i<this.minusList.length; i++){
            strokeWeight(2);
            stroke(0, 0, 255); // rgb for positively charged aa coding codon 0, 0, 255
            push();
            line(this.x + (this.minusList[i])*3, this.y, this.x + (this.minusList[i] + 1)*3, this.y);
            pop();
          }
        }// end if

        // add name of the transcript to the left as text:
        fill(0);
        stroke(235);
        textSize(18);//24
        textAlign(LEFT);
        push();
        text(this.name.substring(1, this.name.length-1), this.x - 115, this.y+8);
        pop();
        // add number of proteins (# termination events that were encountered for this transcript) synthesized:
        // since the beginning of initiation (beginning of simulation):
        fill(0);
        stroke(235);
        textSize(14);//24
        textAlign(LEFT);
        push();
        text(this.readTranslated, this.x - 115, this.y+24);
        // + display the current value of the total proteins counts produced
        // for this transcript and contributed by the other read copies:
        text(this.countTranslated, this.x - 85, this.y+24);
        // display the # of RPF (ribosome protected fragments) for this read
        text("RPF="+str(this.rpfCount), this.x -55, this.y+24);
        pop();

        // add last total elongation time for this transcript:


      } // end displays // end methods
  } // end class lattice

  // Class codon **************************************************************
  //---------------------------------------------------------------------------
  // codon class, properties and methods:
  class Codon{
    constructor(tempTranscript, tempPosition){
      this.type = tempTranscript.sequence.substring(tempPosition, tempPosition+3); //P site codon
      this.belongs = tempTranscript;
      this.position = tempPosition;// index position of codon in the trancript it comes from
      if (this.position <= tempTranscript.mRNAlength - 6){
        this.typeAsite = tempTranscript.sequence.substring(tempPosition+3, tempPosition+6); //A site codon
      }// end if
      if (this.position <= 147){
        this.upstreamWindow = tempTranscript.sequence.substring(0, tempPosition);
      }// end if
      if (this.position >= 148){
        this.upstreamWindow = tempTranscript.sequence.substring(tempPosition - 49*3, tempPosition+3);
      }// end if

    }// end constructor
    // methods:
    gammaDistParams(){
      if (this.position == this.belongs.mRNAlength - 3){
        this.alphaP = 1.0; // termination: the queueing time for termination is exponentially distributed
        this.beta = 100; // and is not a translation limiting factor.
      }// end if
      this.alphaP = 1.0; // termination: the queueing time for termination is exponentially distributed
      this.beta = 100; // and is not a translation limiting factor.
    }// end gammaDistParams method
    setQueueTime(){
      this.alphaP = 1.0; // termination: the queueing time for termination is exponentially distributed
      this.beta = 100; // and is not a translation limiting factor.
      return jStat.gamma.sample(this.alphaP, this.beta)
    } // end setQueueTime
    queueTimeAtPsiteSP(sp){ // sp = species ("yeast" or "coli" or...)
      if (this.position == this.belongs.mRNAlength - 3){
        this.alphaP = 1.0; // termination: the queueing time for termination is exponentially distributed
        this.beta = 100; // and is not a translation limiting factor.
        return jStat.gamma.sample(this.alphaP, this.beta)
        }// end if
      else{
        //console.log('still translocating...with queueuing time set at', HypoExpRandomTime(this.type, this.typeAsite, sp));
        //console.log('codon at P site is currently=', this.type, 'A type=', this.typeAsite);
        //console.log('belongs to transcript:', this.belongs, ' at position ', this.position)
        return HypoExpRandomTime(this.type, this.typeAsite, sp, this.upstreamWindow)
      }// end else
    }// end queueTimeAtPsiteSP method

  }// end class codon

  //*************************************************************************
  // ribosome class, properties and methods:
  //------------------------------------------
  class Ribo{
    constructor(RiboX, RiboY, sp){
      // has a position
      this.x = RiboX;
      this.y = RiboY;
      this.species = sp; // (other possibilities: 'coli', 'human', 'mouse', ...)
      // has a type ('free' or 'initiated' or 'translating')
      this.type = 'free';  // alternative type space: ('free', 'initiated', 'translating')
      //has time events records:
      this.clockTimer = millis();
      this.startedInitiation = 0;
      this.startedTranslocation = 0; // time record of the start of the translocation on a given site (codon)
      this.endedTermination = 0;
      // draw a timer setpoint from an exponential distribution with rate parameter called from initiation rate InitALPHA and InitBETA:
      this.InitiationTimerSP = jStat.gamma.sample(InitALPHA, InitBETA);
      this.TranslocationTimerSP = 0;
      this.CumulatedTranslocationSP = 0; // upon traffic jam, the translocation time setpoints are cumulated here.
      this.translocqueuingTimer = 0;
      this.timeSinceInitiation = 0;
      this.cumulatedTimeOnThisCodon = [];
      this.dwellTimePrevious = 0;
      this.transcriptName = '';
      this.uniqueReadID = '';
      this.goStatus = true;
      this.indexCodon = -3; // indexCodon is the position of the codon in the mRNA sequence (multiple of 3).
    } // first constructor to be overloaded.
      // overloading the constructor: by default only 2 arguments are provided upon instantiation of a new Ribo object: x and y positions.
      // But we could add a third argument specifying which transcript object the ribosome is on. In this case,
      // the extra third argument passed would be the transcript object upon initiation method.
      // end RiboFree class constructor
  print(){
    // this ribosome method prints at the console all the relevent attributes of the selected ribosome
    console.log('type=', this.type);
    console.log('clockTimer=', this.clockTimer);
    console.log('transcriptName=', this.transcriptName);
    console.log('uniqueReadID',this.uniqueReadID);
    console.log('indexCodon=', this.indexCodon);
    //console.log('codon=', Codon(lattice[0], this.indexCodon).type);
  }

  show() {// start show the free ribosome method:
    // has a geometry (attribute)
    // ---- RIBOSOME OBJECT DESIGN: ------------
    // show the ribosome only if they are free or if they are on the first 16 transcripts:
    if (this.type == 'free' || this.y <= 1080){
      strokeWeight(1); //(1)
      stroke(24);
      line(this.x-10, this.y, this.x+20, this.y); // x-15, x+15
      strokeWeight(4);
      stroke(24, 120, 54, 100);
      noStroke();
      // ---------- Large 60S subunit: -----------
      fill(24, 120, 54, 180);
      ellipse(this.x+3, this.y-22, 50, 40); // x+3 remove+3
      // ---------- Small 40S subunit: -----------
      rect(this.x-4, this.y+3, 14, 16); // x-4
      arc(this.x-4, this.y+11, 16, 16, HALF_PI, PI+HALF_PI); //x-4
      arc(this.x+10, this.y+11, 16, 16, PI+HALF_PI, HALF_PI); //x+7
      // -----------protein exit tunnel:----------
      stroke(205, 25, 25, 100);
      line(this.x-20, this.y-31, this.x, this.y-15);
      // ---- show specific information of a chosen transcript unique ID about codon at P sites
      // amino acid at P site, translocation timer set point, real queueing time on previous codon...
      // do it or not (if condition to be implemeneted)...
      var read_select_1 = str(StoredGeneID[1])+"#0";
      var read_select_2 = str(StoredGeneID[2])+"#0";
      // display ribosome footprinted nucleotides (RPF) length = 31 nts. Only the geneID 4 and 5:
      var read_select_3 = str(StoredGeneID[4])+"#0";
      var read_select_4 = str(StoredGeneID[5])+"#0";

      if (this.uniqueReadID==read_select_1 || this.uniqueReadID==read_select_2){//
        fill(0, 0, 255, 125);
        stroke(235);
        textSize(14);
        textAlign(LEFT);
        push();
        // codon at P site and A site:
        var text_codon = str(this.indexCodon);
        for (i_read=0;i_read<lattice.length; i_read++){
          if (lattice[i_read].uniqueID == this.uniqueReadID){
            var codonHereP = new Codon(lattice[i_read], this.indexCodon);
            text_codon = str(codonHereP.type);
            text(text_codon, this.x+25, this.y + 16);
            text(str(codonHereP.typeAsite), this.x+62, this.y + 16);
            // aa decoded:
            var text_aa = decode(codonHereP.type);
            text(text_aa, this.x+36, this.y - 11);
            if (text_aa != 'stop'){
              text_aa = decode(codonHereP.typeAsite);
              text(text_aa, this.x+73, this.y - 11);
              // when codon at P site: what is the set point queueing time
              //if (this.type != 'free'){
              //var text_sp = "sp = "+str(nfs(this.TranslocationTimerSP, 3, 1));
              //text(text_sp, this.x+25, this.y + 32);
              // dwell time previous:
              var text_qt = "last Psite q. time (ms) = "+str(nfs(this.dwellTimePrevious, 3, 1));
              text(text_qt, this.x-120, this.y + 32);
              // show the 31 nucleotides footprint of some RPF read: (some of them! not all!)
              // show codon upstream window:
              //var upseqWindow = codonHereP.upstreamWindow;
              //text(upseqWindow, this.x-120, this.y + 49);
              // display charged list:

              // display Maxwell Boltzmann factor:
              //var textMBf = "Maxwell-Boltzmann factor = "+str(nfs(MaxwellBoltzmannTunnelFactor(upseqWindow), 1, 3));
              //text(textMBf, this.x+109, this.y - 11);
              //console.log('this transcript upstream sign list:');
              //console.log('Maxwell Boltzmann =', MaxwellBoltzmannTunnelFactor(upseqWindow));
            }
          }
        }
        pop();
      }// end if display info on ribosomes of two chosen reads #1 and #2
      if ((this.uniqueReadID==read_select_3 || this.uniqueReadID==read_select_4) && displayRPFStatus=='true'){//
        console.log('+++++   This where we are gonna display the 31 nts    ++++++++++++++');
        fill(0, 0, 255, 125);
        stroke(235);
        textSize(12);//14
        textAlign(CENTER);//LEFT
        push();
        // codon at P site:
        var text_codon = str(this.indexCodon);
        for (i_read=0;i_read<lattice.length; i_read++){
          if (lattice[i_read].uniqueID == this.uniqueReadID){
            var codonHereP = new Codon(lattice[i_read], this.indexCodon);
            text_codon = str(codonHereP.type);
            text(text_codon, this.x, this.y + 31);
            var nt_pos = codonHereP.position;
            //var nt_pos = codon_pos * 3;
            // 31 nts footprint:
            if (nt_pos >= 0 && nt_pos <= lattice[i_read].mRNAlength){
              var nt_pos_inf = max(0, nt_pos - 12);
              var nt_pos_sup = min(nt_pos + 20, lattice[i_read].mRNAlength);
              //var nts_footprint ="5'-" + str(lattice[i_read].sequence.substring(nt_pos - 12, nt_pos + 20)) +"-3'";
              var nts_footprint ="5'-" + str(lattice[i_read].sequence.substring(nt_pos_inf, nt_pos_sup)) +"-3'";
              text(nts_footprint, this.x + 19, this.y + 45); //this.x - 4
              fill(0, 0, 0, 125);
              stroke(235);
              textSize(14);//14
              text('RPF: ', this.x - 145, this.y + 45);
            }
          }
        }//end for
        pop();
      }// end if display RPF 31 nts on these two reads
      console.log('currently the display status of RPF is :', displayRPFStatus);
    }
  }// end show the ribosome

  checkOccupancy(forThisTranscript){
    // this method checks if the next 6 codons from the current P site position of this ribosome are free to allow
    // the next translocation. It returns a boolean (called freeToGo) true if the ribosome is free to translocate or false if not.
    var freeToGo = true;
    if (this.uniqueReadID == forThisTranscript.uniqueID){ // are you on the right trancript ?
      // compare current index position of the ribosome with all occupation sites for this current transcript:
      for (let occupiedByOtherRibo = 0; occupiedByOtherRibo < forThisTranscript.currentFootprint.length; occupiedByOtherRibo ++){
        if ((this.indexCodon < forThisTranscript.currentFootprint[occupiedByOtherRibo]) &&
           (forThisTranscript.currentFootprint[occupiedByOtherRibo] - this.indexCodon <= riboFootprintLength)){
          freeToGo = false;
          break; // There is another ribosome just downstream this one in less than 6 codons (6 = 10 - offset) or (18 nt = 30 - 12 nt).
        }// end if
      } // end for
    }// end if
    return freeToGo;
  }

  translocates(tempTranscript){
    // this ribosome receives in argument a transcript object.
    // Monitor the timer and wait untill the setpoint for translocation is reached. When elapsed time larger than setpoint:
    // Checks if the ribosome is free to go by using the .checkOccupancy method.
    // if yes:
    // 1) this initiated ribosome or translating ribosome can translocates
    // 2) if the ribosome is currently on the last 'STOP' codon, it will terminate in an elapsed time span determined by the termination rate
    // and turn into a 'free' ribosome type.
    // Inquire about the codon type that is currently in tha A site. The A site is one codon downstream the P-site which is known from the current
    // position of the ribosome.

    // We kept track of the name of the current transcript on which this ribosome is on and the translocates method only applies to this transcript:
    // track of the name of transcript being translated.
    // check occupancy status downstream this codon:
    var freeToGo = true;
    this.goStatus=freeToGo;
    if (this.uniqueReadID == tempTranscript.uniqueID){
      // are you on the right trancript ?
      // compare current index position of the ribosome with all occupation sites for this current transcript:
      for (let occupiedByOtherRibo = 0; occupiedByOtherRibo < tempTranscript.currentFootprint.length; occupiedByOtherRibo ++){
        if ((this.indexCodon < tempTranscript.currentFootprint[occupiedByOtherRibo]) &&
           (tempTranscript.currentFootprint[occupiedByOtherRibo] - this.indexCodon <= riboFootprintLength)){
          freeToGo = false;
          this.goStatus = freeToGo;
          break; // There is another ribosome just downstream this one in less than 6 codons (6 = 10 - offset).
        }// end if
      } // end for
      //return freeToGo;
      var PsiteCodon = []; // this is defined as a list but will actually be used as a LIFO stack (push and pop)
      var AsiteCodon = []; // this is defined as a list but will actually be used as a LIFO stack (push and pop)
      //var CumulatedTimeOnThisCodon = []; //$$$$$$$$$$$$$$$$$$$**********************
      // this was defined as a list but will be used as a stack with push and pop...
      // ???? when do you set this.startedTranslocation = millis();
      //var queueingTime = (currentTime - this.startedTranslocation);
      var queueingTime = (millis() - this.startedTranslocation);

      //this.startedTranslocation = millis();
      //this.translocqueuingTimer =  millis();
      if (this.indexCodon <= tempTranscript.mRNAlength - 6){
        // identify and instantiate the codon object at the P site (including A-site context...)
        PsiteCodon.push(new Codon(tempTranscript, this.indexCodon));
        // identify and instantiate the codon object at the A site:
        //AsiteCodon.push(new Codon(tempTranscript, this.indexCodon + 3));
        codonProcessed = new Codon(tempTranscript, this.indexCodon);
        //console.log('transcript name=', tempTranscript.name);
        //console.log('transcript uniqueID', tempTranscript.uniqueID);
        //console.log('tempTranscript length=', tempTranscript.mRNAlength);
        //console.log('codon at current P site:', codonProcessed.type);
        //console.log('codon at current A site:', codonProcessed.typeAsite);
        // draw a timer setpoint for the queueing timeout:
        // queues to translocate:
        //this.TranslocationTimerSP = PsiteCodon.pop().queueTimeAtPsiteSP(this.species);
        this.TranslocationTimerSP = codonProcessed.queueTimeAtPsiteSP(this.species);
        //console.log('queueing sp=', this.TranslocationTimerSP);
        this.CumulatedTranslocationSP += this.TranslocationTimerSP;

        //queueingTime = (currentTime - this.startedTranslocation);
        if (((queueingTime >= (this.TranslocationTimerSP-cycle_time_draw_loop)) && (this.type != 'free')) && (freeToGo)){
          this.type = 'translating';
          // translocate:
          // update transcript data wrt ribosomal occupancy
          // the transcript has now a ribosomal footprint that must be updated:
          for( var i = 0; i < tempTranscript.currentFootprint.length; i++){
            if ( tempTranscript.currentFootprint[i] == this.indexCodon) {
              tempTranscript.currentFootprint[i] = this.indexCodon + 3;
              break;
              //tempTranscript.currentFootprint.splice(i, 1);
            }// end if
          }// end for
          // you have just translocated to the next codon, so you can record
          // the dwell time on the codon you are leaving and
          //save the cumulated SP (in case of traffic jam occured):
          var index_tr = indexOfRead(StoredGeneID, readCount, tempTranscript.uniqueID);
          RiboSeqSPtruth[index_tr][int(this.indexCodon/3)] += this.CumulatedTranslocationSP;
          if (this.indexCodon>=3){
            RiboSeqRRTtruth[index_tr][int(this.indexCodon/3)-1] += this.dwellTimePrevious;
          }
          // update ribosome position:
          this.indexCodon += 3; // upon translocation, the ribosome P site is on the next codon;
          this.x += 3;// x position of the next codon of the current transcript: 0 = A of AUG at P site.
          this.y = tempTranscript.y;// y position of the dito
          // reset the cumulated time on this codon (in case traffic jam occured):
          this.CumulatedTranslocationSP = 0;
          //CumulatedTimeOnThisCodon += queueingTime; // bof ???
          //this.dwellTimePrevious = CumulatedTimeOnThisCodon;
          this.cumulatedTimeOnThisCodon.push(millis());
          //this.dwellTimePrevious = this.cumulatedTimeOnThisCodon.slice(-2)[1] - cumulatedTimeOnThisCodon.slice(-2)[0];
          this.dwellTimePrevious = this.cumulatedTimeOnThisCodon[this.cumulatedTimeOnThisCodon.length - 1] - this.cumulatedTimeOnThisCodon[this.cumulatedTimeOnThisCodon.length - 2];
          //this.dwellTimePrevious = queueingTime;
          // reset timer
          this.startedTranslocation = millis(); // reinitialize the timer
        }
        if (!freeToGo){// there is another ribosome downstream that jams the translocation:
          codonProcessed = new Codon(tempTranscript, this.indexCodon);
          this.TranslocationTimerSP = codonProcessed.queueTimeAtPsiteSP(this.species);
          this.CumulatedTranslocationSP += this.TranslocationTimerSP;
          //this.TranslocationTimerSP=0;
          // reset translocation timer
          this.startedTranslocation =  millis();
        }// end if

      }// end if
      if (this.indexCodon == tempTranscript.mRNAlength - 3) {// the ribosome position is at the stop codon. The next translocation is here a termination step.
        // queues to terminate:
        var finalCodon = new Codon(tempTranscript, this.indexCodon);
        //this.TranslocationTimerSP = finalCodon.setQueueTime();
        this.TranslocationTimerSP = finalCodon.queueTimeAtPsiteSP(this.species);
        this.CumulatedTranslocationSP = this.TranslocationTimerSP;
        if ((queueingTime >= (this.TranslocationTimerSP-cycle_time_draw_loop)) && (this.type != 'free')){
          //CumulatedTimeOnThisCodon += queueingTime; // bof ???
          //this.dwellTimePrevious = CumulatedTimeOnThisCodon;
          // reset the started translocation timer:
          //this.startedTranslocation = 0;
          this.cumulatedTimeOnThisCodon.push(millis());
          //this.dwellTimePrevious = this.cumulatedTimeOnThisCodon.slice(-2)[1] - cumulatedTimeOnThisCodon.slice(-2)[0];
          this.dwellTimePrevious = this.cumulatedTimeOnThisCodon[this.cumulatedTimeOnThisCodon.length - 1] - this.cumulatedTimeOnThisCodon[this.cumulatedTimeOnThisCodon.length - 2];
          this.startedTranslocation = millis();
          // the current transcript is now left by the ribosome and the ribosome attribute for the current transcript will have to become empty:
          //this.type = 'free';
          this.transcriptName = '';
          this.uniqueReadID = '';

          // update transcript data wrt ribosomal occupancy
          // the transcript has now a ribosomal footprint that must be updated:
          for( var i = 0; i < tempTranscript.currentFootprint.length; i++){
            if ( tempTranscript.currentFootprint[i] == this.indexCodon) {
              tempTranscript.currentFootprint.splice(i, 1); // remove 1 element at index i. Correct.
              break;
            }// end if
          }// end for

          // relocate to a free ribosome position and re-order all the free ribosomes at the right place:
          // count the number of free ribosome before this last termination:
          var NfreeRib = 0;
          var NinitiatedRib = 0;
          for (let rib = 0; rib < freeRibosomeInitialNumber; rib +=1){
            if (ribosome[rib].type == 'free'){
              NfreeRib += 1;
            }
            if (ribosome[rib].type != 'free'){
                NinitiatedRib += 1;
            }
          }// end for
          // (free) previously existing ribosomes re-location:
          const riboInterY = 80;
          const riboInterX = 60;
          const riboPerLine = 20;
          var numbRiboLines = int(freeRibosomeInitialNumber/riboPerLine) + 1;
          var freeI = 0;
          for (let nbR=0; nbR < ribosome.length; nbR++){
            if (ribosome[nbR].type == 'free'){
              //ribosome[nbR].x = 100 + freeI*75;
              ribosome[nbR].x = 80+(freeI % riboPerLine)*riboInterX;
              //ribosome[nbR].y = 380;
              i_line = int(freeI/riboPerLine);
              ribosome[nbR].y = 1200 + i_line * riboInterY;
              ribosome[nbR].clockTimer = millis();
              ribosome[nbR].InitiationTimerSP = jStat.gamma.sample(InitALPHA, InitBETA);
              ribosome[nbR].TranslocationTimerSP = 0;
              ribosome[nbR].CumulatedTranslocationSP = 0;
              ribosome[nbR].translocqueuingTimer = 0;
              ribosome[nbR].timeSinceInitiation = 0;
              ribosome[nbR].transcriptName = '';
              ribosome[nbR].uniqueReadID = '';
              ribosome[nbR].indexCodon = -3;
              ribosome[nbR].startedInitiation = millis();
              ribosome[nbR].startedTranslocation = millis();
              ribosome[nbR].dwellTimePrevious = 0;
              ribosome[nbR].cumulatedTimeOnThisCodon = [];
            }// end if
            freeI += 1;
          }// end for
          this.type = 'free';
          this.indexCodon = -3;
          this.transcriptName = ''; // The ribosome is leaving the transcript and forgets the transcript name.
          this.uniqueReadID = '';
          this.x = 80 + ((NfreeRib+1) % riboPerLine)*riboInterX;
          //this.y = 380;
          i_line = int((NfreeRib+1)/riboPerLine);
          this.y = 1200 + i_line * riboInterY;
          this.clockTimer = millis();
          this.InitiationTimerSP = jStat.gamma.sample(InitALPHA, InitBETA);
          this.TranslocationTimerSP = 0;
          this.CumulatedTranslocationSP = 0;
          this.translocqueuingTimer = 0;
          this.timeSinceInitiation = 0;
          this.startedInitiation = millis();
          this.startedTranslocation = millis();
          this.dwellTimePrevious = 0;
          this.cumulatedTimeOnThisCodon = [];
        }// end if (translocation condition)
        // a protein molecule was just produced: update the protein object:
        // or update the count of the complete translation events for this transcript unique read:
        tempTranscript.readTranslated += 1;
      }// end if (for case on stop codon)
    }// end if (for the current transcript name)

  }// end translocates method
  initiates(tempTranscript){ // a free ribosome may initiate on a given available transcript (unique read)
    // a free ribosome receives in argument a randomly chosen transcript (unique read).
    // monitor the timer and wait until time setpoint is reached. When elapsed time larger than setpoint:
    // Checks if the first 10 codons (from 0 to 9 included) of the given transcript are free of any previous P-site occupying ribosomes
    // if yes:
    // 1) this free ribosome becomes an initiated ribosome and this free ribosome disappears (it has switched from population)
    // 2) the new position is the one of the start codon (the P site is at 0 of the transcript)
    // if no:
    // 1) the ribosome stays free
    // 2) its timer is reset to zero.

    var gotFootprint = tempTranscript.currentFootprint;
    //console.log('gotFootprint=', gotFootprint);
    var footprinted = false;
    var listLength = gotFootprint.length;
    for (let psite=0; psite <= listLength; psite +=1){
      if ((gotFootprint[psite] >= 0) && (gotFootprint[psite] < riboFootprintLength)){
        footprinted = true;
        break; // you can get out of the for loop if a footprint was detected
      } // end if
    } // end for loop
    if (((currentTime - this.clockTimer) >= this.InitiationTimerSP) && (this.type=='free') && (footprinted == false)){
      //turn the ribosome into an initiated ribosome:
      this.type = 'initiated';
      this.startedInitiation = millis(); // record time event when initiation started
      this.startedTranslocation = millis(); // record time event when queueing time for next translocation started
      // keep track of the name of the current transcript on which this ribosome is on:
      this.transcriptName = tempTranscript.name; // retrieve name of transcript being intiated or translated
      this.uniqueReadID = tempTranscript.uniqueID;
      this.indexCodon = 0; // upon intiation, the ribosome P site is on codon 0 = 'AUG'
      // update transcript data wrt ribosomal occupancy
      this.x = tempTranscript.x;// x position of the start codon 'AUG' of the current transcript: 0 = A of AUG at P site.
      this.y = tempTranscript.y;// y position of the dito
      // the transcript has now a ribosomal footprint that must be updated:
      if (tempTranscript.currentFootprint[0] == -1){
        tempTranscript.currentFootprint.pop(); // remove the virtual [-1] footprint of the former free ribosome
      }
      tempTranscript.currentFootprint.push(0); // an initiating ribosome (P-site) occupies the start codon, it has footprint 0 to 5 codons included (= 6)
      this.cumulatedTimeOnThisCodon.push(millis());
    }// end if
    this.InitiationTimerSP = jStat.gamma.sample(InitALPHA, InitBETA); // redraw a new random setpoint from the exponential distribution
    this.clockTimer = millis();
  }// end initiate method

  } // end free ribosome class

// function to test local storage remaining memory:
function testLocalStorage() {
  var timeStart = Date.now();
  var timeEnd, countKey, countValue, amountLeft, itemLength;
  var occupied = leftCount = 3; //Shurav's comment on initial overhead
  //create localStorage entries until localStorage is totally filled and browser issues a warning.
  var i = 0;
  while (!error) {
    try {
        //length of the 'value' was picked to be a compromise between speed and accuracy,
        // the longer the 'value' the quicker script and result less accurate. This one is around 2Kb
        localStorage.setItem('testKey' + i, '11111111112222222222333333333344444444445555555555666661111111111222222222233333333334444444444555555555566666');
      } catch (e) {
            var error = e;
      }
      i++;
  }
  //if the warning was issued - localStorage is full.
  if (error) {
    //iterate through all keys and values to count their length
    for (var i = 0; i < localStorage.length; i++) {
      countKey = localStorage.key(i);
      countValue = localStorage.getItem(localStorage.key(i));
      itemLength = countKey.length + countValue.length;
      //if the key is one of our 'test' keys count it separately
      if (countKey.indexOf("testKey") !== -1) {
        leftCount = leftCount + itemLength;
        }
      //count all keys and their values
      occupied = occupied + itemLength;
      }
      ;
      //all keys + values lenght recalculated to Mb
      occupied = (((occupied * 16) / (8 * 1024)) / 1024).toFixed(2);
      //if there are any other keys then our 'testKeys' it will show how much localStorage is left
      amountLeft = occupied - (((leftCount * 16) / (8 * 1024)) / 1024).toFixed(2);
      //iterate through all localStorage keys and remove 'testKeys'
      Object.keys(localStorage).forEach(function(key) {
        if (key.indexOf("testKey") !== -1) {
          localStorage.removeItem(key);
        }
      }); //});
    }
  //calculate execution time
  var timeEnd = Date.now();
  var time = timeEnd - timeStart;
  //create message
  var message = 'Finished in: ' + time + 'ms \n total localStorage: ' + occupied + 'Mb \n localStorage left: ' + amountLeft + "Mb";
  // put the message on console:
  console.log(message);
  //put the message on the screen
  //document.getElementById('scene').innerText = message; //this works with Chrome,Safari, Opera, IE
  //document.getElementById('scene').textContent = message;  //Required for Firefox to show messages
}
