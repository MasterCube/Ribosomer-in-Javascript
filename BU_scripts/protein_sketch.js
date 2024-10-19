// protein_sketch
// WARNING: we move to instance mode because we will create 2 or more canvas_settings
var sketch1 = function(p){
  p.x = 400;
  p.y = 180;
  p.i = 0;
  p.tr = 0;
  p.StoredTimes = [];
  p.i_length = 36;
  p.time_clock = 0;
  p.time_min = 0;
  p.time_sec = 0;
  p.time_text = "";
  p.time_interval_to_display = 1; // 2 seconds...
  //p.timeSeqFlag = [];
  p.condition = false;
  p.lineOne = "";
  p.ord_label_one = "";
  p.abs_label_one = "";
  p.ticklength = 7;
  p.tick_delta_x = 55;
  p.px_tick = [];
  p.iMax = 0; // number of transcripts that were read
  p.iMin = 0; // minimal number of transcripts or numb of bars
  p.numbOfbars = 0;
  p.StoredProtCount = [];
  p.string_array;
  p.arrayOfStrings = [];
  p.protCountMaxNew = 0;
  p.protCountMax = 0;
  p.StoredGeneID =[];
  p.keyOne = "";
  p.textAngle = 45;
  p.pxInf = 35;
  p.pyInf = 580;
  p.pxSup = 1465;
  p.pySup = 50;
  p.py_tick_max = 75;
  p.yMap = 100; // height of bar representing the protein count
  p.setup = function(){
    p.canvasForProteins =  p.createCanvas(1500, 600); // 1920, 938
    p.canvasForProteins.parent('#canvasProtein');
    //------------------------------------------------------------------------
    p.frameRate(50); // use 18 frames/s or 12/s // 6 frame/s = 167 ms between two frames.//12
    p.background(200);
    for (p.i=0; p.i<24; p.i++){
      p.tick_delta_x = 55;
      p.px_tick[p.i] = p.pxInf+115 + p.tick_delta_x * p.i;
    }
    // store 36 time sequence flags for 36 sequence of 3 seconds display interval:
    //for (p.i=0; p.i<36; p.i++){
    //  p.timeSeqFlag[p.i] = 'false';
    //}
  }
  //
  p.draw = function(){
    p.background(255, 255, 255);
    p.fill(0, 90, 90); // orange: 255, 127, 80
    p.stroke(235);
    p.textSize(16);//24
    p.textFont("Helvetica"); //Calibri, Helvetica, Georgia, Courrier new ?
    p.push();
    p.time_clock = p.millis();
    p.time_min = p.floor(p.time_clock/(1000*60));
    p.time_sec = p.int((p.time_clock/1000) % 60);
    p.time_text = "...";

    for (p.i=0; p.i<p.i_length; p.i++){
      p.StoredProtCount[p.i] = []; //remember this is a 2D array instantiation
      if (p.getItem('timePoint#'+p.str(p.i)) != null){
        p.StoredTimes[p.i] = p.getItem('timePoint#'+p.str(p.i));
        //console.log('storedTimes=', p.i, p.StoredTimes[p.i]);
        //p.storeItem('timePoint#'+p.str(p.i), p.StoredTimes[p.i]);
        // for each of the 36 time points:
        // retrieve the protein counts from local storage:
        p.StoredProtCount[p.i] = p.getItem("protCountTimePoint#"+p.str(p.i));
        // if the retrieved field is recognized as a string, convert it into an array of integer numbers!
        // [178,163,187,56,40,52,149,78,196,94,60,232,65,45,60,47,53,61,573,51,48,52,55,78,47,57,291]
        // convert this string into an array of integers...
        if (typeof p.StoredProtCount[p.i] == "string"){
          // remove both brackets:
          p.string_array = p.StoredProtCount[p.i].substring(1, p.StoredProtCount[p.i].length - 1);
          p.StoredProtCount[p.i] = [];
          p.arrayOfStrings = p.string_array.split(/,/);
          for (p.ind=0; p.ind < p.arrayOfStrings.length; p.ind++){
            p.StoredProtCount[p.i][p.ind] = p.int(p.arrayOfStrings[p.ind]);
          }
        }
        p.protCountMaxNew = p.max(p.StoredProtCount[p.i]);
        // for all the 36 time points find the maximum count of proteins:
        if (p.protCountMaxNew > p.protCountMax){
          p.protCountMax = p.protCountMaxNew; //p.protCountMax is the maximum level encountered in the simulation
        }
      }
    }
    p.pop();
    //******
    //  display the graphic box:
    // draw the border of the graphic box:
    p.push();
    p.strokeWeight(4);
    p.stroke(0, 90, 90);
    p.line(p.pxInf, p.pyInf, p.pxSup, p.pyInf);
    p.line(p.pxInf, p.pyInf, p.pxInf, p.pySup);

    p.line(p.pxInf, p.pySup, p.pxSup, p.pySup);
    p.line(p.pxSup, p.pyInf, p.pxSup, p.pySup);
    p.pop();
    // draw the abscissa and ordinate
    p.push();
    p.strokeWeight(2);
    p.stroke(150);
    p.line(p.pxInf+60, p.pyInf-50, p.pxSup-30, p.pyInf-50); // abscissa
    p.line(p.pxInf+60, p.pyInf-50, p.pxInf+60, p.pySup+35); //ordinate
    p.pop();

    // graphic box title and axes labels
    p.lineOne = "Cumulated proteins counts";
    p.textFont("Helvetica"); //Courier or Helvetica or Calibri
    p.textSize(24);
    p.fill(0, 90, 90); //(75, 0, 130) // coral color 255, 127, 80
    p.text(p.lineOne, p.pxInf+10, p.pySup-20);

    // axes labels:
    p.ord_label_one = "protein count";
    p.textFont("Helvetica"); //Courier or Helvetica
    p.textSize(20);
    p.fill(100, 105, 100); //(75, 0, 130)
    p.text(p.ord_label_one, p.pxInf+12, p.pySup+25);

    p.abs_label_one = "transcript ID";
    p.textFont("Helvetica"); //Courier or Helvetica
    p.textSize(20);
    p.fill(100, 105, 100); //(75, 0, 130)
    p.text(p.abs_label_one, p.pxInf + 700, p.pyInf - 25);

    // draw ticks :
    p.py_tick_max = p.map(p.protCountMax, 0, p.protCountMax*1.15, p.pyInf-50, p.pySup+35);
    //console.log('py_tick_max=', p.py_tick_max);
    p.push();
    p.strokeWeight(2);
    p.stroke(150); // stroke(100, 0, 170)
    p.ticklength = 7;
    for (p.i=0; p.i<24; p.i++){
      p.line(p.px_tick[p.i], p.pyInf - 50 - p.ticklength, p.px_tick[p.i], p.pyInf - 50);
    }
    p.line(p.pxInf + 60, p.py_tick_max, p.pxInf + 60 + p.ticklength, p.py_tick_max);
    p.pop();

    // retrieve the 27 stored first transcript ids
    p.iMax = p.getItem("numbLines");
    //p.storeItem("numbLines", p.iMax);
    for (p.i=0; p.i<p.iMax;p.i++){
      p.keyOne = p.getItem(p.str(p.i));
      //p.storeItem(p.str(p.i), p.keyOne);
      // remove the parentheses () around the geneID:
      p.StoredGeneID[p.i] = p.keyOne.substring(1, p.keyOne.length-1);
    }

    // retrieve the protein counts from local storage:
    //p.StoredProtCount[i] = getItem('protCountTimePoint#'+str(i));
    //storeItem('protCountTimePoint#'+str(i), p.StoredProtCount[i]);

    p.numbOfbars=24;
    p.iMin = p.min(p.iMax, p.numbOfbars);
    // Each 3 seconds interval, display sequentially the evolution of 36 time points of protein counts
    // display the protein counts
    //p.pyInf-50, p.pySup+35 (pyInf = 580, pySup = 50)
    for (p.i=0;p.i<p.i_length;p.i++){
      //if (p.time_clock/1000 >= (p.i+1)*p.time_interval_to_display && p.timeSeqFlag[p.i] == 'false'){
        //p.timeSeqFlag[p.i] = 'true';
      p.condition = ((p.i == p.i_length-1) | (p.time_clock/1000 >= (p.i+1)*p.time_interval_to_display && p.time_clock/1000 < (p.i+2)*p.time_interval_to_display));
      if (p.condition){
      //if (p.time_clock/1000 >= (p.i+1)*p.time_interval_to_display && p.time_clock/1000 < (p.i+2)*p.time_interval_to_display){
        // display the time point on the screen:
        p.push();
        p.textFont("Helvetica"); //Courier or Helvetica
        p.textSize(24);
        p.fill(0, 90, 90); //(75, 0, 130)
        p.time_text = 'at simulation time '+ p.str(p.StoredTimes[p.i]);
        if (p.i != p.i_length-1){
          p.text(p.time_text, 415, 30);
        }
        p.time_text = '/ total time '+ p.str(p.StoredTimes[p.i_length-1]);
        p.text(p.time_text, 950, 30);
        p.pop();
        for(p.tr=0; p.tr<p.iMin; p.tr++){// scan and display the protein levels of each geneID for the stored time point
          // scale the appropriate height for the bar representing the given protein level
          p.yMap = p.map(p.StoredProtCount[p.i][p.tr], 0, p.protCountMax*1.15, p.pyInf-50, p.pySup+35);
          // display the # protein number into text:
          p.push();
          p.textFont("Helvetica"); //Courier or Helvetica
          p.textSize(20);
          p.fill(100, 105, 100); //(75, 0, 130)
          p.text(p.StoredProtCount[p.i][p.tr], p.px_tick[p.tr]-35, p.yMap-2);
          p.pop();

          // display bars as lines with appropriate strokeWeight:
          //p.push();
          //p.strokeWeight(16);
          //p.stroke(0, 90, 90, 125); // transparency is 100 between 0-250 (250 is opaque).
          //p.line(p.px_tick[p.tr], p.pyInf-50, p.px_tick[p.tr], p.yMap);
          //p.pop();

          //alternatively display bars as rectangle:
          //...
          p.push();
          p.fill(0, 90, 90, 125);
          p.strokeWeight(0);
          p.stroke(0, 90, 90, 125); // transparency is 100 between 0-250 (250 is opaque).
          p.rect(p.px_tick[p.tr]-10, p.pyInf-50, 20, p.yMap-p.pyInf+50);
          p.pop();

          // add a tick above the bar or rectangle:
          p.push();
          p.strokeWeight(2);
          p.stroke(150); // stroke(100, 0, 170)
          p.ticklength = 7;
          p.line(p.px_tick[p.tr], p.yMap - 2 * p.ticklength, p.px_tick[p.tr], p.yMap);
          p.pop();

          // add a dashed line from the top tick to the geneID labels
          p.push();
          p.strokeWeight(2);
          p.stroke(150, 125); // stroke(100, 0, 170)
          p.drawingContext.setLineDash([8, 16]);
          p.line(p.px_tick[p.tr], p.yMap - 2 * p.ticklength, p.px_tick[p.tr], p.py_tick_max - 1*p.ticklength );
          p.pop();
        }
      }
    }

    // display their labels above the bars:
    for (p.i=0; p.i<p.iMin;p.i++){
      p.push();
      //p.translate(p.px_tick[p.i], p.pyInf - 50 - 2*p.ticklength);
      p.translate(p.px_tick[p.i], p.py_tick_max - 2*p.ticklength);
      p.rotate(p.radians(-50));
      p.text(p.StoredGeneID[p.i], 0, 0);
      //p.text(p.StoredGeneID[p.i], p.px_tick[p.i], p.pyInf - 50 - 2*p.ticklength);
      p.pop();
    }
  }//end of p.draw
}// end var sketch 1

var sketch2 = function(p){
  p.i = 0;
  p.tr = 0;
  p.textAngle = 45;
  p.pxInf = 35;
  p.pyInf = 580;
  p.pxSup = 1465;
  p.pySup = 50;
  p.pySHIFT = p.pyInf + 30;
  p.py_tick_max = 75;
  p.yMap = 100; // height of bar representing the protein count
  p.yMapRPF = 100; // height of bar representing the RPF count on this geneID divided by # transcripts reads copies.
  p.StoredTimes = [];
  p.i_length = 36;
  p.time_clock = 0;
  p.time_min = 0;
  p.time_sec = 0;
  p.time_text = "";
  p.time_interval = 3; // 3 minutes or 5 minutes (time span between two storage time event in ribosomer...)
  p.time_unit = 1;
  p.time_interval_to_display = 1; // 2 seconds...
  //p.timeSeqFlag = [];
  p.condition = false;
  p.lineTwo = "";
  p.lineThree ="";
  p.ord_label_one = "";
  p.abs_label_one = "";
  p.ord_label_two = "";
  p.abs_label_two = "";
  p.ticklength = 7;
  p.tick_delta_x = 55;
  p.px_tick = [];
  p.iMax = 0; // number of transcripts that were read
  p.iMin = 0; // minimal number of transcripts or numb of bars
  p.numbOfbars = 0;
  p.StoredProtCount = [];
  p.protCountMaxNew = 0;
  p.protCountMax = 0;
  p.StoredTransCount = [];
  p.StoredRPFCount = [];
  p.string_array = "";
  p.arrayOfStrings = [];
  p.StoredTEfromProteomic = []; // time points (will be a 2D array below).
  p.StoredTEfromProtavg = []; //averaged over time points for each transcript (is a 1D array).
  p.avg_PROT_over_timePoint = 0;
  p.TEfromProtMax = 0;
  p.TEmaxNew = 0;
  p.TEfromRPFMax = 0;
  p.TEfromRPFmaxNew = 0;
  p.StoredTEfromRPF = []; // time points (will be a 2D array below).
  p.StoredTEfromRPFavg = []; //averaged over time points for each transcript (is a 1D array).
  p.avg_RPF_over_timePoint = 0;
  p.protTransMaxNew = 0;
  p.protTransMax = 0;
  //p.rpfCountMaxNew = 0;
  //p.rpfCountMax = 0;
  p.StoredGeneID = [];
  p.mRNAlength = []; // list of lengths of transcripts in nts
  p.keyOne = "";
  // setup function
  p.setup = function(){
    p.canvasForProteins =  p.createCanvas(1500, 1200); // 1920, 938
    p.canvasForProteins.parent('#canvasProtein');
    //------------------------------------------------------------------------
    p.frameRate(60); // use 18 frames/s or 12/s // 6 frame/s = 167 ms between two frames.//12
    p.background(255, 255, 255);

    for (p.i=0; p.i<24; p.i++){
      p.tick_delta_x = 55;
      p.px_tick[p.i] = p.pxInf+115 + p.tick_delta_x * p.i;
    }

  }
  // draw function
  p.draw = function(){
    p.background(255, 255, 255);
    p.fill(0, 90, 90); // orange: 255, 127, 80
    p.stroke(235);
    p.textSize(16);//24
    p.textFont("Helvetica"); //Calibri, Helvetica, Georgia, Courrier new ?
    p.push();
    p.time_clock = p.millis();
    p.time_min = p.floor(p.time_clock/(1000*60));
    p.time_sec = p.int((p.time_clock/1000) % 60);
    p.time_text = "...";

    for (p.i=0; p.i<p.i_length; p.i++){
      p.StoredProtCount[p.i] = []; //remember this is a 2D array instantiation
      p.StoredRPFCount[p.i] = [];
      p.StoredTEfromProteomic[p.i] = []; //remember this is also a 2D array instantiation
      p.StoredTEfromRPF[p.i] = [];
      if (p.getItem('timePoint#'+p.str(p.i)) != null){
        p.StoredTimes[p.i] = p.getItem('timePoint#'+p.str(p.i));
        //p.storeItem('timePoint#'+p.str(p.i), p.StoredTimes[p.i]);
        // for each of the 36 time points:
        // retrieve the protein counts from local storage:
        p.StoredProtCount[p.i] = p.getItem('protCountTimePoint#'+p.str(p.i));
        // if the retrieved field is recognized as a string, convert it into an array of integer numbers!
        // [178,163,187,56,40,52,149,78,196,94,60,232,65,45,60,47,53,61,573,51,48,52,55,78,47,57,291]
        // convert this string into an array of integers...
        if (typeof p.StoredProtCount[p.i] == "string"){
          // remove both brackets:
          p.string_array = p.StoredProtCount[p.i].substring(1, p.StoredProtCount[p.i].length - 1);
          p.StoredProtCount[p.i] = [];
          p.arrayOfStrings = p.string_array.split(/,/);
          for (p.ind=0; p.ind < p.arrayOfStrings.length; p.ind++){
            p.StoredProtCount[p.i][p.ind] = p.int(p.arrayOfStrings[p.ind]);
          }
        }
        p.StoredRPFCount[p.i] = p.getItem('rpfCountTimePoint#'+p.str(p.i));
        //console.log('p.StoredRPFCount[p.i] as retrived from rpfCountTimePoint#', p.i, p.StoredRPFCount[p.i]);
        //console.log('type upon read is string ? ', typeof p.StoredRPFCount[p.i] == "string");
        if (typeof p.StoredRPFCount[p.i] == "string"){
          // remove the brackets:
          p.string_array = p.StoredRPFCount[p.i].substring(1, p.StoredRPFCount[p.i].length - 1);
          p.StoredRPFCount[p.i] = [];
          p.arrayOfStrings = p.string_array.split(/,/);
          for (p.ind = 0; p.ind < p.arrayOfStrings.length; p.ind++){
            p.StoredRPFCount[p.i][p.ind] = p.int(p.arrayOfStrings[p.ind]);
          }
        }

        p.protCountMaxNew = p.max(p.StoredProtCount[p.i]);
        // for all the 36 time points find the maximum count of proteins:
        if (p.protCountMaxNew > p.protCountMax){
          p.protCountMax = p.protCountMaxNew; //p.protCountMax is the maximum level encountered in the simulation
        }
      }
    }
    p.pop();
    //******
    //  display the graphic box:
    // draw the border of the graphic box:
    p.push();
    p.strokeWeight(4);
    p.stroke(0, 90, 90);
    p.line(p.pxInf, p.pyInf, p.pxSup, p.pyInf);
    p.line(p.pxInf, p.pyInf, p.pxInf, p.pySup);

    p.line(p.pxInf, p.pySup, p.pxSup, p.pySup);
    p.line(p.pxSup, p.pyInf, p.pxSup, p.pySup);
    p.pop();
    // draw the abscissa and ordinate
    p.push();
    p.strokeWeight(2);
    p.stroke(150);
    p.line(p.pxInf+60, p.pyInf-50, p.pxSup-30, p.pyInf-50); // abscissa
    p.line(p.pxInf+60, p.pyInf-50, p.pxInf+60, p.pySup+35); //ordinate
    p.stroke(255, 127, 80);
    p.line(p.pxInf+60+12, p.pyInf-50, p.pxInf+60+12, p.pySup+35); //second ordinate coral

    p.pop();

    // graphic box title and axes labels
    p.lineTwo = "Translation efficiency from proteins counts: TE = # proteins counts per time unit / # transcripts copy number [per hour].";
    p.lineTwoBis = "Translation efficiency expected from Ribo-Seq: TE = # Ribosome protected fragments (RPF) per 1 kb / # transcripts copy number [averaged].";
    p.textFont("Lato"); //Courier or Helvetica or Calibri
    p.textSize(22);//24
    p.fill(0, 90, 90); //(75, 0, 130) // coral color 255, 127, 80
    p.text(p.lineTwo, 40, 17);
    p.fill(255, 127, 80); //(75, 0, 130) // coral color 255, 127, 80
    p.text(p.lineTwoBis, 40, 43);

    // retrieve the 27 stored first transcript ids, mRNA lengths and reads counts:
    p.iMax = p.getItem("numbLines");
    //p.storeItem("numbLines", p.iMax);
    for (p.i=0; p.i<p.iMax;p.i++){
      p.keyOne = p.getItem(p.str(p.i));
      //p.storeItem(p.str(p.i), p.keyOne);
      // remove the parentheses () around the geneID:
      p.StoredGeneID[p.i] = p.keyOne.substring(1, p.keyOne.length-1);
      // retrieve the copy number for each transcript:
      if (p.getItem("nonUniformmRNAabundanceStatus")=='true'){
        p.StoredTransCount[p.i] = p.getItem('readCopiesOf_'+p.str(p.i));
      }
      else{
        p.StoredTransCount[p.i] = 1;
      }
      // retrieve the length of each transcript in nts:
      p.mRNAlength[p.i] = p.getItem(p.keyOne).length;
    }
    // TE as defined by proteomic:
    // and TE as determined from Ribosome protected fragments (RPF) count:
    //--------------------------------------------------------------------
    p.numbOfbars=24;
    p.iMin = p.min(p.iMax, p.numbOfbars);
    for (p.i=0; p.i<p.i_length; p.i++){
      for(p.tr=0; p.tr<p.iMin; p.tr++){// scan, compute and display the TE from prot divided by read counts of each geneID
        // the following quotient is standardized by unit of time interval (3 min rescaled on 1 hour):
        p.time_unit = p.time_interval * (p.i + 1) / 60; // multiple of 3 minutes (or 5 minutes CAUTION !!!!!) rescaled to 1 hour (60 min).
        p.StoredTEfromProteomic[p.i][p.tr] = p.StoredProtCount[p.i][p.tr] / (p.StoredTransCount[p.tr] * p.time_unit);
        // RPF: normalized by mRNA length (rescaled for 1kb) and by read copy number:
        p.StoredTEfromRPF[p.i][p.tr] = p.StoredRPFCount[p.i][p.tr] *1000/ (p.StoredTransCount[p.tr] * p.mRNAlength[p.tr]); // this is instantaneous at time point
      }
      //console.log('TE from proteomic at time point=', p.StoredTEfromProteomic[p.i]);
    }
    // find max value of TE
    p.TEfromProtMax = 0;
    p.TEfromRPFMax = 0;
    for (p.i=0; p.i<p.i_length; p.i++){
      // for all the 36 time points find the maximum count of proteins:
      p.TEmaxNew = p.max(p.StoredTEfromProteomic[p.i]);
      if (p.TEmaxNew > p.TEfromProtMax){
        p.TEfromProtMax = p.TEmaxNew; //p.TEfromProtMax is the maximum TE level encountered in the simulation
      }
      // ... and the maximum count of RPF:
      p.TEfromRPFmaxNew = p.max(p.StoredTEfromRPF[p.i]);
      if (p.TEfromRPFmaxNew > p.TEfromRPFMax){
        p.TEfromRPFMax = p.TEfromRPFmaxNew;
      }
    }
    // when and only when you have all the time points, compute the bulk average of TE from RPFs for each transcript:
    // idem for TE from proteomic:
    if (p.StoredTEfromRPF[p.StoredTimes.length - 1][0] > 0 || p.StoredTEfromRPF[p.StoredTimes.length - 1][p.iMin - 1] > 0){
      for (p.tr=0; p.tr<p.iMin; p.tr++){
        p.avg_RPF_over_timePoint = 0;
        p.avg_PROT_over_timePoint = 0;
        for (p.i=0; p.i<p.i_length; p.i++){
          p.avg_RPF_over_timePoint += p.StoredTEfromRPF[p.i][p.tr];
          p.avg_PROT_over_timePoint += p.StoredTEfromProteomic[p.i][p.tr];
        }//end inner for
        p.StoredTEfromRPFavg[p.tr] = p.avg_RPF_over_timePoint / p.i_length;
        p.StoredTEfromProtavg[p.tr] = p.avg_PROT_over_timePoint / p.i_length;
      }//end outer for
    }//end if
    // final plot scaling preparation
    if (p.StoredTEfromRPF[p.StoredTimes.length - 1][0] > 0 || p.StoredTEfromRPF[p.StoredTimes.length - 1][p.iMin - 1] > 0){
      //p.TEfromRPFMax = p.max(p.StoredTEfromRPF[p.StoredTimes.length-1]);
      //p.TEfromProtMax = p.max(p.StoredTEfromProteomic[p.StoredTimes.length-1]);
      p.TEfromRPFMax = p.max(p.StoredTEfromRPFavg);
      p.TEfromProtMax = p.max(p.StoredTEfromProtavg);
    }
    //---------------------------------------------------------

    //********** plot the TE from proteomic:
    // ... and plot TE from RPFs:
    //--------------------------------------
    p.numbOfbars=24;
    p.iMin = p.min(p.iMax, p.numbOfbars);
    p.W = 20; //width of bars
    p.DX = 40; //x shift to display the RPF TE bars just to the right of the protein TE bars.
    p.DY = 15; //y shift to display coral labels above RPF bars.
    // Each 3 seconds interval, display sequentially the evolution of 36 time points of protein counts
    // display the protein counts
    //p.pyInf-50, p.pySup+35 (pyInf = 580, pySup = 50)
    for (p.i=0;p.i<p.i_length;p.i++){
      //if (p.time_clock/1000 >= (p.i+1)*p.time_interval_to_display && p.timeSeqFlag[p.i] == 'false'){
        //p.timeSeqFlag[p.i] = 'true';
      p.condition = ((p.i == p.i_length-1) || (p.time_clock/1000 >= (p.i+1)*p.time_interval_to_display && p.time_clock/1000 < (p.i+2)*p.time_interval_to_display));
      if (p.condition){
      //if (p.time_clock/1000 >= (p.i+1)*p.time_interval_to_display && p.time_clock/1000 < (p.i+2)*p.time_interval_to_display){
        // display the time point on the screen:
        p.push();
        p.textFont("Helvetica"); //Courier or Helvetica
        p.textSize(20);
        p.fill(0, 90, 90); //(75, 0, 130)
        p.time_text = 'at simulation time '+ p.str(p.StoredTimes[p.i]);
        if (p.i != p.i_length-1){
          p.text(p.time_text, 415, 75);
        }
        p.time_text = '/ total time '+ p.str(p.StoredTimes[p.i_length-1]);
        p.text(p.time_text, 950, 75);
        p.pop();
        for(p.tr=0; p.tr<p.iMin; p.tr++){// scan and display the TE levels of each geneID for the stored time point
          // scale the appropriate height for the bar representing the given TE level
          if (p.i < p.i_length - 1){
            p.yMap = p.map(p.StoredTEfromProteomic[p.i][p.tr], 0, p.TEfromProtMax*1.20, p.pyInf-50, p.pySup+35);
            p.yMapRPF=p.map(p.StoredTEfromRPF[p.i][p.tr], 0, p.TEfromRPFMax*1.25, p.pyInf-50, p.pySup+35);
          }
          if (p.i == p.i_length - 1){
            p.yMap = p.map(p.StoredTEfromProtavg[p.tr], 0, p.TEfromProtMax*1.20, p.pyInf-50, p.pySup+35);
            p.yMapRPF=p.map(p.StoredTEfromRPFavg[p.tr], 0, p.TEfromRPFMax*1.25, p.pyInf-50, p.pySup+35);
          }
          // display the # TE into text:
          p.push();
          p.textFont("Helvetica"); //Courier or Helvetica
          p.textSize(18);
          if (p.i < p.i_length - 1){
            p.fill(100, 105, 100); //(75, 0, 130)
            p.text(p.nfs(p.StoredTEfromProteomic[p.i][p.tr], 2, 1), p.px_tick[p.tr]-41, p.yMap-2);
            p.fill(255, 127, 80); // coral
            p.text(p.nfs(p.StoredTEfromRPF[p.i][p.tr],2, 1), p.px_tick[p.tr]-41+p.DX, p.yMapRPF-2-p.DY);
          }
          if (p.i == p.i_length - 1){
            p.fill(100, 105, 100); //(75, 0, 130)
            p.text(p.nfs(p.StoredTEfromProtavg[p.tr], 2, 1), p.px_tick[p.tr]-41, p.yMap-2);
            p.fill(255, 127, 80); // coral
            p.text(p.nfs(p.StoredTEfromRPFavg[p.tr],2, 1), p.px_tick[p.tr]-41+p.DX, p.yMapRPF-2-p.DY);
          }

          p.pop();

          //display bars as rectangle:
          //...
          p.push();
          p.fill(0, 90, 90, 125);
          p.strokeWeight(2);
          p.stroke(0, 90, 90, 125); // transparency is 100 between 0-250 (250 is opaque).
          // comment out the following if you don't want to display the TE from proteomic:
          p.rect(p.px_tick[p.tr]-10, p.pyInf-50, p.W-2, p.yMap-p.pyInf+50);
          p.pop();
          p.push();
          p.fill(255, 127, 80, 125); // coral bars:
          p.strokeWeight(2);
          p.stroke(255, 127, 80, 200); // coral bars:
          p.rect(p.px_tick[p.tr]+10, p.pyInf-50, p.W-2, p.yMapRPF-p.pyInf+50);
          p.pop();

          // add a tick above the bar or rectangle:
          p.push();
          p.strokeWeight(2);
          p.stroke(150); // stroke(100, 0, 170)
          p.ticklength = 7;
          p.line(p.px_tick[p.tr], p.yMap - 2 * p.ticklength, p.px_tick[p.tr], p.yMap);
          p.stroke(255, 127, 80, 200); // coral bars:
          p.line(p.px_tick[p.tr]+20, p.yMapRPF - 2 * p.ticklength, p.px_tick[p.tr]+20, p.yMapRPF); // tick above coral bar
          p.pop();

          //

          // draw ticks :
          p.py_tick_max = p.map(p.TEfromProtMax, 0, p.TEfromProtMax*1.20, p.pyInf-50, p.pySup+35);
          //p.py_tick_max_RPF = p.map(p.TEfromRPFMax, 0, p.TEfromRPFMax*1.20, p.pyInf-50, p.pySup+35);
          // add a dashed line from the top tick to the geneID labels
          p.push();
          p.strokeWeight(2);
          p.stroke(150, 125); // stroke(100, 0, 170)
          p.drawingContext.setLineDash([8, 16]);
          p.line(p.px_tick[p.tr], p.yMap - 2 * p.ticklength, p.px_tick[p.tr], p.py_tick_max - 1*p.ticklength );
          p.pop();
        }
      }
    }

    // display their labels above the bars:
    for (p.i=0; p.i<p.iMin;p.i++){
      p.push();
      //p.translate(p.px_tick[p.i], p.pyInf - 50 - 2*p.ticklength);
      p.translate(p.px_tick[p.i], p.py_tick_max - 2*p.ticklength);
      p.rotate(p.radians(-50));
      p.fill(100, 105, 100);
      p.text(p.StoredGeneID[p.i], 0, 0);
      //p.text(p.StoredGeneID[p.i], p.px_tick[p.i], p.pyInf - 50 - 2*p.ticklength);
      p.pop();
    }


    //********* plot TE from RPF: add coral ordinate...
    //---------------------------

    //*******************************************************************
    // draw the border of the graphic box:
    p.push();
    p.strokeWeight(4);
    p.stroke(0, 90, 90);
    p.line(p.pxInf, p.pyInf+p.pySHIFT, p.pxSup, p.pyInf+p.pySHIFT);
    p.line(p.pxInf, p.pyInf+p.pySHIFT, p.pxInf, p.pySup+p.pySHIFT);

    p.line(p.pxInf, p.pySup+p.pySHIFT, p.pxSup, p.pySup+p.pySHIFT);
    p.line(p.pxSup, p.pyInf+p.pySHIFT, p.pxSup, p.pySup+p.pySHIFT);
    p.pop();
    // draw the abscissa and ordinate
    p.push();
    p.strokeWeight(2);
    p.stroke(150);
    p.line(p.pxInf+60, p.pyInf-50+p.pySHIFT, p.pxSup-30, p.pyInf-50+p.pySHIFT); // abscissa
    p.line(p.pxInf+60, p.pyInf-50+p.pySHIFT, p.pxInf+60, p.pySup+35+p.pySHIFT); //ordinate
    p.pop();

    // graphic box title and axes labels
    p.lineThree = "Transcripts abundance: reads copy numbers (transcripts counts during translation) --- RNA-Seq within sample";
    p.textFont("Helvetica"); //Courier or Helvetica or Calibri
    p.textSize(24);
    p.fill(0, 90, 90); //(75, 0, 130) // coral color 255, 127, 80
    p.text(p.lineThree, 40, p.pyInf + 60);

    // axes labels:
    p.ord_label_two = "transcript count";
    p.textFont("Helvetica"); //Courier or Helvetica
    p.textSize(20);
    p.fill(100, 105, 100); //(75, 0, 130)
    p.text(p.ord_label_two, p.pxInf+12, p.pySup+25+p.pySHIFT);

    p.abs_label_two = "transcript ID";
    p.textFont("Helvetica"); //Courier or Helvetica
    p.textSize(20);
    p.fill(100, 105, 100); //(75, 0, 130)
    p.text(p.abs_label_two, p.pxInf + 700, p.pyInf-25+p.pySHIFT);

    // find max of reads counts:
    p.protTransMax = p.max(p.StoredTransCount);

    // draw ticks :
    p.py_tick_max = p.map(p.protTransMax, 0, p.protTransMax*1.20, p.pyInf-50+p.pySHIFT, p.pySup+35+p.pySHIFT);
    //console.log('py_tick_max=', p.py_tick_max);
    p.push();
    p.strokeWeight(2);
    p.stroke(150); // stroke(100, 0, 170)
    p.ticklength = 7;
    for (p.i=0; p.i<24; p.i++){
      p.line(p.px_tick[p.i], p.pyInf - 50 - p.ticklength+p.pySHIFT, p.px_tick[p.i], p.pyInf - 50+p.pySHIFT);
    }
    p.line(p.pxInf + 60, p.py_tick_max, p.pxInf + 60 + p.ticklength, p.py_tick_max);
    p.pop();

    p.numbOfbars=24;
    p.iMin = p.min(p.iMax, p.numbOfbars);
    for(p.tr=0; p.tr<p.iMin; p.tr++){// scan and display the read counts of each geneID
      // scale the appropriate height for the bar representing the given transcripts counts
      p.yMap = p.map(p.StoredTransCount[p.tr], 0, p.protTransMax*1.20, p.pyInf-50+p.pySHIFT, p.pySup+35+p.pySHIFT);
      // display the # reads count into text:
      p.push();
      p.textFont("Helvetica"); //Courier or Helvetica
      p.textSize(20);
      p.fill(100, 105, 100); //(75, 0, 130)
      p.text(p.StoredTransCount[p.tr], p.px_tick[p.tr]+11, p.yMap-2);
      p.pop();

      //display bars as rectangle:
      p.push();
      p.fill(90, 90, 90, 125);
      p.strokeWeight(2);
      p.stroke(0, 0, 0, 125); // transparency is 100 between 0-250 (250 is opaque).
      p.rect(p.px_tick[p.tr]-10, p.pyInf-50+p.pySHIFT, 20, p.yMap-p.pyInf+50-p.pySHIFT);
      p.pop();

      // add a tick above the bar or rectangle:
      p.push();
      p.strokeWeight(2);
      p.stroke(150); // stroke(100, 0, 170)
      p.ticklength = 7;
      p.line(p.px_tick[p.tr], p.yMap - 2 * p.ticklength, p.px_tick[p.tr], p.yMap);
      p.pop();

      // add a dashed line from the top tick to the geneID labels
      p.push();
      p.strokeWeight(2);
      p.stroke(150, 125); // stroke(100, 0, 170)
      p.drawingContext.setLineDash([8, 16]);
      p.line(p.px_tick[p.tr], p.yMap - 2 * p.ticklength, p.px_tick[p.tr], p.py_tick_max - 1*p.ticklength);
      p.pop();

      // display their labels above the bars:
      p.push();
      //p.translate(p.px_tick[p.i], p.pyInf - 50 - 2*p.ticklength);
      p.translate(p.px_tick[p.tr], p.py_tick_max - 2*p.ticklength);
      p.rotate(p.radians(-50));
      p.text(p.StoredGeneID[p.tr], 0, 0);
      p.pop();

    }// end of for loop

  }//end of p.draw
}// end var sketch

var myProteinsP5 = new p5(sketch1); // this is a p5 object called myProteinsP5
// this p5 object will have all the variables associated to it and
//the setup or draw functions associated to it.
var myTEp5 = new p5(sketch2);
