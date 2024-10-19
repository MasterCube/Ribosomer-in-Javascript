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
