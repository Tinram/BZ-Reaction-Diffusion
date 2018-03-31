
/**
    * Belousov-Zhabotinsky Reaction Diffusion
    *
    * Converted from Eva Schindling's C++
    *
    * 'p' - process colours
    * 's' - save image
    *
    * @author      Martin Latter
    * @copyright   Martin Latter 26/06/2013
    * @version     4.01, updated to Processing 3.3.6
    * @license     The MIT License (MIT)
    * @link        https://github.com/Tinram/BZ-Reaction-Diffusion.git
*/


/* CONFIGURATION */
int iWidth = 1000, iHeight = 1000;
int iFinalCol = 255; // 255 or 0



boolean bEnd, bEnd2 = false;


float fDiffRateA = 0.34521;          // diffusion rate for chemical A - default: 0.35;
float fDiffRateB = 0.36453;          // diffusion rate for chemical B - default: 0.37;
float fDecayNoise = 0.05;            // noise in decay parameter - default: 0.1
float fDecaySeed = 1;                // default: 1

float fChemAStart = random(40000);   // start concentration for chemical A
float fChemBStart = random(70000);   // start concentration for chemical B

float fMinChemVal = 0.0;             // minimum value for chemical
float fMaxChemVal = 10.0;            // maximum value for chemical

float fGrowth = 16;
float fGrowthX = 0.0;                // vary growth over x-axis

float fChemADiffRateX = 0.15;        // vary diffusion parameter for chemical A over x-axis
float fChemBDiffRateY = 0.5;         // vary diffusion parameter for chemical B over y-axis

float fReactionSpeed = 0.028;        // reaction speed
float fReactionSpeedY = 0.0;         // vary reaction fReactionSpeed over y-axis

float fDecayMean = 12.2;             // mean decay parameter
float fNoiseY = 0.0;                 // vary noise over y-axis

float fFrameTimeStep = 3.0;          // timeStep per frame
float fElapsedTime = 0.0;            // elapsed time

int iUpdatesPerFrame = 1;            // updates to perform per frame


float[][] aChemA;
float[][] aChemB;
float[][] aChemADeriv;
float[][] aChemBDeriv;
float[][] aDecay;


void setup() {

  size(1000, 1000); // iWidth, iHeight

  /* create arrays */
  aChemA = new float[iWidth][iHeight];
  aChemB = new float[iWidth][iHeight];
  aDecay = new float[iWidth][iHeight];
  aChemADeriv = new float[iWidth][iHeight];
  aChemBDeriv = new float[iWidth][iHeight];

  /* initialise arrays */
  for (int x=0; x < iWidth; x++) {
    for (int y = 0; y < iHeight; y++) {
      aChemA[x][y] = fChemAStart;
      aChemB[x][y] = fChemBStart;
      aDecay[x][y] = fDecayMean + fDecayNoise * random(fDecaySeed);
    }
  }
}


void draw() {

  color c;
  float fConc, fConc2, fMix;

  updateChemical();

  loadPixels();

  for (int y = 0; y < iHeight; y++) {
    for (int x = 0; x < iWidth; x++) {
      fConc = (255 * (aChemB[y][x] / fMaxChemVal));
      fConc2 = (255 * (aChemA[y][x] / fMaxChemVal));
      fMix = (fConc2 - fConc);
      c = color(fConc, fMix, fConc2); //  c = color(200, fConc, 200);
      pixels[y * iWidth + x] = c;
    }
  }

  updatePixels();

  // staggered end clear one chemical
  if (bEnd) {

    if (bEnd2) {
      noLoop();
    }

    bEnd2 = true;
  }
}


void updateChemical() {

  float fChemComb, fGrowthComb, fDecayRate, fDiffX, fDiffY, fSpeed;
  int w = iWidth, h = iHeight;

  for (int i = 0; i < iUpdatesPerFrame; i++) {

    fElapsedTime += fFrameTimeStep;

    for (int y = 0; y < h; y++) {

      for (int x = 0; x < w; x++) {

        fChemComb = aChemA[y][x] * aChemB[y][x];
        fGrowthComb = fGrowth + fGrowthX * (x / iWidth);
        fDecayRate = (fNoiseY > 0)? fDecayMean + fNoiseY * (y / iHeight) : aDecay[y][x];
        fDiffX = fDiffRateA + (x / w) * fChemADiffRateX;
        fDiffY = fDiffRateB + (y / h) * fChemBDiffRateY;
        fSpeed = fReactionSpeed + fReactionSpeedY * (y / iHeight);
        aChemADeriv[y][x] = fSpeed * (fGrowthComb - fChemComb) + fDiffX * calcDiffusion(aChemA, x, y);
        aChemBDeriv[y][x] = fSpeed * (fChemComb - aChemB[y][x] - fDecayRate) + fDiffY * calcDiffusion(aChemB, x, y);
      }
    }

    for (int y = 0; y < iHeight; y++) {

      for (int x = 0; x < iWidth; x++) {

        /* update chemicals and constrain concentration to min + max) */
        aChemA[y][x] = constrain(aChemA[y][x] + fFrameTimeStep * aChemADeriv[y][x], fMinChemVal, fMaxChemVal);
        aChemB[y][x] = constrain(aChemB[y][x] + fFrameTimeStep * aChemBDeriv[y][x], fMinChemVal, fMaxChemVal);
      }
    }
  }
}


float calcDiffusion(float[][] aChemical, int x, int y) {

  /* diffusion calculator */

  if (x > 0 && x < iWidth - 1 && y > 0 && y < iHeight - 1) {
    float fVal = aChemical[y][x - 1] + aChemical[y][x + 1] + aChemical[y - 1][x] + aChemical[y + 1][x];
    return -(aChemical[y][x] - fVal / 4);
  } else {
    return 0.0;
  }
}


void keyPressed() {

  if (key == 's') {
    save("BZ_Reaction_Diffusion_" + millis() + ".png");
  }

  if (key == 'p') {

    // make one chemical blank
    for (int x=0; x < iWidth; x++) {
      for (int y = 0; y < iHeight; y++) {
        aChemA[x][y] = iFinalCol;
      }
    }

    bEnd = true;
  }
}
