

void trig_eff_val(){
    float binshiftedpt[20] = {5.72888, 6.23135, 6.7329, 7.23396, 7.7352, 
                            8.23592, 8.73631, 9.23676, 9.73706, 10.8236, 12.8478, 14.876, 16.8769, 18.8868,
                            20.9065, 22.9239, 24.9232, 27.7391, 31.7704, 37.1476};
    TF1 *doubleErfFit = new TF1("sigmoidFit_Sc2", "[0] / (1 + exp(-(x - [1]) / [2]))", 2, 15);
    // GAMMA3
    // doubleErfFit->SetParameters(9.58075e-01, 3.75986e+00, 7.41923e-01); // PbSc
    // doubleErfFit->SetParameters(6.71915e-01, 5.31439e+0, 1.01171e+00 ); // PbGl

    // MINBIAS
    doubleErfFit->SetParameters( 5.24635e-01,  5.25475e+00, 7.89932e-01); // PbGl
    for (int i = 0; i < 20; i++){
        std::cout << binshiftedpt[i] << ", " << doubleErfFit->Eval(binshiftedpt[i]) << std::endl;
    }

    /*
      EXT PARAMETER                                   STEP         FIRST   
      NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
      1  p0           9.58075e-01   8.87647e-03   4.16008e-06  -9.38784e-04
      2  p1           3.75986e+00   2.23323e-01   1.09010e-04   1.84363e-04
      3  p2           7.41923e-01   1.27600e-01   5.26681e-05   1.03341e-04
      FCN=1.78873 FROM MIGRAD    STATUS=CONVERGED     157 CALLS         158 TOTAL
      EDM=2.43656e-07    STRATEGY= 1      ERROR MATRIX ACCURATE 
      EXT PARAMETER                                   STEP         FIRST   
      NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
      1  p0           5.24635e-01   1.36441e-01   5.30124e-05   6.25605e-03
      2  p1           5.25475e+00   5.63642e-01   1.57848e-04  -7.26412e-04
      3  p2           7.89932e-01   2.10229e-01   9.11160e-05   1.36106e-03
      FCN=3.61878 FROM MIGRAD    STATUS=CONVERGED     101 CALLS         102 TOTAL
      EDM=8.76755e-08    STRATEGY= 1      ERROR MATRIX ACCURATE 
      EXT PARAMETER                                   STEP         FIRST   
      NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
      1  p0           6.71915e-01   1.35839e-02   1.03028e-05  -3.09485e-02
      2  p1           5.31439e+00   1.74919e-01   1.58438e-04   2.38564e-03
      3  p2           1.01171e+00   1.15954e-01   9.85904e-05   3.76833e-04
      Info in <TCanvas::Print>: pdf file gamma3_Comb_TrigEff.pdf has been created
    */

}
