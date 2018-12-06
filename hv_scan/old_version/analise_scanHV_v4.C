// algorithm to search the best HV operational point for Modulation experimento in setup CBPF
// 2 criteria: delta Energy and signal/background ratio
// list of points HV (Perci please check the HV input values)

// Massa        Perci
// 0=660V        660V ok
// 1=640V        660V ok
// 2=880V        880V ok 
// 3=760V        760V ok
// 4=720V        720V ok
// 5=450V ok     480V 
// 6=580V        560V ok
// 7=600V        580V ok

Double_t fun_Back(Double_t *x, Double_t *par)
{
    Double_t offset, B;

    offset = x[0];
    B      = par[0]*TMath::Exp(-par[1]*offset);

    return B;
}

Double_t fun_Signal_Back(Double_t *x, Double_t *par)
{
    Double_t S,B;
    Double_t cte = 2.5066297;
    Double_t Norm, offset;

    Norm = 1.0/(cte*par[2]);
    S    = par[0]*Norm*(TMath::Exp(-(par[1]-x[0])*(par[1]-x[0])/(2.0*par[2]*par[2])));  
 
    offset = x[0];
    B      = par[3]*TMath::Exp(-par[4]*offset);

	return S+B;
}

Double_t fun_twoSignal_Back(Double_t *x, Double_t *par)
{
    Double_t S1,S2,B;
    Double_t cte = 2.5066297;
    Double_t Norm1, Norm2;
    Double_t offset;

    Norm1 = 1.0/(cte*par[2]);
    S1    = par[0]*Norm1*(TMath::Exp(-(par[1]-x[0])*(par[1]-x[0])/(2.0*par[2]*par[2])));  
    
    Norm2 = 1.0/(cte*par[5]);
    S2    = par[3]*Norm2*(TMath::Exp(-(par[4]-x[0])*(par[4]-x[0])/(2.0*par[5]*par[5])));  
 
    offset = x[0];
    B      = par[6]*TMath::Exp(-par[7]*offset);

    return S1+S2+B;
}

Double_t fun_threeSignal_Back(Double_t *x, Double_t *par)
{
    Double_t S1,S2,S3,B;
    Double_t cte = 2.5066297;
    Double_t Norm1, Norm2, Norm3, offset;

    Norm1 = 1.0/(cte*par[2]);
    S1    = par[0]*Norm1*(TMath::Exp(-(par[1]-x[0])*(par[1]-x[0])/(2.0*par[2]*par[2])));  
    
    Norm2 = 1.0/(cte*par[5]);
    S2    = par[3]*Norm2*(TMath::Exp(-(par[4]-x[0])*(par[4]-x[0])/(2.0*par[5]*par[5])));  
    
    Norm3 = 1.0/(cte*par[8]);
    S3    = par[6]*Norm3*(TMath::Exp(-(par[7]-x[0])*(par[7]-x[0])/(2.0*par[8]*par[8])));  
 
    offset = x[0];
    B      = par[9]*TMath::Exp(-par[10]*offset);

    return S1+S2+S3+B;
}


Float_t ErrDiv(Float_t a, Float_t da, Float_t b, Float_t db){

    return (1/b)*(1/b)*sqrt((b*b*da*da) + (a*a*db*db));
}

void result(int ichannel=0, int nbin  = 1000, int nbin_fun = 200) 
{
   
    int model = 0;
    TString source_name= "";

    if ((ichannel==0) || (ichannel==1) || (ichannel==4) || (ichannel==5) )
    {
        model = 1;
        source_name="{}^{60}Co"; 
    } else {
        model = 0;
        source_name="{}^{137}Cs"; 
    }
    

    ////////// detector 0 : Co 
    /*

    const int dim = 10; 
    Char_t namedata[dim][200] = {"HVCalibration/Grupo1-2/processed/mx_b_20170519_2044/calibration/mx_b_20170519_2044_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170519_2101/calibration/mx_b_20170519_2101_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170519_2119/calibration/mx_b_20170519_2119_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170519_2135/calibration/mx_b_20170519_2135_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170519_2149/calibration/mx_b_20170519_2149_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_1922/calibration/mx_b_20170523_1922_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_1958/calibration/mx_b_20170523_1958_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_2016/calibration/mx_b_20170523_2016_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_2042/calibration/mx_b_20170523_2042_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_2059/calibration/mx_b_20170523_2059_000000.root"};

    Char_t namedata_back[dim][300] = {"HVCalibration/Background/processed/mx_b_20170606_1922/calibration/mx_b_20170606_1922_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_1950/calibration/mx_b_20170606_1950_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2045/calibration/mx_b_20170606_2045_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2111/calibration/mx_b_20170606_2111_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2146/calibration/mx_b_20170606_2146_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2203/calibration/mx_b_20170606_2203_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2221/calibration/mx_b_20170606_2221_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2240/calibration/mx_b_20170606_2240_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2256/calibration/mx_b_20170606_2256_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2313/calibration/mx_b_20170606_2313_000000.root"};


    
    float hv[dim]   = {520,540,560,580,600,640,660,680,700,720};
    TString name_latex3 = "HV_{scionix} = 620 V";
    float minbin[dim] = {0,0,0,0,0,0,0,0,0,0};
    float maxbin[dim] = {100,120,160,200,300,500,800,800,800,800};

    const int dim_use = 10; 
    int useful[dim_use] = {0,1,2,3,4,5,6,7,8,9};

    float minfun[dim_use] = {30,42,57,77,100,170,220,280,360,460};
    float maxfun[dim_use] = {45,62,82,105,142,235,305,383,490,635};

    /// gauss1:N1,media,width, gauss2:N2,media,width, back:N3,tail(exponential)
    float parinit0[dim_use] = {1000,20000,19000,15000,14000,10000,10000,10000,10000,10000};
    float parinit1[dim_use] = {35,47,65,85,114,195,250,317,400,500};
    float parinit2[dim_use] = {0.9,2,2,2,2,3,3,4,4,4};
    float parinit3[dim_use] = {1000,20000,19000,15000,10000,10000,10000,10000,10000,10000};
    float parinit4[dim_use] = {41,55,74,98,130,220,280,360,450,570};
    float parinit5[dim_use] = {0.9,2,2,2,2,3,3,4,4,4};
    float parinit6[dim_use] = {200,1000,1000,1000,1000,1000,1000,1000,1000,1000};
    float parinit7[dim_use] = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001};
 
    */







    //// detector 1 : Co 
    /*
    const int dim = 10; 
    Char_t namedata[dim][200] = {"HVCalibration/Grupo1-2/processed/mx_b_20170519_2044/calibration/mx_b_20170519_2044_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170519_2101/calibration/mx_b_20170519_2101_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170519_2119/calibration/mx_b_20170519_2119_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170519_2135/calibration/mx_b_20170519_2135_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170519_2149/calibration/mx_b_20170519_2149_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_1922/calibration/mx_b_20170523_1922_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_1958/calibration/mx_b_20170523_1958_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_2016/calibration/mx_b_20170523_2016_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_2042/calibration/mx_b_20170523_2042_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_2059/calibration/mx_b_20170523_2059_000000.root"};

    Char_t namedata_back[dim][300] = {"HVCalibration/Background/processed/mx_b_20170606_1922/calibration/mx_b_20170606_1922_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_1950/calibration/mx_b_20170606_1950_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2045/calibration/mx_b_20170606_2045_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2111/calibration/mx_b_20170606_2111_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2146/calibration/mx_b_20170606_2146_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2203/calibration/mx_b_20170606_2203_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2221/calibration/mx_b_20170606_2221_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2240/calibration/mx_b_20170606_2240_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2256/calibration/mx_b_20170606_2256_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2313/calibration/mx_b_20170606_2313_000000.root"};

    float hv[dim]     = {540,560,580,600,620,660,680,700,720,740};
    TString name_latex3 = "HV_{scionix} = 640 V";
    float minbin[dim] = {0,0,0,0,0,0,0,0,0,0};
    float maxbin[dim] = {100,120,160,200,300,500,800,800,800,800};

    const int dim_use = 9; 
    int useful[dim_use] = {0,1,2,3,4,5,6,7,8};

    float minfun[dim_use] = {40,55,77,102,136,228,290,360,460};
    float maxfun[dim_use] = {57,77,102,136,177,300,380,480,610};

    float parinit0[dim_use] = {2000,2000,2000,2000,2000,2000,2000,2000,2000};
    float parinit1[dim_use] = {47,63,85,113,148,250,318,400,505};
    float parinit2[dim_use] = {1,2,2,2,2,3,3,4,4};
    float parinit3[dim_use] = {2000,2000,2000,2000,2000,2000,2000,2000,2000};
    float parinit4[dim_use] = {53,72,97,128,168,280,360,455,570};
    float parinit5[dim_use] = {1,2,2,2,2,3,3,4,4};
    float parinit6[dim_use] = {200,1000,1000,1000,1000,1000,1000,1000,1000};
    float parinit7[dim_use] = {0.001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001};
 
    */







    // detector 2 : Cs 
    // using 50 bins for fit
    /*
    const int dim = 10; 
    Char_t namedata[dim][200] = {"HVCalibration/Grupo1-2/processed/mx_b_20170519_2044/calibration/mx_b_20170519_2044_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170519_2101/calibration/mx_b_20170519_2101_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170519_2119/calibration/mx_b_20170519_2119_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170519_2135/calibration/mx_b_20170519_2135_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170519_2149/calibration/mx_b_20170519_2149_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_1922/calibration/mx_b_20170523_1922_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_1958/calibration/mx_b_20170523_1958_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_2016/calibration/mx_b_20170523_2016_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_2042/calibration/mx_b_20170523_2042_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_2059/calibration/mx_b_20170523_2059_000000.root"};

    Char_t namedata_back[dim][300] = {"HVCalibration/Background/processed/mx_b_20170606_1922/calibration/mx_b_20170606_1922_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_1950/calibration/mx_b_20170606_1950_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2045/calibration/mx_b_20170606_2045_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2111/calibration/mx_b_20170606_2111_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2146/calibration/mx_b_20170606_2146_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2203/calibration/mx_b_20170606_2203_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2221/calibration/mx_b_20170606_2221_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2240/calibration/mx_b_20170606_2240_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2256/calibration/mx_b_20170606_2256_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2313/calibration/mx_b_20170606_2313_000000.root"};

    
    float hv[dim]     = {700,720,740,760,780,820,840,860,880,900};
    TString name_latex3 = "HV_{scionix} = 800 V";
    float minbin[dim] = {0,0,0,0,0,0,0,0,0,0};
    float maxbin[dim] = {90,90,90,90,90,90,90,90,90,90};

    const int dim_use = 6; 
    int useful[dim_use] = {4,5,6,7,8,9};

    float minfun[dim_use] = {16,24,30,36,44,52};
    float maxfun[dim_use] = {24,35,44,50,60,72};

    float parinit0[dim_use] = {4000,4000,5000,4000,4000,4000};
    float parinit1[dim_use] = {20,29,38,42,50,62};
    float parinit2[dim_use] = {1,1,1.5,1,1,1};
    float parinit3[dim_use] = {100,100,500,100,100,100};
    float parinit4[dim_use] = {0.0001,0.0001,0.000001,0.0001,0.0001,0.0001};
    
    float parinit5[dim_use] = {0,0,0,0,0,0};
    float parinit6[dim_use] = {1,1,1,1,1,1};
    float parinit7[dim_use] = {1,1,1,1,1,1};

    */





    // detector 3 : Cs 
    // /*
    const int dim = 9;
    Char_t namedata[dim][200] = {"HVCalibration/Grupo1-2/processed/mx_b_20170519_2101/calibration/mx_b_20170519_2101_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170519_2119/calibration/mx_b_20170519_2119_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170519_2135/calibration/mx_b_20170519_2135_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170519_2149/calibration/mx_b_20170519_2149_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_1922/calibration/mx_b_20170523_1922_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_1958/calibration/mx_b_20170523_1958_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_2016/calibration/mx_b_20170523_2016_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_2042/calibration/mx_b_20170523_2042_000000.root",
                                 "HVCalibration/Grupo1-2/processed/mx_b_20170523_2059/calibration/mx_b_20170523_2059_000000.root"};
    Char_t namedata_back[dim][300] = {//"HVCalibration/Background/processed/mx_b_20170606_1922/calibration/mx_b_20170606_1922_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_1950/calibration/mx_b_20170606_1950_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2045/calibration/mx_b_20170606_2045_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2111/calibration/mx_b_20170606_2111_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2128/calibration/mx_b_20170606_2128_000000.root",
                                      //"HVCalibration/Background/processed/mx_b_20170606_2146/calibration/mx_b_20170606_2146_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2203/calibration/mx_b_20170606_2203_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2221/calibration/mx_b_20170606_2221_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2240/calibration/mx_b_20170606_2240_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2256/calibration/mx_b_20170606_2256_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2313/calibration/mx_b_20170606_2313_000000.root"};
 
    float hv[dim]     = {620.0,640.0,660.0,680.0,720.0,740.0,760.0,780,800.0};
    TString name_latex3 = "HV_{scionix} = 700 V";
    float minbin[dim] = {0,0,0,0,0,0,0,0,0};
    float maxbin[dim] = {40,50,60,80,140,160,200,250,300};
    
    const int dim_use = 8; 
    int useful[dim_use] = {1,2,3,4,5,6,7,8};
    
    float minfun[dim_use] = {16,20,26,42,52,65,75,95};
    float maxfun[dim_use] = {24,30,38,60,75,90,115,140};
    
    float parinit0[dim_use] = {300,500,500,500,500,500,500,500};
    float parinit1[dim_use] = {18,25,31,50,62,76,95,118};
    float parinit2[dim_use] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
    float parinit3[dim_use] = {100,100,100,100,100,100,100,100};
    float parinit4[dim_use] = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001};
    float parinit5[dim_use],parinit6[dim_use],parinit7[dim_use];
  
    // */






    // detector 4 : Co 
    /*

    // background sample is anomalous (too small) for this detector

    const int dim = 11;
    Char_t namedata[dim][200] = {
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_1956/calibration/mx_b_20170524_1956_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2013/calibration/mx_b_20170524_2013_000000.root",   
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2043/calibration/mx_b_20170524_2043_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2059/calibration/mx_b_20170524_2059_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2130/calibration/mx_b_20170524_2130_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2152/calibration/mx_b_20170524_2152_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2215/calibration/mx_b_20170524_2215_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2235/calibration/mx_b_20170524_2235_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2249/calibration/mx_b_20170524_2249_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2309/calibration/mx_b_20170524_2309_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2329/calibration/mx_b_20170524_2329_000000.root"};

    Char_t namedata_back[dim][300] = {"HVCalibration/Background/processed/mx_b_20170606_1922/calibration/mx_b_20170606_1922_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_1950/calibration/mx_b_20170606_1950_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2045/calibration/mx_b_20170606_2045_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2111/calibration/mx_b_20170606_2111_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2128/calibration/mx_b_20170606_2128_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2146/calibration/mx_b_20170606_2146_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2203/calibration/mx_b_20170606_2203_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2221/calibration/mx_b_20170606_2221_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2240/calibration/mx_b_20170606_2240_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2256/calibration/mx_b_20170606_2256_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2313/calibration/mx_b_20170606_2313_000000.root"};
 
    float hv[dim]     = {560,580,600,620,640,660,680,700,720,740,760};
    TString name_latex3 = "HV_{scionix} = 660 V";
    float minbin[dim] = {0,0,0,0,0,0,0,0,0,0,0};
    float maxbin[dim] = {30,30,30,30,30,50,70,100,100,140,200};

    const int dim_use = 4; 
    int useful[dim_use] = {6,7,8,9};

    float minfun[dim_use] = {27,36,38,52};
    float maxfun[dim_use] = {40,50,52,70};

    float parinit0[dim_use] = {500,500,500,500};
    float parinit1[dim_use] = {31,40,42,57};
    float parinit2[dim_use] = {0.5,1,1,1};
    float parinit3[dim_use] = {500,500,500,500};
    float parinit4[dim_use] = {35,46,47,65};
    float parinit5[dim_use] = {0.5,1,1,1};
    float parinit6[dim_use] = {40,20,20,20};
    float parinit7[dim_use] = {0.0001,0.0001,0.0001,0.0001};
 
    */






    // detector 5 : Co 

    /*
    const int dim = 11;
    Char_t namedata[dim][200] = {
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_1956/calibration/mx_b_20170524_1956_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2013/calibration/mx_b_20170524_2013_000000.root",   
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2043/calibration/mx_b_20170524_2043_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2059/calibration/mx_b_20170524_2059_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2130/calibration/mx_b_20170524_2130_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2152/calibration/mx_b_20170524_2152_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2215/calibration/mx_b_20170524_2215_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2235/calibration/mx_b_20170524_2235_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2249/calibration/mx_b_20170524_2249_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2309/calibration/mx_b_20170524_2309_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2329/calibration/mx_b_20170524_2329_000000.root"};

    Char_t namedata_back[dim][300] = {"HVCalibration/Background/processed/mx_b_20170606_1922/calibration/mx_b_20170606_1922_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_1950/calibration/mx_b_20170606_1950_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2045/calibration/mx_b_20170606_2045_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2111/calibration/mx_b_20170606_2111_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2128/calibration/mx_b_20170606_2128_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2146/calibration/mx_b_20170606_2146_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2203/calibration/mx_b_20170606_2203_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2221/calibration/mx_b_20170606_2221_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2240/calibration/mx_b_20170606_2240_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2256/calibration/mx_b_20170606_2256_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2313/calibration/mx_b_20170606_2313_000000.root"};
 
    float hv[dim]     = {370,390,410,430,450,470,490,510,530,550,570};
    TString name_latex3 = "HV_{scionix} = 470 V";
    float minbin[dim] = {0,0,0,0,0,0,0,0,0,0,0};
    float maxbin[dim] = {30,30,30,50,50,50,50,50,50,50,50};

    const int dim_use = 3; 
    int useful[dim_use] = {3,4,5};

    float minfun[dim_use] = {13,20,30};
    float maxfun[dim_use] = {21,29,40};

    float parinit0[dim_use] = {1000,2000,2000,};
    float parinit1[dim_use] = {15,23,33};
    float parinit2[dim_use] = {1,1,1};
    float parinit3[dim_use] = {2000,2000,2000};
    float parinit4[dim_use] = {18,26,37};
    float parinit5[dim_use] = {1,1,1};
    float parinit6[dim_use] = {1000,1000,1000};
    float parinit7[dim_use] = {0.0001,0.0001,0.0001};
 
    */







    // detector 6 : Cs 
    /*
    const int dim = 11;
    Char_t namedata[dim][200] = {
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_1956/calibration/mx_b_20170524_1956_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2013/calibration/mx_b_20170524_2013_000000.root",   
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2043/calibration/mx_b_20170524_2043_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2059/calibration/mx_b_20170524_2059_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2130/calibration/mx_b_20170524_2130_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2152/calibration/mx_b_20170524_2152_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2215/calibration/mx_b_20170524_2215_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2235/calibration/mx_b_20170524_2235_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2249/calibration/mx_b_20170524_2249_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2309/calibration/mx_b_20170524_2309_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2329/calibration/mx_b_20170524_2329_000000.root"};

    Char_t namedata_back[dim][300] = {"HVCalibration/Background/processed/mx_b_20170606_1922/calibration/mx_b_20170606_1922_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_1950/calibration/mx_b_20170606_1950_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2045/calibration/mx_b_20170606_2045_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2111/calibration/mx_b_20170606_2111_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2128/calibration/mx_b_20170606_2128_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2146/calibration/mx_b_20170606_2146_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2203/calibration/mx_b_20170606_2203_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2221/calibration/mx_b_20170606_2221_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2240/calibration/mx_b_20170606_2240_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2256/calibration/mx_b_20170606_2256_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2313/calibration/mx_b_20170606_2313_000000.root"};
 
    float hv[dim]     = {380,400,420,440,460,480,500,520,540,560,580};
    TString name_latex3 = "HV_{scionix} = 480 V";
    float minbin[dim] = {0,0,0,0,0,0,0,0,0,0,0};
    float maxbin[dim] = {30,30,30,40,40,50,70,100,150,200,300};

    const int dim_use = 6; 
    int useful[dim_use] = {5,6,7,8,9,10};
    float minfun[dim_use] = {15,22,32,52,70,100};
    float maxfun[dim_use] = {24,34,44,70,95,125};

    float parinit0[dim_use] = {500,500,500,500,500,500};
    float parinit1[dim_use] = {19,27,38,62,84,112};
    float parinit2[dim_use] = {0.5,0.5,0.5,0.5,0.5,0.5};
    float parinit3[dim_use] = {100,100,100,100,100,100};
    float parinit4[dim_use] = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001};
    float parinit5[dim_use],parinit6[dim_use],parinit7[dim_use];
  
    */






    // detector 7 : Cs 
    /*
    const int dim = 11;
    Char_t namedata[dim][200] = {
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_1956/calibration/mx_b_20170524_1956_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2013/calibration/mx_b_20170524_2013_000000.root",   
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2043/calibration/mx_b_20170524_2043_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2059/calibration/mx_b_20170524_2059_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2130/calibration/mx_b_20170524_2130_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2152/calibration/mx_b_20170524_2152_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2215/calibration/mx_b_20170524_2215_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2235/calibration/mx_b_20170524_2235_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2249/calibration/mx_b_20170524_2249_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2309/calibration/mx_b_20170524_2309_000000.root",
                                 "HVCalibration/Grupo3-4/processed/mx_b_20170524_2329/calibration/mx_b_20170524_2329_000000.root"};

    Char_t namedata_back[dim][300] = {"HVCalibration/Background/processed/mx_b_20170606_1922/calibration/mx_b_20170606_1922_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_1950/calibration/mx_b_20170606_1950_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2045/calibration/mx_b_20170606_2045_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2111/calibration/mx_b_20170606_2111_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2128/calibration/mx_b_20170606_2128_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2146/calibration/mx_b_20170606_2146_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2203/calibration/mx_b_20170606_2203_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2221/calibration/mx_b_20170606_2221_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2240/calibration/mx_b_20170606_2240_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2256/calibration/mx_b_20170606_2256_000000.root",
                                      "HVCalibration/Background/processed/mx_b_20170606_2313/calibration/mx_b_20170606_2313_000000.root"};
 
    float hv[dim]       = {400,420,440,460,480,500,520,540,560,580,600};
    TString name_latex3 = "HV_{scionix} = 500 V";
    float minbin[dim] = {10,10,10,10,10,10,10,10,10,10,10};
    float maxbin[dim] = {30,30,40,50,40,50,70,100,150,200,300};

    const int dim_use = 7; 
    int useful[dim_use] = {4,5,6,7,8,9,10};
    float minfun[dim_use] = {13,16,24,32,42,60,82};
    float maxfun[dim_use] = {17,25,34,44,60,80,105};

    float parinit0[dim_use] = {300,500,500,500,500,500,500};
    float parinit1[dim_use] = {15,20,30,37,50,70,90};
    float parinit2[dim_use] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5};
    float parinit3[dim_use] = {100,100,100,100,100,100,100};
    float parinit4[dim_use] = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001};
    float parinit5[dim_use],parinit6[dim_use],parinit7[dim_use];
  
    */

//////// COMMON //////////////////////////////////////////
//////////////////////////////////////////////////////////

    TString name_spectrum             = "scan HV : detector ";
    TString name_spectrum_2           = "scan HV : detector ";
    TString name_spectrum_pdf         = "figuras/plots";
    TString name_spectrum_2_pdf       = "figuras/plots";
    TString name_spectrum_back        = "scan HV : background";
    TString name_HVscan_criteria1_pdf = "figuras/graph_criteria1_";
    TString name_HVscan_criteria1     = "scan HV : energy resolution criteria : detector ";
    TString name_HVscan_criteria2_pdf = "figuras/graph_criteria2_";
    TString name_HVscan_criteria2     = "scan HV : signal/background criteria : detector ";
    TString name_latex1               = "#frac{#Delta E}{E} using" + source_name;
    TString name_latex2               = "#frac{signal}{back} using" + source_name;
    
    name_detector = to_string(ichannel);
    name_spectrum += name_detector; 
    name_spectrum_2 += name_detector; 
    name_spectrum_pdf += name_detector;
    name_spectrum_pdf += ".pdf";
    name_spectrum_2_pdf += name_detector;
    name_spectrum_2_pdf += ".pdf";
    name_HVscan_criteria1_pdf += name_detector;
    name_HVscan_criteria1_pdf += ".pdf";
    name_HVscan_criteria2_pdf += name_detector;
    name_HVscan_criteria2_pdf += ".pdf";
    name_HVscan_criteria1 += name_detector;
    name_HVscan_criteria2 += name_detector;

    int channel,error;
    double time;
    float integral;

    float delta_sens[dim],significance[dim];
    float edelta_sens[dim];
    float esignificance[dim];
    float hv_use[dim_use],ehv_use[dim_use];
    float mean_sig,emean_sig,width_sig,ewidth_sig,yield_sig,eyield_sig;
    double yield_back,eyield_back;
    int nevents;

    TF1 *fun[dim_use];
    char namefun[20];
    double dx;

    TFile *file1,*file2;
    TTree *j1,*j2;
    TH1 *h[dim],*ht[dim_use],*hb[dim_use];
    char nameh[20],namel[20],namebh[20];
    TCanvas *l,*c,*d,*m;
    TLatex la0,la1,la2,la3,la4;

    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    int rest = dim%3;
    int num  = dim/3;
    if (rest!=0){
        num = num+1;
    } 
    rest = dim_use%3;
    int num_use  = dim_use/3;
    if (rest!=0){
        num_use = num_use+1;
    } 
    for (int k=0; k<dim_use; k++){
        edelta_sens[k] = 0.0;
        esignificance[k] = 0.0;
        ehv_use[k] = 0.0;
    }
        
    c = new TCanvas("c",name_spectrum,200,10,900,1300);
    c->Divide(num,3);

    for (int i=0; i<dim; i++){

        file1 = new TFile(namedata[i],"UPDATE");
        j1 = (TTree*)file1->Get("T"); 
        sprintf(nameh,"h%d",i); 
        sprintf(namel,"HV=%3.0f",hv[i]);
        h[i] = new TH1F(nameh, "integral", nbin, minbin[i], maxbin[i]);
        h[i]->SetMinimum(0);

        j1->SetBranchAddress("integral",&integral);
        j1->SetBranchAddress("channel",&channel); 
        j1->SetBranchAddress("time",&time); 
        j1->SetBranchAddress("error",&error); 

        nevents = j1->GetEntries();
        for (int k=0; k<nevents; k++){
            j1->GetEntry(k);
            if ((channel==ichannel)&&(error==0)){
                h[i]->Fill(integral*1000000000);
            }
        }
        c->cd(i+1);
        h[i]->SetMarkerColor(4);
        h[i]->SetMarkerStyle(21);
        h[i]->SetLineWidth(2);
        h[i]->SetLineColor(1);
        h[i]->Draw();
        h[i]->GetYaxis()->SetTitle("counts");
        h[i]->GetXaxis()->SetTitle("integral (V*s)");
        la0.SetNDC(kTRUE);
        la0.SetTextSize(0.08);
        la0.DrawLatex(0.65,0.8,namel);  
    }

    d = new TCanvas("d",name_spectrum_2,200,10,900,1300);
    d->Divide(num_use,3);

    for (int i=0; i<dim_use; i++){
        file1 = new TFile(namedata[useful[i]],"UPDATE");
        j1 = (TTree*)file1->Get("T"); 
        sprintf(nameh,"ht%d",i); 
        sprintf(namel,"HV=%3.0f",hv[useful[i]]);
        ht[i] = new TH1F(nameh, "integral",nbin_fun,minfun[i],maxfun[i]);
        ht[i]->SetMinimum(0);
        j1->SetBranchAddress("integral",&integral);
        j1->SetBranchAddress("channel",&channel); 
        j1->SetBranchAddress("time",&time); 
        j1->SetBranchAddress("error",&error); 
        nevents = j1->GetEntries();
        double maxTime = 0.0;
        for (int k=0; k<nevents; k++){
            j1->GetEntry(k);
            if ((channel==ichannel)&&(error==0)){
                if (maxTime<time) maxTime = time; 
                ht[i]->Fill(integral*1000000000);
            }
        }
        hv_use[i] = hv[useful[i]];
        dx = (maxfun[i]-minfun[i])/nbin_fun;
        sprintf(namefun,"fun%d",i); 
        if (model==1){
          fun[i] = new TF1(namefun,fun_twoSignal_Back,minfun[i],maxfun[i],8);
          fun[i]->SetParameter(0,parinit0[i]/dx);   
          fun[i]->SetParameter(1,parinit1[i]);
          fun[i]->SetParameter(2,parinit2[i]);
          fun[i]->SetParameter(3,parinit3[i]/dx);
          fun[i]->SetParameter(4,parinit4[i]);
          fun[i]->SetParameter(5,parinit5[i]);
          fun[i]->SetParameter(6,parinit6[i]);
          fun[i]->SetParameter(7,parinit7[i]);
          /*
          fun[i]->FixParameter(0,parinit0[i]/dx);   
          fun[i]->FixParameter(1,parinit1[i]);
          fun[i]->FixParameter(2,parinit2[i]);
          fun[i]->FixParameter(3,parinit3[i]/dx);
          fun[i]->FixParameter(4,parinit4[i]);
          fun[i]->FixParameter(5,parinit5[i]);
          fun[i]->FixParameter(6,parinit6[i]);
          fun[i]->FixParameter(7,parinit7[i]);
          */
        }
        if (model==0){
          fun[i] = new TF1(namefun,fun_Signal_Back,minfun[i],maxfun[i],5);
          fun[i]->SetParameter(0,parinit0[i]/dx);   
          fun[i]->SetParameter(1,parinit1[i]);
          fun[i]->SetParameter(2,parinit2[i]);
          fun[i]->SetParameter(3,parinit3[i]);
          fun[i]->SetParameter(4,parinit4[i]);
          //if (i==4) fun[i]->FixParameter(0,parinit0[i]/dx);
          //if (i==4) fun[i]->FixParameter(1,parinit1[i]);
          //if (i==4) fun[i]->FixParameter(2,parinit2[i]);
          //if (i==4) fun[i]->FixParameter(3,parinit3[i]);
          //if (i==4) fun[i]->FixParameter(4,parinit4[i]);           
        }     

        d->cd(i+1);
        ht[i]->SetMarkerColor(4);
        ht[i]->SetMarkerStyle(21);
        ht[i]->SetLineWidth(2);
        ht[i]->SetLineColor(1);

        ht[i]->Fit(fun[i],"R");
        ht[i]->GetYaxis()->SetTitle("counts");
        ht[i]->GetXaxis()->SetTitle("integral (V*s)");
        la0.SetNDC(kTRUE);
        la0.SetTextSize(0.08);
        la0.DrawLatex(0.65,0.8,namel);

        yield_sig  = fun[i]->GetParameter(0); 
        mean_sig   = fun[i]->GetParameter(1); 
        width_sig  = fun[i]->GetParameter(2); 
        eyield_sig = fun[i]->GetParError(0); 
        emean_sig  = fun[i]->GetParError(1); 
        ewidth_sig = fun[i]->GetParError(2);         
        delta_sens[i]  = 100*width_sig/mean_sig; 
        edelta_sens[i] = 100*ErrDiv(width_sig,ewidth_sig,mean_sig,emean_sig);

        file2 = new TFile(namedata_back[useful[i]],"UPDATE");
        j2 = (TTree*)file2->Get("T"); 
        sprintf(namebh,"hb%d",i); 
        int nbin_back = 6*width_sig/dx;
        float minback = mean_sig-3*width_sig;
        float maxback = mean_sig+3*width_sig;
        hb[i] = new TH1F(namebh, "integral",nbin_back,minback,maxback);
        hb[i]->SetMinimum(0);
        j2->SetBranchAddress("integral",&integral);
        j2->SetBranchAddress("channel",&channel); 
        j2->SetBranchAddress("time",&time); 
        j2->SetBranchAddress("error",&error); 

        nevents = j2->GetEntries();
        double maxTime_back = 0.0;
        for (int k=0; k<nevents; k++){
            j2->GetEntry(k);
            if ((channel==ichannel)&&(error==0)){
                if (maxTime_back<time) maxTime_back = time; 
                hb[i]->Fill(integral*1000000000);
            }
        } 
        hb[i]->Scale(maxTime/maxTime_back);
        yield_back = hb[i]->IntegralAndError(1,nbin_back,eyield_back);

        cout << " yields (sig & back) = " << yield_sig << " " << yield_back << endl;      

        //hb[i]->SetMarkerColor(1);
        //hb[i]->SetMarkerStyle(21);
        //hb[i]->SetLineWidth(0);
        //hb[i]->SetLineColor(1);
        //hb[i]->SetFillColor(4);
        hb[i]->Draw("SAME");

        significance[i]  = yield_sig/(dx*yield_back); 
        esignificance[i] = ErrDiv(yield_sig,eyield_sig,yield_back*dx,eyield_back*dx);  
    }

    d->SaveAs(name_spectrum_2_pdf);

    /// criteria 1
    l = new TCanvas("l",name_HVscan_criteria1,200,10,600,600);
    TGraphErrors *g1 = new TGraphErrors(dim_use,hv_use,delta_sens,ehv_use,edelta_sens);
    g1->SetTitle("");
    g1->SetMarkerColor(4);
    g1->SetMarkerStyle(21);
    g1->GetYaxis()->SetTitle("deltaE/E (%)");
    g1->GetXaxis()->SetTitle("HV");
    l->cd(1);
    g1->Draw("AP");
    //g1->Fit("expo");
    la1.SetNDC(kTRUE);
    la1.DrawLatex(0.3,0.8,name_latex1);
    la2.SetNDC(kTRUE);
    la2.DrawLatex(0.3,0.7,name_latex3);

    l->SaveAs(name_HVscan_criteria1_pdf);

    /// criteria 2
    m = new TCanvas("m",name_HVscan_criteria2,200,10,600,600);
    TGraphErrors *g2 = new TGraphErrors(dim_use,hv_use,significance,ehv_use,esignificance);
    g2->SetTitle("");
    g2->SetMarkerColor(4);
    g2->SetMarkerStyle(21);
    g2->GetYaxis()->SetTitle("signal/background");
    g2->GetXaxis()->SetTitle("HV");
    m->cd(1);
    g2->Draw("AP");
    //g2->Fit("expo");
    la3.SetNDC(kTRUE);
    la3.DrawLatex(0.3,0.8,name_latex2);
    la4.SetNDC(kTRUE);
    la4.DrawLatex(0.3,0.7,name_latex3);

    m->SaveAs(name_HVscan_criteria2_pdf);

}
