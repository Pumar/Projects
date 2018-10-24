#include <iostream>
#include <cmath>
using namespace std;

void h12ascii (TH1* h)
{
   Int_t n = h->GetNbinsX();
   
   for (Int_t i=1; i<=n; i++) {
      printf("%g %g\n",
             h->GetBinLowEdge(i)+h->GetBinWidth(i)/2,
             h->GetBinContent(i));
   }
}
