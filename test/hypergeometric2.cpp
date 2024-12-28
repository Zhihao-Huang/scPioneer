#include <Rcpp.h>
#include <regex>
using namespace Rcpp;
using namespace std;

double ramanujan(double x){
  return x*log(x) - x + log(x*(1 + 4*x*(1+2*x)))/6 + log(3.141593)/2;
}
// [[Rcpp::export]]
double bignchoosek(double n, double k){
  if (n == k) {
    return 1;
  }else if(k == 0){
    return 0;
  }else{
    return exp(ramanujan(n) - ramanujan(k) - ramanujan(n-k));
  }
}
// Calculate P value of hypergeometric test.
// [[Rcpp::plugins(cpp17)]]
//' @useDynLib angrycell, .registration = TRUE
//' @importFrom Rcpp evalCpp
//' @export
// [[Rcpp::export]]
NumericVector get_Terms_P_cpp2(CharacterVector candidate, double N, CharacterVector terms) {
  double n = candidate.size();
  int T = terms.size();
  NumericVector pvalue(T);
  for (int i=0; i<T; ++i){
    //split data
    auto s = as<string>(terms[i]);
    const char * oneterm = s.c_str();
    Rcpp::CharacterVector res;
    std::regex re(",");
    auto it = std::cregex_token_iterator(oneterm, oneterm + std::strlen(oneterm), re, -1);
    auto end = std::cregex_token_iterator();
    while(it != end) {
      res.push_back(*it);
      ++it;
    }
    double M = res.size();
    //Calculate share number k.
    int k = 0;
    for (CharacterVector::iterator it = candidate.begin(); it != candidate.end(); ++it) {
      if (std::find(res.begin(), res.end(), *it) != res.end()) {
        k += 1;
      }
    }
    double pv = 0;
    //calculate number of set; Similar to R function choose()
    double y1 = M;
    double z1 = k;
    double p1 = bignchoosek(y1,z1);
    double y2 = N-M;
    double z2 = n-k;
    double p2 = bignchoosek(y2,z2);
    double y3 = N;
    double z3 = n;
    double p3 = bignchoosek(y3,z3);
    pv = p1 * p2 / p3;
    pvalue[i] = pv;
    Rcout << "The value of p1 : " << p1 << "\n";
    Rcout << "The value of p2 : " << p2 << "\n";
    Rcout << "The value of p3 : " << p3 << "\n";
  }
  return pvalue;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
get_Terms_P_cpp2(candidate1, 13000, db$geneSymbol)
bignchoosek(200,50)
bignchoosek(200,200)
bignchoosek(200,0)
*/



