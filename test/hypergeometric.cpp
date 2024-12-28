#include <Rcpp.h>
#include <regex>
using namespace Rcpp;
using namespace std;


// Calculate P value of hypergeometric test.
// [[Rcpp::plugins(cpp17)]]
//' @useDynLib angrycell, .registration = TRUE
//' @importFrom Rcpp evalCpp
//' @export
// [[Rcpp::export]]
NumericVector get_Terms_P_cpp(CharacterVector candidate, double N, CharacterVector terms) {
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
    double a = 1;
    double b = 1;
    for (double z=1; z<k+1; ++z) {
      a *= z;
      b *= M-z+1;
    }
    double p1 = b/a;
    double a2 = 1;
    double b2 = 1;
    for (double z=1; z<n-k+1; ++z) {
      a2 *= z;
      b2 *= N-M-z+1;
    }
    double p2 = b2/a2;
    double a3 = 1;
    double b3 = 1;
    for (double z=1; z<n+1; ++z) {
      a3 *= z;
      b3 *= N-z+1;
    }
    double p3 = b3/a3;
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
get_Terms_P_cpp(candidate, 13000, db$geneSymbol)
*/

