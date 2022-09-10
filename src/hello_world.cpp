#include <Rcpp.h>
using namespace Rcpp;

// Hello world by Rcpp.
// [[Rcpp::export]]
Rcpp::String rcpp_hello_world()
{
  Rcpp::String s("Data science is fantastic!");
  return(s);
}
