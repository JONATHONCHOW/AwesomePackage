#include <Rcpp.h>
using namespace Rcpp;

// Hello world by Rcpp.
// [[Rcpp::export]]
String rcpp_hello_world()
{
  String s("Data science is fantastic!");
  return(s);
}
