#include <Rcpp.h>
#include <map>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector idUtility_rcpp(double bID, double dID, double aID, double n, 
                             const LogicalMatrix& ok, NumericVector group,
                             Nullable<NumericMatrix> ID_ = R_NilValue) {
  
  // allocate output
  const int xrows = ok.nrow();
  const int xcols = ok.ncol();
  NumericVector out(xrows*xcols);

  
  Rcpp::Function rownames("rownames");
  
  // None in front return -Inf for cells blocked by objects
  if ((Rf_isNull(ID_)) | (is_false(any(ok)))){
    
    for(int row = 0; row < xrows; row++) {
      for(int col = 0; col < xcols; col++) {
        
        if (ok(row,col)){
          out(row*xcols+col) = 0;
        } else{
          out(row*xcols+col) = R_NegInf;
        }
      }
    }
    return out;
  
  } else { 
    
    NumericMatrix ID(ID_); //re-initialise ID because it is not NULL
    
    Rcpp::Function rownames("rownames");
    Rcpp::CharacterVector ID_row_names = rownames(ID);
    
    CharacterVector group_names = group.names();
    group_names.erase(0); // remove 1st element
    
    const int ID_xrows = ID.nrow();
    NumericVector bID_2 = (ID_xrows);
    
    // Group dependent b, bigger for outgroup by dID
    for(int row = 0; row < ID_xrows; row++) {
      if(ID_row_names[row]==group_names(row)){
        bID_2(row)= bID;
        
      } else {
        
        bID_2(row)= bID+dID;
      }
    }
    
    // Object or pedestrian clash
    const int ok_xcols = ok.ncol();
    const int ok_xrows = ok.nrow();
    NumericMatrix out(ok);
    bool boolean_value;
    
    for (int i = 0; i < ok_xrows; i++) {  
      for (int j = 0; j < ok_xcols; j++) {
        boolean_value = ok(i,j) & (ID(1,j)>0);
        if(boolean_value){
          out(i,j) = 0;
        } else{
          out(i,j) = R_NegInf;
        }
      }
    }
    
    // Repulsion
    if(bID !=0){
      // out[ok] <- -apply(bID / (ID[, ok, drop = FALSE]^p["aID"]), 2, sum) // more info about ID
    }
    return out;
  }
  
}
