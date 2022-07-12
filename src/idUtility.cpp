#include <Rcpp.h>
#include <map>

using namespace Rcpp;

//' Compute IDUtility
//'
//' Inter-personal distance utility for cell 1..33. b parameter divided by
//' sum over power of distances between bodies for cell to all inFront peds.
//' 
//' @param bID Numeric scalar.
//' @param dID Numeric scalar.
//' @param aID Numeric scalar.
//' @param n Numeric scalar indexing the subject in the state.
//' @param ok Logical matrix indicating if cells are blocked.
//' @param group Named numeric scalar with group indices for each pedestrian.
//' @param ID_ Numeric matrix of the type NULL - if not NULL, ID_ is a Numeric matrix of predicted distances from the subject to other pedestrians in the front.
//' @returns Numeric vector with interpersonal distance utilities for each cell.
//' @rdname ID
//' @aliases ID
//' @export
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
    
    NumericMatrix ID(ID_); //re-initialise ID if it is not NULL
    
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
    
    // Repulsion
    NumericVector out_as_vec(33);
    if(bID !=0){
      const int ID_xcols = ID.ncol();
      NumericMatrix out_temp(ID);
      
      for (int col = 0; col < ID_xcols; col++){
        for(int row = 0; row < ID_xrows; row++){
          out_temp(row,col) = -(bID_2[row] / (std::pow(ID(row, col), aID)));
        }
        out_as_vec[col] = sum(out_temp(_, col));
      }
      return out_as_vec;
      
    } else {
      
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
          } else {
            out(i,j) = R_NegInf;
          }
        }
      }
      
      // Return a flattened array of the out matrix (i.e., as_vector(out))
      for (int i = 0; i < ok_xrows; i++) { 
        for (int j = 0; j < ok_xcols; j++) {
          out_as_vec(i*ok_xrows+j) = out(i,j);
        }
      }
      return out_as_vec;
    } 
  }
}
