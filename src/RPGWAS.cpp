
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "RPGWAS.hpp"

namespace RPGWAS {

RPGWASClass::RPGWASClass(arma::vec t_Tarvec,                    
                         arma::vec t_Riskvec, 
                         arma::mat t_designMat,
                         Rcpp::DataFrame t_GRM,
                         arma::vec t_resid,
                         double t_lambda,
                         arma::vec t_gammas,
                         double t_gamma_riskVec,
                         arma::vec t_beta_null,
                         arma::vec t_resid_risk,
                         arma::vec t_t0,
                         arma::vec t_resid_unrelated_outliers,  
                         double t_sum_R_nonOutlier, 
                         double t_R_GRM_R_nonOutlier,    
                         double t_R_GRM_R_TwoSubjOutlier,   
                         double t_R_GRM_R,   
                         arma::vec t_MAF_interval, 
                         Rcpp::List t_TwoSubj_list,   
                         Rcpp::List t_ThreeSubj_list,   
                         double t_SPA_Cutoff,   
                         double t_zeta,            
                         double t_tol)
{ 
  m_Tarvec = t_Tarvec;                    
  m_Riskvec = t_Riskvec; 
  m_designMat = t_designMat;
  m_GRM = t_GRM;
  m_resid = t_resid;
  m_lambda = t_lambda;
  m_gammas = t_gammas;
  m_gamma_riskVec = t_gamma_riskVec;
  m_beta_null = t_beta_null;  
  m_resid_risk = t_resid_risk;
  m_t0 = t_t0;
  m_resid_unrelated_outliers = t_resid_unrelated_outliers;
  m_sum_unrelated_outliers2 = sum(t_resid_unrelated_outliers % t_resid_unrelated_outliers);
  m_sum_R_nonOutlier = t_sum_R_nonOutlier;
  m_R_GRM_R_nonOutlier = t_R_GRM_R_nonOutlier;
  m_R_GRM_R_TwoSubjOutlier = t_R_GRM_R_TwoSubjOutlier;
  m_R_GRM_R = t_R_GRM_R;
  m_MAF_interval = t_MAF_interval;
  m_TwoSubj_list = t_TwoSubj_list;
  m_ThreeSubj_list = t_ThreeSubj_list;
  m_SPA_Cutoff = t_SPA_Cutoff;
  m_zeta = t_zeta;
  m_tol = t_tol;
}
}