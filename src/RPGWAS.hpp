
#ifndef RPGWAS_HPP
#define RPGWAS_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <chrono>

#include "UTIL.hpp"

namespace RPGWAS{

class RPGWASClass
{
private:

  ////////////////////// -------------------- members ---------------------------------- //////////////////////
  
  arma::vec m_Tarvec;                    // Target phenotype
  arma::vec m_Riskvec;                   // Risk phenotype
  arma::mat m_designMat;                 // covariate matrix
  Rcpp::DataFrame m_GRM;                 // Genetic Relatedness Matrix
  arma::vec m_resid;                     // residuals
  double m_lambda;                       // precalculated lambda
  arma::vec m_gammas;                    // results of fit null logistic regression (tarVec ~ designMat + riskVec)
  double m_gamma_riskVec;                // results of fit null logistic regression (tarVec ~ designMat + riskVec)
  arma::mat m_inv_tX_X;                  // prepare for fit full linear regression (riskVec ~ designMat + geno)
  arma::mat m_inv_tX_X_tX;               // prepare for fit full linear regression (riskVec ~ designMat + geno)
  arma::vec m_t0;                        // fitted values of gammas %*% designMat
  arma::vec m_resid_unrelated_outliers;  // unrelated outlier residuals
  double m_sum_R_nonOutlier;             // sum of non-outlier residuals
  double m_sum_unrelated_outliers2;      // sum of squares of unrelated outlier residuals
  double m_R_GRM_R_nonOutlier;           // residuals x GRM x residuals for non-outlier families
  double m_R_GRM_R_TwoSubjOutlier;       // residuals x GRM x residuals for outlier families (n = 2)
  double m_R_GRM_R;                      // residuals x GRM x residuals
  arma::vec m_MAF_interval;              // MAF interval divides the MAFs into several intervals
  Rcpp::List m_TwoSubj_list;             // List of residuals and IBD probabilities in outlier families (n = 2)
  Rcpp::List m_ThreeSubj_list;           // List of residuals and Chow-Liu tree in outlier families (n > 2)
  
  double m_SPA_Cutoff;                   // cutoff of standardized score to use normal approximation or SPA
  double m_zeta;                         // initial saddle point for negative side, default is zero
  double m_tol;                          // accuracy of Newton's methods, default 1e-4 for beta; 1e-5 for tau 
  
public:
  
  RPGWASClass(arma::vec t_Tarvec,
              arma::vec t_Riskvec,
              arma::mat t_designMat,
              Rcpp::DataFrame t_GRM,
              arma::vec t_resid,
              double t_lambda,
              arma::vec t_gammas,
              double t_gamma_riskVec,
              arma::mat t_inv_tX_X,
              arma::mat t_inv_tX_X_tX,
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
              double t_tol);
  
  // function to fit the full linear regression (Riskvec ~ designMat + Gvec)
  Rcpp::List fit_full(arma::vec t_Gvec) 
  {
    
    arma::mat X = join_horiz(arma::ones<arma::vec>(m_Riskvec.n_rows), m_designMat); 
    arma::mat B = trans(X) * t_Gvec;
    double D = as_scalar(trans(t_Gvec) * t_Gvec);
    
    arma::mat inv_A_B = m_inv_tX_X_tX * t_Gvec;
    double M22_value = 1.0 / (D - as_scalar(trans(B) * inv_A_B));
    arma::mat M22 = arma::mat(1, 1, arma::fill::zeros);
    M22(0, 0) = M22_value;
    
    arma::mat M11 = m_inv_tX_X + inv_A_B * M22 * trans(inv_A_B);
    arma::mat M12 = - M11 * B / D;
    
    arma::mat M = arma::join_vert(
      arma::join_horiz(M11, M12),
      arma::join_horiz(trans(M12), M22)
    );
    
    arma::mat Z = join_horiz(X, t_Gvec);
    arma::vec est = M * trans(Z) * m_Riskvec;
    
    double sigma2_sq = as_scalar(trans(m_Riskvec - Z * est) * (m_Riskvec - Z * est)) / (m_Riskvec.n_rows - Z.n_cols);
    double var_geno = sigma2_sq * M(M.n_rows - 1, M.n_cols - 1);
    
    double score = est(est.n_elem - 1) * est(est.n_elem - 1) / var_geno;
    double pval = R::pchisq(score, 1, 0, 0);
    
    return Rcpp::List::create(
      Rcpp::Named("betas") = est.subvec(0, est.n_elem - 2),
      Rcpp::Named("beta.geno") = est(est.n_elem - 1), 
      Rcpp::Named("sigma2_sq") = sigma2_sq,
      Rcpp::Named("pval") = pval
    );
  }
  
  // function to generate new residuals using updated Y_hat
  arma::vec calculate_new_Residuals(arma::vec t_Gvec) 
  {
    arma::mat X = arma::join_horiz(arma::ones(m_Riskvec.n_rows), m_designMat);
    Rcpp::List fit_results = fit_full(t_Gvec);
    
    arma::vec new_Riskvec = m_Riskvec - Rcpp::as<double>(fit_results["beta.geno"]) * t_Gvec; // wait for checking
    
    arma::uvec obs_Y = arma::find_finite(m_Tarvec);
    arma::uvec miss_Y = arma::find_nonfinite(m_Tarvec);
    
    arma::vec t = m_t0 + m_gamma_riskVec * new_Riskvec.elem(obs_Y);
    arma::vec f = arma::normpdf(t);
    arma::vec F = arma::normcdf(t);
    arma::vec deno = F % (1 - F);
    
    arma::uvec indices1 = arma::find(F == 1);
    deno.elem(indices1) = arma::normcdf(-t.elem(indices1)) + 1e-200;
    
    arma::uvec indices2 = arma::find(F == 0);
    deno.elem(indices2) = deno.elem(indices2) + 1e-200;
    
    double sigma2_sq = Rcpp::as<double>(fit_results["sigma2_sq"]);
    arma::vec betas = arma::join_vert(Rcpp::as<arma::vec>(fit_results["betas"]), Rcpp::as<arma::vec>(fit_results["beta.geno"]));
    
    double ia2b2 = - m_gamma_riskVec * arma::sum(arma::pow(f, 2) % arma::pow(t_Gvec.elem(obs_Y), 2) / deno);
    double ib2b2 = arma::as_scalar(arma::sum(arma::pow(t_Gvec, 2)) / Rcpp::as<double>(fit_results["sigma2_sq"])) - (m_gamma_riskVec * ia2b2);
    double lambda = ia2b2 / ib2b2;
    
    arma::vec R_alpha2 = f % (m_Tarvec.elem(obs_Y) - F) / deno;
    arma::vec R1_beta2 = (1 / sigma2_sq) * (m_Riskvec - arma::join_horiz(X, t_Gvec) * betas);
    arma::vec R2_beta2 = - m_gamma_riskVec * R_alpha2;
    
    arma::vec resid_obs = R_alpha2 - lambda * (R1_beta2.elem(obs_Y) + R2_beta2);
    arma::vec resid_miss = -lambda * R1_beta2.elem(miss_Y);
    
    // Fill the Resid column in the DataFrame
    arma::vec resid(m_Riskvec.n_rows);
    for (unsigned int i = 0; i < obs_Y.n_elem; i++) {
      resid[obs_Y(i)] = resid_obs(i);
    }
    for (unsigned int i = 0; i < miss_Y.n_elem; i++) {
      resid[miss_Y(i)] = resid_miss(i);
    }
    
    // Return only the SubjID and Residual columns
    return resid;
  }
  
  // The MGF and its first and second derivative MGF of G (genotype)
  arma::mat MGF_cpp(double t, 
                    const Rcpp::List update_ThreeSubj_list,
                    double MAF)
  {
    // Unrelated subjects.
    arma::vec lambda = arma::exp(t * m_resid_unrelated_outliers);
    
    arma::vec alpha = 1 - MAF + MAF * lambda; 
    arma::vec alpha_1 = MAF * m_resid_unrelated_outliers % lambda; 
    arma::vec alpha_2 = m_resid_unrelated_outliers % alpha_1;
    
    arma::vec M_G0_all = alpha % alpha;
    arma::vec M_G1_all = 2 * alpha % alpha_1;
    arma::vec M_G2_all = 2 * (alpha_1 % alpha_1 + alpha % alpha_2);
    
    // Two related subjects in a family.
    int n1 = m_TwoSubj_list.length();
    if (n1 != 0)
    {
      for (int i = 0; i < n1; i++)
      {
        Rcpp::List TwoSubj_list_temp = m_TwoSubj_list[i];
        // arma::vec Resid = Rcpp::as<arma::vec>(TwoSubj_list_temp["Resid"]);
        // arma::vec Rho = Rcpp::as<arma::vec>(TwoSubj_list_temp["Rho"]);
        arma::vec Resid = TwoSubj_list_temp["Resid"];
        arma::vec Rho = TwoSubj_list_temp["Rho"];
        
        arma::vec temp = (1 - Rho) * MAF * (1 - MAF);
        
        double R1 = Resid[0]; double etR1 = exp(t * R1);
        double R2 = Resid[1]; double etR2 = exp(t * R2);
        double Rsum = R1 + R2;
        
        arma::vec midterm1 = etR1 * temp;
        arma::vec midterm2 = etR2 * temp;
        arma::vec midterm3 = etR1 * etR2 * (MAF - temp);
        
        arma::vec M_G0 = midterm1 + midterm2 + midterm3 - temp + 1 - MAF;
        arma::vec M_G1 = R1 * midterm1 + R2 * midterm2 + Rsum * midterm3;
        arma::vec M_G2 = R1*R1 * midterm1 + R2*R2 * midterm2 + Rsum*Rsum * midterm3;
        
        M_G0_all = arma::join_cols(M_G0_all, M_G0);
        M_G1_all = arma::join_cols(M_G1_all, M_G1);
        M_G2_all = arma::join_cols(M_G2_all, M_G2);
      }
    }
    
    // Three above Related Subjects.
    int n2 = update_ThreeSubj_list.length();
    if (n2 != 0)
    {
      for (int i = 0; i < n2; i++)
      {
        Rcpp::List ThreeSubj_list_temp = update_ThreeSubj_list[i];
        // arma::vec stand_S = Rcpp::as<arma::vec>(ThreeSubj_list_temp["stand.S"]);
        // arma::vec arr_prob = Rcpp::as<arma::vec>(ThreeSubj_list_temp["arr.prob"]);
        arma::vec stand_S = ThreeSubj_list_temp["stand.S"];
        arma::vec arr_prob = ThreeSubj_list_temp["arr.prob"];
        
        arma::vec midterm0 = exp(t * stand_S) % arr_prob;
        arma::vec midterm1 = stand_S % midterm0;
        arma::vec midterm2 = stand_S % midterm1;
        
        M_G0_all = arma::join_cols(M_G0_all, arma::vec{arma::accu(midterm0)});
        M_G1_all = arma::join_cols(M_G1_all, arma::vec{arma::accu(midterm1)});
        M_G2_all = arma::join_cols(M_G2_all, arma::vec{arma::accu(midterm2)});
      }
    }
    
    return arma::join_rows(M_G0_all, M_G1_all, M_G2_all);
  }
  
  // Newton's method to get the saddle point
  double fastgetroot_cpp(const Rcpp::List update_ThreeSubj_list,
                         double Score,
                         double MAF,
                         double init_t,
                         double tol,
                         int maxiter = 50)
  {
    double t = init_t;
    arma::vec MGF0; arma::vec MGF1; arma::vec MGF2;
    double CGF1 = 0; double CGF2 = 0;
    double diff_t = R_PosInf;
    int iter;
    
    double mean = 2 * MAF * m_sum_R_nonOutlier;
    double var = 2 * MAF * (1 - MAF) * m_R_GRM_R_nonOutlier;
    
    for (iter = 1; iter < maxiter; iter++)
    {
      double old_t = t;
      double old_diff_t = diff_t;
      double old_CGF1 = CGF1;
      
      arma::mat MGF_all = MGF_cpp(t, update_ThreeSubj_list, MAF);
      
      MGF0 = MGF_all.col(0);
      MGF1 = MGF_all.col(1);
      MGF2 = MGF_all.col(2);
      
      arma::vec temp = MGF1 / MGF0;
      CGF1 = arma::accu(temp) + mean + var * t - Score;
      CGF2 = arma::accu(MGF2 / MGF0) - arma::accu(temp % temp) + var;
      
      diff_t = - CGF1/CGF2;
      
      // std::cout << "iter:\t" << iter << std::endl;
      // std::cout << "t:\t" << t << std::endl;
      // std::cout << "CGF1:\t" << CGF1 << std::endl;
      // std::cout << "CGF2:\t" << CGF2 << std::endl;
      // std::cout << "diff_t:\t" << diff_t << std::endl;
      // std::cout << std::endl;
      
      if (std::isnan(diff_t) || std::isinf(CGF2))
      {
        t = t / 2;
        diff_t = std::min(std::abs(t), 1.0) * arma::sign(Score);
        continue;
      }
      
      if (std::isnan(old_CGF1) || (arma::sign(old_CGF1) != 0 && arma::sign(CGF1) != arma::sign(old_CGF1))) 
      {
        if (std::abs(diff_t) < tol) 
        {
          t = old_t + diff_t;
          break;
        } else {
          while (std::abs(old_diff_t) > tol && std::abs(diff_t) > std::abs(old_diff_t) - tol) 
          {
            diff_t = diff_t / 2;
          }
          t = old_t + diff_t;
          continue;
        }
      }
      
      if (arma::sign(Score) != arma::sign(old_t + diff_t) && 
          (arma::sign(old_CGF1) == 0 || arma::sign(CGF1) == arma::sign(old_CGF1))) 
      {
        while (arma::sign(Score) != arma::sign(old_t + diff_t)) 
        {
          diff_t = diff_t / 2;
        }
        t = old_t + diff_t;
        continue;
      }
      
      t = old_t + diff_t;
      if (std::abs(diff_t) < tol) break;
    }
    
    // return Rcpp::List::create(Named("root") = t,
    //                           Named("iter") = iter);
    
    return t;
  }
  
  // function to get one side p value
  double GetProb_SPA(const Rcpp::List update_ThreeSubj_list,
                     double Score,
                     double MAF,
                     bool lower_tail,
                     double zeta,
                     double tol)
  {
    zeta = fastgetroot_cpp(update_ThreeSubj_list, Score, MAF, zeta, tol);
    
    arma::mat MGF_all = MGF_cpp(zeta, update_ThreeSubj_list, MAF);
    
    arma::vec MGF0 = MGF_all.col(0);
    arma::vec MGF1 = MGF_all.col(1);
    arma::vec MGF2 = MGF_all.col(2);
    
    double mean = 2 * MAF * m_sum_R_nonOutlier;
    double var = 2 * MAF * (1 - MAF) * m_R_GRM_R_nonOutlier;
    
    arma::vec temp = MGF1 / MGF0;
    double CGF0 = arma::accu(log(MGF0)) + mean * zeta + 0.5 * var * zeta * zeta;
    double CGF2 = arma::accu(MGF2 / MGF0) - arma::accu(temp % temp) + var;
    
    double w = arma::sign(zeta) * sqrt(2 * (zeta * Score - CGF0));
    double v = zeta * sqrt(CGF2);
    
    double u = w + 1/w * log(v/w);
    double pval = R::pnorm(u, 0, 1, lower_tail, false);
    
    // std::cout << "zeta:\t" << zeta << std::endl;
    // std::cout << "p value:\t" << pval << std::endl;
    // std::cout << std::endl;
    
    return pval;
  }
  
  // function to get two side p values
  double getMarkerPval(arma::vec t_Gvec,
                       double t_altFreq,
                       double& t_zScore,
                       double& t_hwepval,
                       double t_hwepvalCutoff)
  {
    // updated on 2023-05-23 to get hwe pvalue
    gethwepval(t_Gvec, t_hwepval, t_hwepvalCutoff);
    
    double MAF = std::min(t_altFreq, 1 - t_altFreq);
    
    // std::cout << arma::mean(t_GVec) << std::endl;
    // std::cout << t_altFreq << std::endl;
    
    double G_var = 2 * MAF * (1 - MAF);
    double Score, Score_var;
    
    Rcpp::List fit_results = fit_full(t_Gvec);
    
    if (Rcpp::as<double>(fit_results["pval"]) > 1e-3)
    {
      Score = sum(t_Gvec % m_resid) - mean(t_Gvec) * sum(m_resid);
      Score_var = G_var * m_R_GRM_R;
      t_zScore = Score/sqrt(Score_var);
      
    } 
    else{
      arma::vec new_resid = calculate_new_Residuals(t_Gvec);
      Score = sum(t_Gvec % new_resid) - mean(t_Gvec) * sum(new_resid);

      Rcpp::IntegerVector indice1 = m_GRM["indice1"];
      Rcpp::IntegerVector indice2 = m_GRM["indice2"];
      arma::uvec indices1 = Rcpp::as<arma::uvec>(indice1) - 1;  // R to C++ index (0-based)
      arma::uvec indices2 = Rcpp::as<arma::uvec>(indice2) - 1;  // R to C++ index (0-based)

      arma::vec pos1 = new_resid.elem(indices1);
      arma::vec pos2 = new_resid.elem(indices2);
      arma::vec Values = Rcpp::as<arma::vec>(m_GRM["Value"]);
      arma::vec Cov = arma::abs(Values) % pos1 % pos2;

      double R_GRM_R = sum(Cov);
      Score_var = G_var * R_GRM_R;
      t_zScore = Score/sqrt(Score_var);
    }

    if (std::abs(t_zScore) <= m_SPA_Cutoff)
    {
      // std::cout << "t_zScore:" << t_zScore << std::endl;
      double pval = R::pnorm(std::abs(t_zScore), 0, 1, false, false);
      
      return 2 * pval;
    }
    
    int order2 = arma::index_max(m_MAF_interval >= MAF);
    int order1 = order2 - 1;
    
    // std::cout << order2 << " " << order1 << "\t";
    
    // if (MAF <= m_MAF_interval[0] || MAF > 0.5)
    // {
    //   Rcpp::stop("Minor allele frequency is out of MAF_interval, MAF is\t", MAF);
    // }
    
    double MAF_ratio = (m_MAF_interval[order2] - MAF)/(m_MAF_interval[order2] - m_MAF_interval[order1]);
    
    double Var_ThreeOutlier = 0;
    
    int n1 = m_ThreeSubj_list.length();
    Rcpp::List update_ThreeSubj_list(n1);
    
    if (n1 != 0)
    {
      for (int i = 0; i < n1; i++)
      {
        Rcpp::List ThreeSubj_list_temp = m_ThreeSubj_list[i];
        
        // arma::vec CLT_temp1(243, arma::fill::zeros);
        // arma::vec CLT_temp2(243, arma::fill::zeros);
        // arma::vec stand_S(243, arma::fill::zeros);
        
        arma::mat CLT_temp =  ThreeSubj_list_temp["CLT"];
        arma::vec stand_S = ThreeSubj_list_temp["stand.S"];
        arma::vec CLT_temp1 = CLT_temp.col(order1);
        arma::vec CLT_temp2 = CLT_temp.col(order2);
        
        arma::vec arr_prob = MAF_ratio * CLT_temp1 + (1 - MAF_ratio) * CLT_temp2;
        
        update_ThreeSubj_list[i] = Rcpp::List::create(Rcpp::Named("stand.S") = stand_S,
                                                      Rcpp::Named("arr.prob") = arr_prob);
        
        arma::vec temp1 = stand_S % arr_prob;
        
        double temp2 = arma::accu(temp1);
        double Var_ThreeOutlier_temp = arma::accu(stand_S % temp1) - temp2 * temp2;
        Var_ThreeOutlier = Var_ThreeOutlier + Var_ThreeOutlier_temp;
      }
    }
    
    double Var_nonOutlier = G_var * m_R_GRM_R_nonOutlier;

    double Var_unrelated_outliers = G_var * m_sum_unrelated_outliers2;
    double Var_TwoOutlier = G_var * m_R_GRM_R_TwoSubjOutlier;
    
    double EmpVar = Var_nonOutlier + Var_unrelated_outliers + Var_TwoOutlier + Var_ThreeOutlier;
    
    double Var_Ratio = Score_var / EmpVar;
    double Score_adj = Score / sqrt(Var_Ratio);
    
    double zeta1 = std::abs(Score_adj) / Score_var; zeta1 = std::min(zeta1, 1.2);
    double zeta2 = - std::abs(m_zeta);
    
    double pval1 = GetProb_SPA(update_ThreeSubj_list, std::abs(Score_adj), MAF, false, zeta1, 1e-4);
    double pval2 = GetProb_SPA(update_ThreeSubj_list, -std::abs(Score_adj), MAF, true, zeta2, m_tol);
    double pval = pval1 + pval2;
    
    // std::cout << "pval:" << pval << "\t";
    return pval;
  }
  
};

}

#endif
