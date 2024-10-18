#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  // Input Data
  
  DATA_INTEGER(NSA1); 

  DATA_STRUCT(spde,spde_t); // INLA SPDE object (components of precision matrix)
  DATA_SPARSE_MATRIX(A); // INLA SPDE projection matrix: mesh to pixels [dim: NSA1 x nMesh]

  DATA_ARRAY(npos); // [dim: NSA1 x 3]
  DATA_ARRAY(nneg); // [dim: NSA1 x 3]

  DATA_IVECTOR(sa1_ses); // [dim: NSA1]
  DATA_IVECTOR(sa1_ra); // [dim: NSA1]
  
  // Parameters
  
  PARAMETER(log_range_young);
  PARAMETER(log_sd_young);
  PARAMETER(log_range_old);
  PARAMETER(log_sd_old);
  
  PARAMETER(log_sd_ses_young);
  PARAMETER(log_sd_ses_old);
  PARAMETER(log_sd_ra_young);
  PARAMETER(log_sd_ra_old);

  PARAMETER(intercept_young);
  PARAMETER(intercept_old);
  PARAMETER(intercept_middle);
  
  PARAMETER(logit_middle_age_prop);
  
  PARAMETER_VECTOR(ses_effects_young); // [dim: 10]
  PARAMETER_VECTOR(ra_effects_young); // [dim: 5]
  PARAMETER_VECTOR(ses_effects_old); // [dim: 10]
  PARAMETER_VECTOR(ra_effects_old); // [dim: 5]
  
  PARAMETER_VECTOR(field_young); // [dim: nMesh]
  PARAMETER_VECTOR(field_old); // [dim: nMesh]
  
  // Parameter Transforms
  
  Type range_young = exp(log_range_young);
  Type kappa_young = 2.8284/range_young;
  Type sd_young = exp(log_sd_young);
  Type range_old = exp(log_range_old);
  Type kappa_old = 2.8284/range_old;
  Type sd_old = exp(log_sd_old);
  
  Type sd_ses_young = exp(log_sd_ses_young);
  Type sd_ses_old = exp(log_sd_ses_old);
  Type sd_ra_young = exp(log_sd_ra_young);
  Type sd_ra_old = exp(log_sd_ra_old);

  Type middle_age_prop = Type(1.0)/(Type(1.0)+exp(-logit_middle_age_prop));
  
  // Priors
  
  Type nll = 0.0;
  
  nll -= dnorm(log_sd_young,Type(-1),Type(0.5),true);
  nll -= dnorm(log_range_young,Type(1),Type(0.5),true);
  SparseMatrix<Type> Q_young = Q_spde(spde,kappa_young);
  nll += SCALE(GMRF(Q_young),sd_young)(field_young);

  nll -= dnorm(log_sd_old,Type(-1),Type(0.5),true);
  nll -= dnorm(log_range_old,Type(1),Type(0.5),true);
  SparseMatrix<Type> Q_old = Q_spde(spde,kappa_old);
  nll += SCALE(GMRF(Q_old),sd_old)(field_old);

  nll -= dnorm(log_sd_ses_young,Type(-1),Type(0.5),true);
  nll += SCALE(AR1(Type(0.9)),sd_ses_young)(ses_effects_young);
  nll -= dnorm(log_sd_ses_old,Type(-1),Type(0.5),true);
  nll += SCALE(AR1(Type(0.9)),sd_ses_old)(ses_effects_old);

  nll -= dnorm(log_sd_ra_young,Type(-1),Type(0.5),true);
  nll += SCALE(AR1(Type(0.9)),sd_ra_young)(ra_effects_young);
  nll -= dnorm(log_sd_ra_old,Type(-1),Type(0.5),true);
  nll += SCALE(AR1(Type(0.9)),sd_ra_old)(ra_effects_old);

  nll -= dnorm(logit_middle_age_prop,Type(0),Type(1),true);

  // Algebra

  vector<Type> baseline_field_young(NSA1);
  vector<Type> full_ses_effects_young(NSA1);
  for (int i=0; i<NSA1; i++) {full_ses_effects_young[i] = ses_effects_young[sa1_ses[i]];}
  vector<Type> full_ra_effects_young(NSA1);
  for (int i=0; i<NSA1; i++) {full_ra_effects_young[i] = ra_effects_young[sa1_ra[i]];}
  baseline_field_young = A*field_young;
  baseline_field_young = baseline_field_young.array() + intercept_young + full_ses_effects_young.array() + full_ra_effects_young.array();

  vector<Type> baseline_field_old(NSA1);
  vector<Type> full_ses_effects_old(NSA1);
  for (int i=0; i<NSA1; i++) {full_ses_effects_old[i] = ses_effects_old[sa1_ses[i]];}
  vector<Type> full_ra_effects_old(NSA1);
  for (int i=0; i<NSA1; i++) {full_ra_effects_old[i] = ra_effects_old[sa1_ra[i]];}
  baseline_field_old = A*field_old;
  baseline_field_old = baseline_field_old.array() + intercept_old + full_ses_effects_old.array() + full_ra_effects_old.array();

  matrix<Type> joint_field(NSA1,3);
  joint_field.col(0) = baseline_field_young.array();
  joint_field.col(2) = baseline_field_old.array();
  joint_field.col(1) = baseline_field_young.array()*middle_age_prop + baseline_field_old.array()*(Type(1.0)-middle_age_prop) + intercept_middle;

  // Likelihood

  vector<Type> npos_alt(NSA1);
  vector<Type> ntot_alt(NSA1);
  vector<Type> field_alt(NSA1);
  for (int i=0; i<3; i++) {
    npos_alt = npos.col(i);
    ntot_alt = npos.col(i)+nneg.col(i);
    field_alt = joint_field.col(i);
    nll -= dbinom_robust(npos_alt,ntot_alt,field_alt,true).sum();
  }

  // Reporting

  REPORT(baseline_field_young);
  REPORT(baseline_field_old);
  REPORT(joint_field);

  return nll;
}
