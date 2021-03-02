simul <- "
	// compute transmission rate
double beta_I = exp(season1*b1_I+season2*b2_I+season3*b3_I+season4*b4_I+season5*b5_I+season6*b6_I+season4*bf_I*flooding);
double lambda_I = beta_I*((I1_I + I2_I + I3_I)/popI);

//gamma white noise
double dW_I = rgammawn(sigPRO_I,dt);

double sprob = exp(-delta*dt); // survival probability
double rprob = 1-exp(-gamma*dt); // prob. of recovery
double rprob3 = 1-exp(-(gamma*eta)*dt); // prob. of recovery

// INNER
double dBM_I = ((delta*popI)+dpopdtI)*dt;
double dMS1_I = M_I*(1-exp(-omeMS*dt));
double dS1I1_I = S1_I*(1-exp(-lambda_I*dW_I));
double dS2I2_I = S2_I*(1-exp(-lambda_I*sigma*dW_I));
double dS3I3_I = S3_I*(1-exp(-lambda_I*sigma*sigma*dW_I));
double dI1S2_I = I1_I*rprob;
double dI2S3_I = I2_I*rprob;
double dI3S3_I = I3_I*rprob3;

M_I  += dBM_I-dMS1_I;
S1_I += dMS1_I-dS1I1_I;
I1_I += dS1I1_I-dI1S2_I;
S2_I += dI1S2_I-dS2I2_I;
I2_I += dS2I2_I-dI2S3_I;
S3_I += dI2S3_I-dS3I3_I+dI3S3_I;
I3_I += dS3I3_I-dI3S3_I;

M_I  *= sprob;
S1_I *= sprob;
I1_I *= sprob;
S2_I *= sprob;
I2_I *= sprob;
S3_I *= sprob;
I3_I *= sprob;

cases_I += (rho/2)*(dS1I1_I+dS2I2_I);
W_I += (dW_I-dt)/sigPRO_I;
"
