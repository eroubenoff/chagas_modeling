simul <- "
// compute enso
double ensoSig = A*((tan(((h*3.14159)/2)*(enso_jan/3)))/tan((h*3.14159)/2));

// compute transmission rate
double beta_I = exp(season1*b1_I+season2*b2_I+season3*b3_I+season4*b4_I+season5*b5_I+season6*b6_I+(season4*b4E_I + season5*b5E_I)*ensoSig);

//gamma white noise
double dW_I = rgammawn(sigPRO,dt);

// compute transition numbers
double dBS_I = (delta*popI+dpopdtI)*dt;
double dSI_I = exp(-nu*(L - 1995))*((I_I/popI)*beta_I + omega)*S_I*dW_I;
double dIR_I = muIR*I_I*dt;
double dRS_I = muRS*R_I*dt;
double dSD_I = delta*S_I*dt;
double dID_I = delta*I_I*dt;      
double dRD_I = delta*R_I*dt;  

// compute equations
S_I += dBS_I - dSI_I - dSD_I + dRS_I;
I_I += dSI_I - dIR_I - dID_I;
R_I += dIR_I - dRS_I - dRD_I;
if(t < 1996) {cases_I += rho*dSI_I;} 
else {cases_I += (rho/2)*dSI_I;}
W_I += (dW_I-dt)/sigPRO;

// track errors
if (S_I < 0.0) { err -= S_I; S_I=0.0; }
if (I_I < 0.0) { err -= I_I; I_I=0.0; }
if (R_I < 0.0) { err -= R_I; R_I=0.0; }
"
