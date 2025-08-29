**Description of model data files:**

**Integral projection model data**

*b*: boundary points of discretized crab size in integral projection
model

*y*: midpoint of discretized crab size of integral projection model

**Time series data (D1)**

*totalt*: total number of time points in each year *i* in time series
data

*totalo*: total number of observations (traps) in year *i*, at time *t*
in time series data

*soak_days*: number of soak days for each trap in time *t*, at trap *j*,
at year *y* in time series data

*f_index*: binary indicator of Fukui traps at time *t*, in trap *j*, at
year *y* in time series data

*m_index*: binary indicator of Minnow traps at time *t*, in trap *j*, at
year *y* in time series data

*s_index*: binary indicator of Shrimp traps at time *t*, in trap *j*, at
year *y* in time series data

*recruit_intro*: binary data structure indicating the time point when
recruits enter the process model at year *y*, corresponding to t = 6

*index_frac*: date index (i.e., D(t)) for each time *t* and year *y*,
expressed as a fraction of a year (i.e., 0 =
beginning of April, 0.5 = beginning of October)

*n_cap*: total count of crabs removed, $n_{t,i,y}^R$, at time *t*, site
*i*, size *x* in time series data

*counts*: count of crabs removed, $n_{t,j,i,y}^C$, at time *t*, trap
*j*, year *y*, size *x* in time series data


**Size-at-age data (D2)**

*growth_data*: size (CW), age, and year observed of each crab 


**Mark-recapture data (D3)**

*roche_mc_mark*: count of crabs marked at time *t* of size *x* in mark-recapture data

*roche_mc_catch*: count of crabs recaptured at time *t* of size *x* in mark-recapture data

*f_index_mc_rc*: binary indicator of Fukui traps at obs *j*, time *t* in
mark-recapture data

*m_index_mc_rc*: binary indicator of Minnow traps at obs *j*, time *t* in
mark-recapture data

*s_index_mc_rc*: binary indicator of Shrimp traps at obs *j*, time *t* in
mark-recapture data

*soak_days_mc_rc*: number of soak days for each trap *j*, time *t* in mark-recapture
data

*mc_index_rc*: date index (i.e., D(t)) for the mark time point and
recapture time point, expressed as a fraction of a year (i.e., 0 =
beginning of April, 0.5 = beginning of October)

*roche_mc_totalo*: total number of observations (traps) at time *t* in mark-recapture data


**Alternative mark-recapture data (D3)**

*m2_mc_sl*: count of recaptured marked crabs of size *x* in alternative mark-recapture
data

*n1_mc_sl*: count of marked crabs of size *x* in alternative mark-recapture data

*mc_index_sl*: date index (i.e., D(t)) for the mark time point and
recapture time point, expressed as a fraction of a year (i.e., 0 =
beginning of April, 0.5 = beginning of October) in alternative mark-recapture dataset

*totalo_mc_sl*: total number of traps set in recapture period in alternative mark-recapture dataset

*f_index_mc_sl*: binary indicator of Fukui traps at obs *j* during recapture period in alternative
mark-recapture dataset

*m_index_mc_sl*: binary indicator of Minnow traps at obs *j* during recapture period in alternative
mark-recapture dataset

*s_index_mc_sl*: binary indicator of Shrimp traps at obs *j* during recapture period in alternative
mark-recapture dataset

*soak_days_mc_sl*: number of soak days for each trap *j* during recapture period in alternative
mark-recapture dataset

