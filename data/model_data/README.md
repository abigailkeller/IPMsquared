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
at site *i* in time series data

*f_index*: binary indicator of Fukui traps at time *t*, in trap *j*, at
site *i* in time series data

*m_index*: binary indicator of Minnow traps at time *t*, in trap *j*, at
site *i* in time series data

*s_index*: binary indicator of Shrimp traps at time *t*, in trap *j*, at
site *i* in time series data

*recruit_intro*: binary data structure indicating the time point when
recruits enter the process model at site *i*, corresponding to t = 6

*index_frac*: date index (i.e., D(t)) for each time *t* and site *i*,
expressed as a fraction of a year (i.e., 0.5 = beginning of July)

*n_cap*: total count of crabs removed, $n_{t,i,y}^R$, at time *t*, site
*i*, size *y* in time series data

*counts*: count of crabs removed, $n_{t,j,i,y}^C$, at time *t*, trap
*j*, site *i*, size *y* in time series data


**Size-at-age data (D2)**

*growth_data*: size (CW), age, and year observed of each crab 


**Batch mark-recapture data (D3)**

*f_index_mc*: binary indicator of Fukui traps at obs *s* in
mark-recapture data

*m_index_mc*: binary indicator of Minnow traps at obs *s* in
mark-recapture data

*s_index_mc*: binary indicator of Shrimp traps at obs *s* in
mark-recapture data

*m2_mc*: count of recaptured marked crabs of size *y* in mark-recapture
data

*n1_mc*: count of marked crabs of size *y* in mark-recapture data

*soak_days_mc*: number of soak days for each trap *s* in mark-recapture
data

*mc_index*: date index (i.e., D(t)) for the mark time point and
recapture time point, expressed as a fraction of a year (i.e., 0.5 =
beginning of July)

*totalo_mc*: total number of observations (traps) in mark-recapture data
