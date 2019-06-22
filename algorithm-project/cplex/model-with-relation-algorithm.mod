/*********************************************
 * OPL 12.9.0.0 Model
 * Author: aweftr
 * Creation Date: Jun 22, 2019 at 3:20:31 PM
 *********************************************/
int I0 = ...;
int J0 = ...;
int k = ...;
float c = ...;
//float b = ...;
int Ij[1..J0] = ...;
//int mySeed;
float time_m[1..(k-1)] = ...;
float time_c[k..I0][1..J0] = ...;
float time_exp[1..I0] = ...;
/*
execute{
var now = new Date();
mySeed = Opl.srand(Math.round(now.getTime()/1000));
}

int time_m[1] = 8 + rand(tabSize);
*/

range process_range = 1..I0;
range customer_range = 1..J0;

float C_c[process_range][customer_range] = ...;
float C_m[process_range] = ...;

float C_c_ext[process_range][customer_range] = ...;
float C_m_ext[process_range] = ...;

float P_c[process_range][customer_range] = ...;
float P_m[process_range] = ...;

float service_time_dis[process_range][1..2] = ...;

float U_c[process_range][customer_range] = ...;
float U_m[process_range] = ...;

float T_c1[process_range][customer_range]= ...;
float T_m1[process_range]= ...;
float T_c2[process_range][customer_range]= ...;
float T_m2[process_range]= ...;

float Tj_exp[customer_range] = ...;

dvar float T_c_ext_dvar[process_range][customer_range];
dvar float T_m_ext_dvar[process_range];

minimize
  (1 - (((sum (i in 1..(k-1))
    ((1 - abs((time_exp[i] - time_m[i]) / time_exp[i])) * ((time_m[i] * C_m[i]) / (time_m[i] * C_m[i] + abs(T_m_ext_dvar[i] * C_m_ext[i]))))) +
  (sum (i in k..I0)
    sum (j in 1..J0)
      ((1 - abs((time_exp[i] - time_c[i][j]) / time_exp[i])) * ((time_c[i][j] * C_c[i][j]) / (time_c[i][j] * C_c[i][j] + abs(time_exp[i] * C_c_ext[i][j])))))) /
  (k - 1 + (sum (j in 1..J0) (Ij[j] - k + 1))))) * 0.5 +
  ((sum (j in 1..J0)
    (abs((Tj_exp[j] - ((sum (i in 1..(k-1)) (time_m[i] + T_m_ext_dvar[i])) + (sum (i in k..Ij[j]) (time_c[i][j] + T_c_ext_dvar[i][j])))) / Tj_exp[j]))) / 3) * 0.5;

subject to{
  /*forall(j in 1..J0)
    total_time_ct:
    (sum (i in 1..(k-1)) (time_m[i] + T_m_ext_dvar[i]) + sum (i in k..Ij[j]) (time_c[i][j] + T_c_ext_dvar[i][j])) <= Tj_exp[j] * (1 + b);
*/
  forall(i in 1..(k-1))
    mass_time_delay_ct:
    T_m1[i + 1] <= time_exp[i] - time_m[i] - T_m_ext_dvar[i];

  forall(i in 1..(k-1))
    mass_time_ahead_ct:
    T_m2[i + 1] >= time_exp[i] - time_m[i] - T_m_ext_dvar[i];

  forall(i in k..(I0 - 1))
    forall(j in 1..J0)
    cus_time_delay_ct:
    T_c1[i + 1][j] <= time_exp[i] - time_c[i][j] - T_c_ext_dvar[i][j];

  forall(i in k..(I0 - 1))
    forall(j in 1..J0)
    cus_time_ahead_ct:
    T_c2[i + 1][j] >= time_exp[i] - time_c[i][j] - T_c_ext_dvar[i][j];

  forall(i in 1..(k-1))
    mass_satisfication_ct:
    U_m[i] <= (1 - abs((time_exp[i] - time_m[i]) / time_exp[i])) * ((time_m[i] * C_m[i]) / (time_m[i] * C_m[i] + abs(T_m_ext_dvar[i] * C_m_ext[i])));

  forall(i in k..I0)
    forall(j in 1..J0)
    cus_satisfication_ct:
    U_c[i][j] <= (1 - abs((time_exp[i] - time_c[i][j]) / time_exp[i])) * ((time_c[i][j] * C_c[i][j]) / (time_c[i][j] * C_c[i][j] + abs(T_c_ext_dvar[i][j] * C_c_ext[i][j])));

  Z_ct:
    ((sum (i in 1..(k-1))
    (time_m[i] * C_m[i] + abs(T_m_ext_dvar[i]) * C_m_ext[i])) * 3 +
  (sum (i in k..I0)
  	(sum (j in 1..J0)
      (time_c[i][j] * C_c[i][j] + abs(T_c_ext_dvar[i][j]) * C_c_ext[i][j]))) +
  (sum (i in 1..(k-1))
    (abs(time_exp[i] - time_m[i] - T_m_ext_dvar[i]) * P_m[i])) * 3 +
  (sum (i in k..I0)
    (sum (j in 1..J0)
      (abs(time_exp[i] - time_c[i][j] - T_c_ext_dvar[i][j]) * P_c[i][j])))) <=
  (sum (i in 1..(k-1))
    (time_m[i] * C_m[i]) * 3 +
    sum (i in k..I0)
      sum(j in 1..J0)
        (time_c[i][j] * C_c[i][j])) * (1 + c);

}
