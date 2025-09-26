import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt


def get_L2_errors(phi_coarse, phi_medium, phi_fine, hx_list, ht_list):
  """ Computes discrete L2 norm between two fields defined on the same grid. """

  hx_coarse, hx_medium, hx_fine = hx_list
  ht_coarse, ht_medium, ht_fine = ht_list

  # ratio in x grids ( should be 2)
  rx_cm = int(round(hx_coarse/hx_medium))
  rx_cf = int(round(hx_coarse/hx_fine))
  # print("rx_cm = ", rx_cm)
  # print("rx_cf = ", rx_cf)

  # ratio in t grids ( should be 1)
  # rt_cm = int(round(ht_coarse/ht_medium))
  # rt_cf = int(round(ht_coarse/ht_fine))
  # print("rt_cm = ", rt_cm)
  # print("rt_cf = ", rt_cf)

  # print("phi_coarse.shape = ", phi_coarse.shape)
  # print("phi_medium.shape = ", phi_medium.shape)
  # print("phi_medium[::rt_cm,::rx_cm].shape = ", phi_medium[:,::rx_cm].shape)
  
  delta_cm = (phi_coarse - phi_medium[:,::rx_cm])
  delta_mf = (phi_medium[:,::rx_cm] - phi_fine[:,::rx_cf])

  # print('rx_cm, rx_cf, rt_cm, rt_cf:')
  # print(rx_cm, rx_cf, rt_cm, rt_cf)

  L2_cm = np.sqrt(np.sum( delta_cm**2, axis=1) * hx_coarse )
  L2_mf = np.sqrt(np.sum( delta_mf**2, axis=1) * hx_coarse )
  return L2_cm, L2_mf


def get_L1_errors(phi_coarse, phi_medium, phi_fine, hx_list, ht_list):
  """ Computes discrete L2 norm between two fields defined on the same grid. """

  hx_coarse, hx_medium, hx_fine = hx_list
  ht_coarse, ht_medium, ht_fine = ht_list

  # ratio in x grids ( should be 2)
  rx_cm = int(round(hx_coarse/hx_medium))
  rx_cf = int(round(hx_coarse/hx_fine))

  # ratio in t grids ( should be 1)
  # rt_cm = int(round(ht_coarse/ht_medium))
  # rt_cf = int(round(ht_coarse/ht_fine))
  
  delta_cm = (phi_coarse - phi_medium[:,::rx_cm])
  delta_mf = (phi_medium[:,::rx_cm] - phi_fine[:,::rx_cf])

  # print('rx_cm, rx_cf, rt_cm, rt_cf:')
  # print(rx_cm, rx_cf, rt_cm, rt_cf)

  L1_cm = np.sum( np.abs(delta_cm), axis=1) * hx_coarse
  L1_mf = np.sum( np.abs(delta_mf), axis=1) * hx_coarse
  return L1_cm, L1_mf


def get_convergence_Q(L2_coarse_medium, L2_medium_fine, tol = 1e-14):
  return L2_coarse_medium/(L2_medium_fine + tol)


def get_convergence_Q_log2(L2_coarse_medium, L2_medium_fine, tol = 1e-14):
  return np.log2(L2_coarse_medium/(L2_medium_fine + tol))


def plot_convergence(problem, t_array, hx_list, ht_list, Q_list, var_name_list, var_label_list, frame, etaovers, ep_values, delta, Q_type):
  gct_dir = f"{base_dir}/gct/"
  if not os.path.exists(gct_dir):
    os.makedirs(gct_dir)

  hx_c, hx_m, hx_f = hx_list
  ht_c, ht_m, ht_f = ht_list
  
  plt.figure()

  for i in range(len(Q_list)):
    plt.plot(t_array, Q_list[i], label=var_label_list[i])

  plt.xlabel(r"$t$ (1/GeV)")
  plt.ylabel(r"Q(t)")

  # ymin = np.min(np.minimum(Q_list[0], Q_list[1])) - 0.2 * np.min(np.minimum(Q_list[0], Q_list[1]))
  # ymax = np.max(np.maximum(Q_list[0], Q_list[1])) + 0.2 * np.max(np.maximum(Q_list[0], Q_list[1]))
  # ymin = 
  # plt.ylim(ymin, ymax)

  plt.legend()

  plt.savefig(gct_dir + Q_type + "_" +problem+ "_" + ep_values +"_frame=" + str(frame) + "_delta=" + str(delta) +"_hx="+str(hx_c)+","+str(hx_m)+","+str(hx_f)+"_etas="+str(etaovers)+"_"+"_".join(var_name_list)+".pdf")


def info_from_unique(unique): 
  sep_all = unique.split("_")
  problem = sep_all[1]
  only_data = sep_all[2:]
  info = {  }
  info["problem"] = problem
  for data in only_data:
    split = data.split("=")
    info[split[0]] = split[1]
  print(info)
  return info


def get_data_for_conv_test(i, full_path_list, info_list):
  return {
    "t": np.loadtxt(f"{full_path_list[i]}/t.txt"),
    "x": np.loadtxt(f"{full_path_list[i]}/x.txt"),
    "T00": np.loadtxt(f"{full_path_list[i]}/T00.txt"),
    "T0x": np.loadtxt(f"{full_path_list[i]}/T0x.txt"),
    "ep": np.loadtxt(f"{full_path_list[i]}/ep.txt"),
    "v": np.loadtxt(f"{full_path_list[i]}/v.txt"),
    "hx": float(info_list[i]["hx"]),
    "ht": float(info_list[i]["ht"]),
  }


def info_from_unique(unique):
  sep_all = unique.split("_")
  problem = sep_all[1]
  only_data = sep_all[2:]
  info = {  }
  info["problem"] = problem
  for data in only_data:
    split = data.split("=")
    info[split[0]] = split[1]
  print(info)
  return info


################
##### MAIN #####
################

# YOU HAVE TO CHANGE THIS MANUALLY
# CONVENTION: [coarse sim, medium sim, fine sim]

base_dir = ""
unique_list = [ # list of directories where data was saved in "unique" format (see C++ code)
    "",
    "",
    "",
  ]
problem = ""
frame = 1
etaovs = 0.08
delta = 5
ep_values = "(epL,epR)=(1.48,1.06)"


# EVERYTHING BELOW THIS COMMENT IS AUTOMATIC

full_path_list = [f"{base_dir}/{u}" for u in unique_list]
info_list = [info_from_unique(u) for u in unique_list]

coarse = get_data_for_conv_test(0, full_path_list, info_list)
medium = get_data_for_conv_test(1, full_path_list, info_list)
fine = get_data_for_conv_test(2, full_path_list, info_list)
hx_list = [coarse["hx"], medium["hx"], fine["hx"]]
ht_list = [coarse["ht"], medium["ht"], fine["ht"]]

# ---------

T00_L2_cm, T00_L2_mf = get_L2_errors(coarse["T00"], medium["T00"], fine["T00"], hx_list, ht_list)
T0x_L2_cm, T0x_L2_mf = get_L2_errors(coarse["T0x"], medium["T0x"], fine["T0x"], hx_list, ht_list)
ep_L2_cm, ep_L2_mf = get_L2_errors(coarse["ep"], medium["ep"], fine["ep"], hx_list, ht_list)
v_L2_cm, v_L2_mf = get_L2_errors(coarse["v"], medium["v"], fine["v"], hx_list, ht_list)

Q_T00_log2_L2 = get_convergence_Q_log2(T00_L2_cm, T00_L2_mf)
Q_T0x_log2_L2 = get_convergence_Q_log2(T0x_L2_cm, T0x_L2_mf)
Q_ep_log2_L2_cm = get_convergence_Q_log2(ep_L2_cm, ep_L2_mf)
Q_v_log2_L2_cm = get_convergence_Q_log2(v_L2_cm, v_L2_mf)

Q_T00_L2 = get_convergence_Q(T00_L2_cm, T00_L2_mf)
Q_T0x_L2 = get_convergence_Q(T0x_L2_cm, T0x_L2_mf)
Q_ep_L2 = get_convergence_Q(ep_L2_cm, ep_L2_mf)
Q_v_L2 = get_convergence_Q(v_L2_cm, v_L2_mf)

Q_list_log2_L2 = [Q_T00_log2_L2, Q_T0x_log2_L2, Q_ep_log2_L2_cm, Q_v_log2_L2_cm]
Q_list_L2 = [Q_T00_L2, Q_T0x_L2, Q_ep_L2, Q_v_L2]

var_name_list = ["T00","T0x","ep","v"]
var_label_list = [r"$T^{00}(t,x)$",r"$T^{0x}(t,x)$",r"$\varepsilon(t,x)$",r"$v(t,x)$"]
t_array = coarse["t"]

# plotting
plot_convergence(problem, t_array, hx_list, ht_list, Q_list_L2, var_name_list, var_label_list, frame, etaovs, ep_values, delta, "L2")
plot_convergence(problem, t_array, hx_list, ht_list, Q_list_log2_L2, var_name_list, var_label_list, frame, etaovs, ep_values, delta, "L2_log2")
# ---------

T00_L1_cm, T00_L1_mf = get_L1_errors(coarse["T00"], medium["T00"], fine["T00"], hx_list, ht_list)
T0x_L1_cm, T0x_L1_mf = get_L1_errors(coarse["T0x"], medium["T0x"], fine["T0x"], hx_list, ht_list)
ep_L1_cm, ep_L1_mf = get_L1_errors(coarse["ep"], medium["ep"], fine["ep"], hx_list, ht_list)
v_L1_cm, v_L1_mf = get_L1_errors(coarse["v"], medium["v"], fine["v"], hx_list, ht_list)

Q_T00_log2_L1 = get_convergence_Q_log2(T00_L1_cm, T00_L1_mf)
Q_T0x_log2_L1 = get_convergence_Q_log2(T0x_L1_cm, T0x_L1_mf)
Q_ep_log2_L1_cm = get_convergence_Q_log2(ep_L1_cm, ep_L1_mf)
Q_v_log2_L1_cm = get_convergence_Q_log2(v_L1_cm, v_L1_mf)

Q_T00_L1 = get_convergence_Q(T00_L1_cm, T00_L1_mf)
Q_T0x_L1 = get_convergence_Q(T0x_L1_cm, T0x_L1_mf)
Q_ep_L1 = get_convergence_Q(ep_L1_cm, ep_L1_mf)
Q_v_L1 = get_convergence_Q(v_L1_cm, v_L1_mf)

Q_list_log2_L1 = [Q_T00_log2_L1, Q_T0x_log2_L1, Q_ep_log2_L1_cm, Q_v_log2_L1_cm]
Q_list_L1 = [Q_T00_L1, Q_T0x_L1, Q_ep_L1, Q_v_L1]

var_name_list = ["T00","T0x","ep","v"]
var_label_list = [r"$T^{00}(t,x)$",r"$T^{0x}(t,x)$",r"$\varepsilon(t,x)$",r"$v(t,x)$"]
t_array = coarse["t"]

# plotting
plot_convergence(problem, t_array, hx_list, ht_list, Q_list_L1, var_name_list, var_label_list, frame, etaovs, ep_values, delta, "L1")
plot_convergence(problem, t_array, hx_list, ht_list, Q_list_log2_L1, var_name_list, var_label_list, frame, etaovs, ep_values, delta, "L1_log2")