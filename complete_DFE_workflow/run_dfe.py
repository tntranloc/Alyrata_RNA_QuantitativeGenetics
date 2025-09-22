import numpy as np
import sys
import pickle
import scipy
#import matplotlib.pyplot as plt
#import seaborn as sns
import os
import operator
import pandas as pd
from tqdm import tqdm
from numpy import logical_and, logical_not
from scipy.special import gammaln
import scipy.stats.distributions
import scipy.integrate
import scipy.optimize
from dadi import Numerics, Inference, Misc
from dadi.Spectrum_mod import Spectrum
import dadi
from scipy.stats import gamma
from copy import deepcopy
from scipy.optimize import basinhopping

# Set up plotting parameters (currently commented out)
# sns.set_theme(style="whitegrid", context="paper")
# plt.rcParams.update({
#     "figure.dpi": 600,
#     "savefig.dpi": 600,
#     "grid.alpha": 0.3,
#     "grid.linestyle": "--",
#     "axes.edgecolor": "#444444",
#     "axes.linewidth": 0.8,
#     "font.family": "serif"
# })

# Import dadi functions
import dadi
from dadi import logging
logging.basicConfig()

# Define demographic parameters and other constants
Nanc = 2846301.0 / 2
Nbot = 1708210.0 / 2
SPN2 = 13218.0 / 2
PLN2 = 898319.0 / 2
SPN1 = 13218.0 / 2
PLN1 = 101220.0 / 2
TPLSP = 26902.0
TPLAUS = 372701.0
M_SPPL = 0.415
M_PLSP = 0.036
T_ISOL = 24726
TB = (TPLAUS - TPLSP) / (2 * Nanc)
TF = (TPLSP - T_ISOL) / (2 * Nanc)
TI = T_ISOL / (2 * Nanc)
nuB = Nbot / Nanc
SPnuF2 = SPN2 / Nanc
PLnuF2 = PLN2 / Nanc
SPnuF1 = SPN1 / Nanc
PLnuF1 = PLN1 / Nanc
SP_PL_theta_ratio = 0.9023831

# Define the spectra class and necessary functions
class spectra:
    def __init__(self, params, ns, demo_sel_func, pts=500, pts_l=None, Npts=500, n=20170., int_breaks=None,
                 int_bounds=(1e-4, 1000.), mp=False, echo=False, cpus=None):
        self.ns = ns
        self.spectra = []
        if int_breaks is not None:
            numbreaks = len(int_breaks)
            stepint = Npts / (numbreaks - 1)
            self.gammas = []
            for i in reversed(range(0, numbreaks - 1)):
                self.gammas = np.append(
                    self.gammas, -np.logspace(np.log10(int_breaks[i + 1]),
                                              np.log10(int_breaks[i]),
                                              stepint))
        else:
            self.gammas = -np.logspace(np.log10(int_bounds[1]),
                                       np.log10(int_bounds[0]), Npts)
        if pts_l is None:
            self.pts = pts
            self.pts_l = [self.pts, self.pts + (self.pts / 5), self.pts + (self.pts / 5) * 2]
        else:
            self.pts_l = pts_l
        self.params = tuple(params)
        for ii, gamma in enumerate(tqdm(self.gammas)):
            self.spectra.append(top_level_extrap_func(tuple(params) + (gamma,), self.ns, self.pts_l))
            if echo:
                print(f'{ii}: {gamma}')
        self.neu_spec = top_level_extrap_func(tuple(params) + (0,), self.ns, self.pts_l)
        self.extrap_x = self.spectra[0].extrap_x
        self.spectra = np.array(self.spectra)

    def integrate(self, params, sel_dist, theta):
        sel_args = (self.gammas,) + tuple(params)
        weights = sel_dist(*sel_args)
        weight_neu, err_neu = scipy.integrate.quad(sel_dist, self.gammas[-1], 0, args=tuple(params))
        pops = len(self.neu_spec.shape)
        if pops == 1:
            integrated = self.neu_spec * weight_neu + Numerics.trapz(
                weights[:, np.newaxis] * self.spectra, self.gammas, axis=0)
        elif pops == 2:
            integrated = self.neu_spec * weight_neu + Numerics.trapz(
                weights[:, np.newaxis, np.newaxis] * self.spectra, self.gammas, axis=0)
        elif pops == 3:
            integrated = self.neu_spec * weight_neu + Numerics.trapz(
                weights[:, np.newaxis, np.newaxis, np.newaxis] * self.spectra, self.gammas, axis=0)
        else:
            raise IndexError("Must have one to three populations")
        integrated_fs = Spectrum(integrated, extrap_x=self.extrap_x)
        return integrated_fs * theta

def top_level_extrap_func(params, ns, pts_l):
    return three_epoch_growth_sel_new(params, ns, pts_l[0])

def three_epoch_growth_sel_new(params, ns, pts):
    nuB, nuF2, nuF1, TB, TF, TI, pi, gamma = params
    T1 = (TB + TF + TI) * pi
    T2 = (TB + TF + TI) * (1 - pi)
    xx = dadi.Numerics.default_grid(pts)
    min_gamma = -500
    if gamma < min_gamma:
        phi = dadi.PhiManip.phi_1D(xx, gamma=min_gamma)
    else:
        phi = dadi.PhiManip.phi_1D(xx, gamma=gamma)
    if nuB * gamma < min_gamma:
        phi = dadi.Integration.one_pop(phi, xx, T1, nuB, gamma=min_gamma)
    else:
        phi = dadi.Integration.one_pop(phi, xx, T1, nuB, gamma=gamma)
    if gamma == 0:
        nu_func = lambda t: nuF2 * (nuF1 / nuF2) ** (t / TI)
    else:
        nu_func = lambda t: np.where(gamma * nuF2 * (nuF1 / nuF2) ** (t / T2) < min_gamma,
                                     np.abs(min_gamma / gamma),
                                     nuF2 * (nuF1 / nuF2) ** (t / T2))
    phi = dadi.Integration.one_pop(phi, xx, T2, nu_func, gamma=gamma)
    fs = dadi.Spectrum.from_phi(phi, (28,), (xx,))
    return fs

def three_epoch_growth_sel_both_new(params, ns, pts):
    min_gamma = -500
    nuB, PLnuF2, PLnuF1, SPnuF2, SPnuF1, TB, TF, TI, pi, gamma = params
    T1 = (TB + TF + TI) * pi
    T2 = (TB + TF + TI) * (1 - pi)
    xx = dadi.Numerics.default_grid(pts)
    if gamma < min_gamma:
        phi = dadi.PhiManip.phi_1D(xx, gamma=min_gamma)
    else:
        phi = dadi.PhiManip.phi_1D(xx, gamma=gamma)
    if nuB * gamma < min_gamma:
        phi_2 = dadi.Integration.one_pop(phi, xx, T1, nuB, gamma=min_gamma)
    else:
        phi_2 = dadi.Integration.one_pop(phi, xx, T1, nuB, gamma=gamma)
    phi_3 = dadi.PhiManip.phi_1D_to_2D(xx, phi_2)
    if gamma == 0:
        nu_func1 = lambda t: PLnuF2 * (PLnuF1 / PLnuF2) ** (t / T2)
        nu_func2 = lambda t: SPnuF2 * (SPnuF1 / SPnuF2) ** (t / T2)
    else:
        nu_func1 = lambda t: np.where(gamma * PLnuF2 * (PLnuF1 / PLnuF2) ** (t / T2) < min_gamma,
                                     np.abs(min_gamma / gamma),
                                     PLnuF2 * (PLnuF1 / PLnuF2) ** (t / T2))
        nu_func2 = lambda t: np.where(gamma * SPnuF2 * (SPnuF1 / SPnuF2) ** (t / T2) < min_gamma,
                                     np.abs(min_gamma / gamma),
                                     SPnuF2 * (SPnuF1 / SPnuF2) ** (t / T2))
    phi_4 = dadi.Integration.two_pops(phi_3, xx, T2, nu1=nu_func1, nu2=nu_func2, gamma1=gamma, gamma2=gamma)
    fs = dadi.Spectrum.from_phi(phi_4, (28, 28), (xx, xx), pop_ids=["PL", "SP"])
    return fs

def gamma_dist(mgamma, alpha, beta):
    return scipy.stats.distributions.gamma.pdf(-mgamma, alpha, scale=beta)

def three_epoch_growth_neut_SP_opt(params, pts):
    nuB, PLnu, SPnuF2, SPnuF1, pi = params
    xx = dadi.Numerics.default_grid(1000)
    T1 = (TB + TF + TI) * pi
    T2 = (TB + TF + TI) * (1 - pi)
    phi = dadi.PhiManip.phi_1D(xx)
    phi_2 = dadi.Integration.one_pop(phi, xx, T1, nuB)
    phi_3_PL = dadi.Integration.one_pop(phi_2, xx, T2, PLnu)
    nu_funcSP = lambda t: SPnuF2 * (SPnuF1 / SPnuF2) ** (t / T2)
    phi_4_SP = dadi.Integration.one_pop(phi_2, xx, T2, nu_funcSP)
    fs_PL = dadi.Spectrum.from_phi(phi_3_PL, (28,), (xx,))
    fs_SP = dadi.Spectrum.from_phi(phi_4_SP, (28,), (xx,))
    return fs_PL, fs_SP

def obj_funct_SP_opt(params, data_PL_4fold, data_SP_4fold):
    demo_params = np.concatenate([10**np.array([params[0], params[1], params[2], params[3]]), [params[4]]])
    sfs_PL, sfs_SP = three_epoch_growth_neut_SP_opt(demo_params, 400)
    sfs_PL_fold = sfs_PL.fold() * 10**params[5]
    sfs_SP_fold = sfs_SP.fold() * 10**params[5] * SP_PL_theta_ratio
    result = -Inference.ll(sfs_PL_fold, data_PL_4fold) - Inference.ll(sfs_SP_fold, data_SP_4fold)
    return result

def optimize_parameters(data_PL_4fold, data_SP_4fold, start_params):
    optimized_result = basinhopping(
        obj_funct_SP_opt,
        start_params,
        T=100,
        stepsize=5,
        interval=5,
        disp=True,
        niter=10,
        minimizer_kwargs={
            "method": "L-BFGS-B",
            "bounds": [(-1, 1), (-3, 1), (-4, 3), (-4, 3), (0.001, 0.999), (4, np.log10(5e5))],
            "options": {"disp": False, "maxcor": 20},
            "args": (data_PL_4fold, data_SP_4fold)  # Pass the data as arguments
        }
    )
    return optimized_result.x

def process_bootstrap(start_idx, end_idx, folder_path):
    ns_SP = np.array([28])
    ns_PL = np.array([28])

    for boot_idx in range(start_idx, end_idx + 1):
        print(f"Processing bootstrap {boot_idx}...")
        data_SP_0fold = dadi.Spectrum.from_file(f'{folder_path}/SP_VA_0fold_boot_{boot_idx}.sfs')
        data_PL_0fold = dadi.Spectrum.from_file(f'{folder_path}/PL_VA_0fold_boot_{boot_idx}.sfs')
        data_SP_4fold = dadi.Spectrum.from_file(f'{folder_path}/SP_VA_4fold_boot_{boot_idx}.sfs')
        data_PL_4fold = dadi.Spectrum.from_file(f'{folder_path}/PL_VA_4fold_boot_{boot_idx}.sfs')

        # Optimize parameters using basinhopping
        start_sfs_PL, _ = three_epoch_growth_neut_SP_opt([nuB, PLnuF2, SPnuF2, SPnuF1, (TB)/(TF+TI+TB)], 1000)
        start_theta = dadi.Inference.optimal_sfs_scaling(start_sfs_PL, data_PL_4fold)
        start_params = [np.log10(nuB), np.log10(PLnuF2), np.log10(SPnuF2), np.log10(SPnuF1), 0.1, np.log10(start_theta)]
        params = optimize_parameters(data_PL_4fold, data_SP_4fold, start_params)

        nuB_opt, PLnu_opt, SPnuF2_opt, SPnuF1_opt, _, theta_PL_opt = 10**np.array(params)
        opt = params[4]
        T1 = (TB + TF + TI) * opt
        T2 = (TB + TF + TI) * (1 - opt)
        demo_params_opt_PL = [nuB_opt, PLnu_opt, PLnu_opt, TB, TF, TI, opt]
        demo_params_opt_SP = [nuB_opt, SPnuF2_opt, SPnuF1_opt, TB, TF, TI, opt]
        demo_params_opt = [nuB_opt, PLnu_opt, PLnu_opt, SPnuF2_opt, SPnuF1_opt, TB, TF, TI, opt]

        theta_PL = theta_PL_opt
        theta_SP = theta_PL_opt * SP_PL_theta_ratio
        theta_SP_opt_ns = theta_SP * 2.76
        theta_PL_opt_ns = theta_PL * 2.76

        spectra_SP_opt = spectra(demo_params_opt_SP, ns_SP, three_epoch_growth_sel_new, pts=1000,
                                  int_bounds=(0.01, 1e5), Npts=500, echo=False, mp=True)
        spectra_PL_opt = spectra(demo_params_opt_PL, ns_PL, three_epoch_growth_sel_new, pts=1000,
                                  int_bounds=(0.01, 1e5), Npts=500, echo=False, mp=True)

        SP_shape_set = np.logspace(start=-1, stop=1.0, num=250)
        SP_scale_fit_set = np.ones_like(SP_shape_set)
        SP_llhood_fit_set = np.ones_like(SP_shape_set)
        for ii, shape in enumerate(SP_shape_set):
            def obj_funct(scale):
                sfs_SP = spectra_SP_opt.integrate([shape, scale], gamma_dist, 1).fold()
                return -Inference.ll_multinom(sfs_SP, data_SP_0fold)
            fit = scipy.optimize.minimize(obj_funct, [SP_scale_fit_set[ii - 1]], method="SLSQP",
                                          bounds=[(np.exp(-8), np.exp(20))])
            SP_scale_fit_set[ii] = fit["x"][0]
            SP_llhood_fit_set[ii] = fit["fun"]

        SP_shape_set_t = np.logspace(start=-2, stop=0, num=250)
        SP_scale_fit_set_t = np.ones_like(SP_shape_set_t)
        SP_llhood_fit_set_t = np.ones_like(SP_shape_set_t)
        for ii, shape in enumerate(SP_shape_set_t):
            def obj_funct(scale):
                sfs_SP = spectra_SP_opt.integrate([shape, np.exp(scale)], gamma_dist, theta_SP_opt_ns).fold()
                return -Inference.ll(sfs_SP, data_SP_0fold)
            fit = scipy.optimize.minimize(obj_funct, [np.log(100)], method="SLSQP", bounds=[(-8, 20)])
            SP_scale_fit_set_t[ii] = np.exp(fit["x"][0])
            SP_llhood_fit_set_t[ii] = fit["fun"]

        PL_shape_set = np.logspace(start=-1, stop=1, num=250)
        PL_scale_fit_set = np.ones_like(PL_shape_set)
        PL_llhood_fit_set = np.ones_like(PL_shape_set)
        for ii, shape in enumerate(PL_shape_set):
            def obj_funct(scale):
                sfs_PL = spectra_PL_opt.integrate([shape, scale], gamma_dist, 1).fold()
                return -Inference.ll_multinom(sfs_PL, data_PL_0fold)
            fit = scipy.optimize.minimize(obj_funct, [PL_scale_fit_set[ii - 1]], method="SLSQP",
                                          bounds=[(np.exp(-8), np.exp(20))])
            PL_scale_fit_set[ii] = fit["x"][0]
            PL_llhood_fit_set[ii] = fit["fun"]

        PL_shape_set_t = np.logspace(start=-2, stop=0, num=250)
        PL_scale_fit_set_t = np.ones_like(PL_shape_set_t)
        PL_llhood_fit_set_t = np.ones_like(PL_shape_set_t)
        for ii, shape in enumerate(PL_shape_set_t):
            def obj_funct(scale):
                sfs_PL = spectra_PL_opt.integrate([shape, np.exp(scale)], gamma_dist, theta_PL_opt_ns).fold()
                return -Inference.ll(sfs_PL, data_PL_0fold)
            fit = scipy.optimize.minimize(obj_funct, [np.log(100)], method="L-BFGS-B", bounds=[(-8, 20)])
            PL_scale_fit_set_t[ii] = np.exp(fit["x"][0])
            PL_llhood_fit_set_t[ii] = fit["fun"]

        shape_set = np.logspace(start=-1, stop=1.0, num=250)
        scale_fit_set = np.ones_like(shape_set)
        llhood_fit_set = np.ones_like(shape_set)
        for ii, shape in enumerate(shape_set):
            def obj_funct(scale):
                sfs_SP = spectra_SP_opt.integrate([shape, scale], gamma_dist, 1).fold()
                sfs_PL = spectra_PL_opt.integrate([shape, scale], gamma_dist, 1).fold()
                return -Inference.ll_multinom(sfs_SP, data_SP_0fold) - Inference.ll_multinom(sfs_PL, data_PL_0fold)
            fit = scipy.optimize.minimize(obj_funct, [scale_fit_set[ii - 1]], method="SLSQP",
                                          bounds=[(np.exp(-8), np.exp(20))])
            scale_fit_set[ii] = fit["x"][0]
            llhood_fit_set[ii] = fit["fun"]

        shape_set_t = np.logspace(start=-2, stop=0, num=250)
        scale_fit_set_t = np.ones_like(shape_set_t)
        llhood_fit_set_t = np.zeros_like(shape_set_t)
        for ii, shape in enumerate(shape_set_t):
            def obj_funct(scale):
                sfs_SP = spectra_SP_opt.integrate([shape, np.exp(scale)], gamma_dist, theta_SP_opt_ns).fold()
                sfs_PL = spectra_PL_opt.integrate([shape, np.exp(scale)], gamma_dist, theta_PL_opt_ns).fold()
                return -Inference.ll(sfs_SP, data_SP_0fold) - Inference.ll(sfs_PL, data_PL_0fold)
            fit = scipy.optimize.minimize(obj_funct, [np.log(100)], method="SLSQP", bounds=[(-8, 20)])
            scale_fit_set_t[ii] = np.exp(fit["x"][0])
            llhood_fit_set_t[ii] = fit["fun"]

        gamma_set = np.logspace(-1, 4, 100)
        spectra_full_sel = []
        for ii, gamma_tmp in enumerate(tqdm(gamma_set)):
            spectra_full_sel.append(three_epoch_growth_sel_both_new(demo_params_opt + [-gamma_tmp], 1, 100))

        F_gamma_both = gamma.cdf(gamma_set, a=shape_set_t[np.argmin(llhood_fit_set_t)],
                                 scale=scale_fit_set_t[np.argmin(llhood_fit_set_t)])
        p_gamma_both = F_gamma_both - np.concatenate([[0], F_gamma_both[0:-1]])
        p_gamma_both[-1] += 1 - p_gamma_both[-1]

        PL_cDFE_set = []
        SP_cDFE_set = []
        for ii in range(14):
            cDFE_PL_tmp = np.zeros_like(gamma_set)
            cDFE_SP_tmp = np.zeros_like(gamma_set)
            for jj, gamma_tmp in enumerate(gamma_set):
                cDFE_PL_tmp[jj] = p_gamma_both[jj] * spectra_full_sel[jj].marginalize([1]).fold()[ii + 1]
                cDFE_SP_tmp[jj] = p_gamma_both[jj] * spectra_full_sel[jj].marginalize([0]).fold()[ii + 1]
            PL_cDFE_set.append(deepcopy(cDFE_PL_tmp / np.sum(cDFE_PL_tmp)))
            SP_cDFE_set.append(deepcopy(cDFE_SP_tmp / np.sum(cDFE_SP_tmp)))

        gamma_bins = np.digitize(gamma_set, np.array([-1, 0.1, 1, 10, 100, 1e6]), right=True)
        bins = np.arange(1, 6)
        PL_bin_prob_set = []
        SP_bin_prob_set = []
        for ii in range(14):
            PL_bin_probs = [0] * 5
            SP_bin_probs = [0] * 5
            for jj, bin_ind in enumerate(bins):
                PL_bin_probs[jj] = np.sum(PL_cDFE_set[ii][gamma_bins == bin_ind])
                SP_bin_probs[jj] = np.sum(SP_cDFE_set[ii][gamma_bins == bin_ind])
            PL_bin_probs = np.array(PL_bin_probs)
            SP_bin_probs = np.array(SP_bin_probs)
            PL_bin_probs[np.array(PL_bin_probs) < 0] = 0
            SP_bin_probs[np.array(SP_bin_probs) < 0] = 0
            PL_bin_prob_set.append(deepcopy(PL_bin_probs) / np.sum(PL_bin_probs))
            SP_bin_prob_set.append(deepcopy(SP_bin_probs) / np.sum(SP_bin_probs))

        PL_bin_prob_set = np.array(PL_bin_prob_set)
        SP_bin_prob_set = np.array(SP_bin_prob_set)

        PL_fractions = np.sum(PL_bin_prob_set, axis=0)
        SP_fractions = np.sum(SP_bin_prob_set, axis=0)

        PL_fractions_norm = PL_fractions / np.sum(PL_fractions)
        SP_fractions_norm = SP_fractions / np.sum(SP_fractions)

        # Save results to CSV
        results = {
            "Fitness Effect Class": [
                "(gamma <= 0.1)",
                "(0.1 < gamma <= 1)",
                "(1 < gamma <= 10)",
                "(10 < gamma < 100)",
                "(gamma > 100)"
            ],
            "PL Fractions (Raw)": PL_fractions,
            "PL Fractions (Normalised)": PL_fractions_norm,
            "SP Fractions (Raw)": SP_fractions,
            "SP Fractions (Normalised)": SP_fractions_norm
        }

        df = pd.DataFrame(results)
        df.to_csv(f'{folder_path}/results_boot_{boot_idx}.csv', index=False)
        print(f"Results for bootstrap {boot_idx} saved to CSV.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process bootstrap SFS results.")
    parser.add_argument("--start", type=int, required=True, help="Start index of bootstrap files.")
    parser.add_argument("--end", type=int, required=True, help="End index of bootstrap files.")
    parser.add_argument("--folder", type=str, required=True, help="Folder path containing bootstrap files.")
    args = parser.parse_args()

    process_bootstrap(args.start, args.end, args.folder)
