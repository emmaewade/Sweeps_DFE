"""
Uses fitdadi to infer DFE.

"""


import sys
import os
import csv
import json
import logging
import time
import argparse
import warnings
import re

import ast
import numpy as np
import nlopt
import dadi
import dadi.DFE as DFE
from dadi.DFE import *
import scipy.stats.distributions
import scipy.integrate
import scipy.optimize


class ArgumentParserNoArgHelp(argparse.ArgumentParser):
    """Like *argparse.ArgumentParser*, but prints help when no arguments."""

    def error(self, message):
        """Print error message, then help."""
        sys.stderr.write('error: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


class DFEInference():
    """Wrapper class to allow functions to reference each other."""

    def ExistingFile(self, fname):
        """Return *fname* if existing file, otherwise raise ValueError."""
        if os.path.isfile(fname):
            return fname
        else:
            raise ValueError("%s must specify a valid file name" % fname)

    def neugamma(self, xx, params):
        """Return Neutral + Gamma distribution (unused while the NEUGAMMA
        inference blocks below are commented out)."""
        pneu, alpha, beta = params
        xx = np.atleast_1d(xx)
        out = (1-pneu)*DFE.PDFs.gamma(xx, (alpha, beta))
        # Assume gamma < 1e-4 is essentially neutral
        out[np.logical_and(0 <= xx, xx < 1e-4)] += pneu/1e-4
        # Reduce xx back to scalar if it's possible
        return np.squeeze(out)

    def fitDFEParser(self):
        """Return *argparse.ArgumentParser* for ``fitdadi_infer_DFE.py``."""
        parser = ArgumentParserNoArgHelp(
            description=(
                'Given inferred demographic model parameters, this script '
                'infers a maximum-likelihood distriubtion of fitness '
                'effects.'),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument(
            'syn_input_sfs', type=self.ExistingFile,
            help=('Synonynomous site-frequency spectrum from which the '
                  'demographic parameters were inferred.'))
        parser.add_argument(
            'nonsyn_input_sfs', type=self.ExistingFile,
            help=('Nonsynonynomous site-frequency spectrum from which the '
                  'distribution of fitness effects should be inferred.'))
        parser.add_argument(
            'input_demography', type=self.ExistingFile,
            help=('Inferred demographic parameters output by Dadi.'))
        parser.add_argument(
            'input_model', type=str,
            help=('Input demographic model.'))
        parser.add_argument(
            '--initial_alpha', type=float,
            dest='initial_alpha',
            help=('Float value for initial alpha of gamma-DFE.'),
            default=0.2)
        parser.add_argument(
            '--initial_beta', type=float,
            dest='initial_beta',
            help=('Float value for initial beta of gamma-DFE.'),
            default=0.1)
        parser.add_argument(
            '--L_syn', type=float,
            dest='L_syn',
            help=('Sum of allele counts for synonymous SFS.'),
            default=8058343.0)
        parser.add_argument(
            '--fold', dest='fold', action='store_true',
            help=('Fold the input SFSs and the outputted model SFS.'))
        parser.add_argument(
            '--no-fold', dest='fold', action='store_false',
            help=('Use unfolded SFSs throughout.'))
        parser.set_defaults(fold=True)
        parser.add_argument(
            '--likelihood', choices=['poisson', 'multinom'], default='poisson',
            help=('poisson: theta_nonsyn fixed at theta_syn * ns_s_ratio, '
                  'model SFS not rescaled. multinom: model SFS rescaled to '
                  'the observed total (theta_nonsyn ignored).'))
        parser.add_argument(
            '--ns_s_ratio', type=float, required=True,
            help=('Ratio of the nonsynonymous to synonymous mutational '
                  'target within the analysed window (theta_nonsyn / '
                  'theta_syn).'))
        parser.add_argument(
            '--pts', type=int, nargs=3, default=[1200, 1600, 2000],
            help=('dadi grid sizes used to build the selection cache.'))
        parser.add_argument(
            '--shape_bounds', type=float, nargs=2, default=[1e-3, 100.],
            help=('Lower/upper bound on the gamma shape parameter.'))
        parser.add_argument(
            '--scale_bounds', type=float, nargs=2, default=[1e-2, 1e6],
            help=('Lower/upper bound on the gamma scale parameter, in '
                  'population-scaled (2 * Na * |s|) units.'))
        parser.add_argument(
            '--cpus', type=int,
            default=int(os.environ.get('NSLOTS', '1')),
            help=('CPUs for building the selection cache.'))
        parser.add_argument(
            'outprefix', type=str,
            help='The file prefix for the output files')
        return parser

    def main(self):
        """Execute main function."""
        # Parse command line arguments
        parser = self.fitDFEParser()
        args = vars(parser.parse_args())
        prog = parser.prog

        # Assign arguments
        syn_input_sfs = args['syn_input_sfs']
        nonsyn_input_sfs = args['nonsyn_input_sfs']
        outprefix = args['outprefix']
        input_demography = args['input_demography']
        input_model = args['input_model']
        initial_alpha = args['initial_alpha']
        initial_beta = args['initial_beta']
        L_syn = args['L_syn']
        fold = args['fold']
        likelihood = args['likelihood']
        use_multinom = likelihood == 'multinom'
        ns_s_ratio = args['ns_s_ratio']
        pts_l = list(args['pts'])
        cpus = args['cpus']
        lower_bound = [args['shape_bounds'][0], args['scale_bounds'][0]]
        upper_bound = [args['shape_bounds'][1], args['scale_bounds'][1]]

        # Numpy options
        np.set_printoptions(linewidth=np.inf)

        # create output directory if needed
        outdir = os.path.dirname(args['outprefix'])
        if outdir:
            if not os.path.isdir(outdir):
                if os.path.isfile(outdir):
                    os.remove(outdir)
                os.mkdir(outdir)

        # Output files: logfile, *demography.txt, inferred_DFE.txt
        # Remove output files if they already exist
        underscore = '' if args['outprefix'][-1] == '/' else '_'
        inferred_DFE = \
           '{0}{1}inferred_DFE.txt'.format(
                args['outprefix'], underscore)
        logfile = '{0}{1}log.log'.format(args['outprefix'], underscore)
        # CSV summary output, one row per DFE model fitted.
        dfe_results_csv = \
            '{0}{1}dfe_results.csv'.format(args['outprefix'], underscore)
        # File paths for saved empirical SFS copies.
        empirical_syn_sfs_file = \
            '{0}{1}empirical_syn_sfs.fs'.format(args['outprefix'], underscore)
        empirical_nonsyn_sfs_file = \
            '{0}{1}empirical_nonsyn_sfs.fs'.format(
                args['outprefix'], underscore)
        to_remove = [logfile, inferred_DFE, dfe_results_csv]
        for f in to_remove:
            if os.path.isfile(f):
                os.remove(f)

        # Set up to log everything to logfile.
        logging.shutdown()
        logging.captureWarnings(True)
        logging.basicConfig(
            format='%(asctime)s - %(levelname)s - %(message)s',
            level=logging.INFO)
        logger = logging.getLogger(prog)
        warning_logger = logging.getLogger("py.warnings")
        logfile_handler = logging.FileHandler(logfile)
        logger.addHandler(logfile_handler)
        warning_logger.addHandler(logfile_handler)
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s')
        logfile_handler.setFormatter(formatter)
        logger.setLevel(logging.INFO)

        # print some basic information
        logger.info('Beginning execution of {0} in directory {1}\n'.format(
            prog, os.getcwd()))
        logger.info('Progress is being logged to {0}\n'.format(logfile))
        logger.info('Parsed the following arguments:\n{0}\n'.format(
            '\n'.join(['\t{0} = {1}'.format(*tup) for tup in args.items()])))
        logger.info('Fold setting: {0}\n'.format(fold))
        logger.info('Likelihood: {0}; theta_nonsyn/theta_syn ratio: {1}; '
                    'pts: {2}; cpus: {3}\n'.format(
                        likelihood, ns_s_ratio, pts_l, cpus))
        logger.info('Gamma bounds (shape, scaled scale): lower {0}, '
                    'upper {1}\n'.format(lower_bound, upper_bound))

        # Construct initial Spectrum object from input sfs's.
        syn_data = dadi.Spectrum.from_file(syn_input_sfs)
        nonsyn_data = dadi.Spectrum.from_file(nonsyn_input_sfs)

        # Fold the input SFSs if requested. Guarded so an already-folded
        # input isn't folded a second time, which dadi disallows.
        if fold:
            if not syn_data.folded:
                syn_data = syn_data.fold()
            if not nonsyn_data.folded:
                nonsyn_data = nonsyn_data.fold()

        # Save copies of the empirical SFS actually used for this
        # inference run, alongside the DFE results, for provenance.
        syn_data.to_file(empirical_syn_sfs_file)
        nonsyn_data.to_file(empirical_nonsyn_sfs_file)
        logger.info('Saved empirical synonymous SFS to {0}\n'.format(
            empirical_syn_sfs_file))
        logger.info('Saved empirical nonsynonymous SFS to {0}\n'.format(
            empirical_nonsyn_sfs_file))
        logger.info('Empirical SFS folded: syn = {0}, nonsyn = {1}'.format(
            syn_data.folded, nonsyn_data.folded))

        syn_ns = syn_data.sample_sizes # Number of samples.
        nonsyn_ns = nonsyn_data.sample_sizes # Number of samples

        logger.info('Loading input demography.')
        with open(input_demography, 'r') as f:
            lines = [line for line in f]

        for line in lines:
            if 'Best fit parameters:' in line:
                # Parse every number inside the brackets. The previous
                # r"\d+\.\d+" pattern misread numpy's scientific notation
                # (e.g. "4.29e+01" -> 4.29) and found nothing in "1.e+02".
                bracket = line[line.index('[') + 1:line.index(']')]
                demog_params = [float(param) for param in
                                bracket.replace(',', ' ').split()]
            if 'Optimal value of theta_syn:' in line:
                theta_syn = str(line.split(': ')[1])
                theta_syn = theta_syn[0:-2]
                theta_syn = float(theta_syn)
        logger.info('Input demographic parameters are: ' +
                    str(demog_params) + '.')
        n_expected = {'two_epoch': 2, 'three_epoch': 4}.get(input_model)
        if n_expected is not None and len(demog_params) != n_expected:
            logger.error('Expected {0} demographic parameters for {1}, '
                         'parsed {2} from {3}.'.format(
                             n_expected, input_model, demog_params,
                             input_demography))
            sys.exit(1)
        logger.info('Input theta_syn is: ' + str(theta_syn) + '.')
        theta_nonsyn = theta_syn * ns_s_ratio
        logger.info('theta_nonsyn = {0} (theta_syn {1} x ratio {2}).'.format(
            theta_nonsyn, theta_syn, ns_s_ratio))
        if min(pts_l) <= nonsyn_ns[0]:
            logger.error('Grid pts {0} must exceed the sample size {1}.'.format(
                pts_l, nonsyn_ns[0]))
            sys.exit(1)
        logger.info('Generating spectrum from input demography.')
        # input_model = 'two_epoch'
        if input_model == 'two_epoch':
            print('Demographic parameters are: ' + str(demog_params) + '.')
            spectra = DFE.Cache1D(demog_params, nonsyn_ns, DFE.DemogSelModels.two_epoch_sel,
                                  pts=pts_l, gamma_bounds=(1e-5, 500), gamma_pts=100,
                                  verbose=True, cpus=cpus)
        elif input_model == 'three_epoch':
            spectra = DFE.Cache1D(demog_params, nonsyn_ns, DFE.DemogSelModels.three_epoch_sel,
                                  pts=pts_l, gamma_bounds=(1e-5, 500), gamma_pts=100,
                                  verbose=True, cpus=cpus)
        else:
            spectra = DFE.Cache1D(demog_params, nonsyn_ns, dadi.Demographics1D.snm,
                                  pts=pts_l, gamma_bounds=(1e-5, 500), gamma_pts=100,
                                  verbose=True, cpus=cpus)

        # Diagnostic for NaN entries in the cached selection spectra.
        # One NaN gamma bin makes the integrated model NaN for every DFE
        # parameter set. Wrapped in try/except so a change in dadi's
        # attribute names can't crash the run.
        try:
            cache_arr = np.asarray(spectra.spectra, dtype=float)[:, 1:-1]
            bad_bins = np.isnan(cache_arr).any(axis=1)
            logger.info('{0} of {1} cached gamma bins contain NaN.'.format(
                int(bad_bins.sum()), len(bad_bins)))
            if bad_bins.any():
                logger.warning('Gammas with NaN cache entries: {0}'.format(
                    np.asarray(spectra.gammas)[bad_bins]))
        except Exception as exc:
            logger.warning('Could not check cache for NaN: {0}'.format(exc))

        logger.info('Fitting DFE.')
        # sel_params = [0.2, 10.]
        initial_guesses = []
        initial_guesses.append([initial_alpha, initial_beta])
        initial_guesses.append([0.2, 1.])
        initial_guesses.append([0.2, 10.])
        initial_guesses.append([0.2, 100.])
        initial_guesses.append([0.2, 1000.])
        initial_guesses.append([0.2, 10000])
        initial_guesses.append([0.2, 100000])
        initial_guesses.append([0.2, 1000000])
        initial_guesses.append([1., 0.1])
        initial_guesses.append([1., 1.])
        initial_guesses.append([1., 10.])
        initial_guesses.append([1., 100.])
        initial_guesses.append([1., 1000.])
        initial_guesses.append([1., 10000])
        initial_guesses.append([1., 100000])
        initial_guesses.append([1., 1000000])
        initial_guesses.append([5., 0.1])
        initial_guesses.append([5., 1.])
        initial_guesses.append([5., 10.])
        initial_guesses.append([5., 100.])
        initial_guesses.append([5., 1000.])
        initial_guesses.append([5., 10000.])
        initial_guesses.append([5., 100000.])
        initial_guesses.append([5., 1000000.])
        initial_guesses.append([10., 100000.])

        max_ll = -np.inf
        best_model = None
        best_params = None

        for i in range(25):
            sel_params = initial_guesses[i]
            # Perturb each start (previously computed but never used) and
            # keep it inside the bounds.
            p0 = dadi.Misc.perturb_params(
                sel_params, lower_bound=lower_bound,
                upper_bound=upper_bound)
            # Optimization itself runs on the unfolded cache -- spectra.integrate
            # is called internally by optimize_log_fmin on every iteration, and
            # is left untouched here regardless of `fold`.
            popt = dadi.Inference.optimize_log_fmin(p0, nonsyn_data,
                spectra.integrate, pts=None,
                func_args=[DFE.PDFs.gamma, theta_nonsyn],
                lower_bound=lower_bound, upper_bound=upper_bound,
                verbose=len(sel_params), maxiter=50,
                multinom=use_multinom)
            model_sfs = spectra.integrate(
                popt, None, DFE.PDFs.gamma, theta_nonsyn, None)
            # Fold here, once per start, before scaling/likelihood -- this is
            # outside optimize_log_fmin's internal objective evaluations, not
            # inside them.
            if fold:
                if not model_sfs.folded:
                    model_sfs = model_sfs.fold()
            if use_multinom:
                model_sfs = dadi.Inference.optimally_scaled_sfs(
                    model_sfs, nonsyn_data)
                this_ll = dadi.Inference.ll_multinom(model_sfs, nonsyn_data)
            else:
                # Poisson: model SFS stays at theta_nonsyn scale.
                this_ll = dadi.Inference.ll(model_sfs, nonsyn_data)
            print(this_ll)
            logger.info('Gamma DFE start {0}: params = {1}, ll = {2}'.format(
                i, popt, this_ll))
            if this_ll > max_ll:
                best_model = model_sfs
                max_ll = this_ll
                best_params = popt

        if best_model is None:
            logger.error('All gamma DFE starts returned NaN likelihoods. '
                         'Check the cache NaN diagnostic above.')
            sys.exit(1)

        logger.info('best_model folded: {0}'.format(best_model.folded))

        # ------------------------------------------------------------------
        # NEUGAMMA: neutral + gamma DFE inference (commented out)
        # ------------------------------------------------------------------
        # logger.info('Fitting Neutral + Gamma DFE')
        #
        # ng_initial_guesses = []
        # ng_initial_guesses.append([0.15, 0.2, 0.1])
        # ng_initial_guesses.append([0.15, 0.2, 1.])
        # ng_initial_guesses.append([0.15, 0.2, 10.])
        # ng_initial_guesses.append([0.15, 0.2, 100.])
        # ng_initial_guesses.append([0.15, 0.2, 1000.])
        # ng_initial_guesses.append([0.15, 0.2, 10000])
        # ng_initial_guesses.append([0.15, 0.2, 100000])
        # ng_initial_guesses.append([0.15, 0.2, 1000000])
        # ng_initial_guesses.append([0.15, 1., 0.1])
        # ng_initial_guesses.append([0.15, 1., 1.])
        # ng_initial_guesses.append([0.15, 1., 10.])
        # ng_initial_guesses.append([0.15, 1., 100.])
        # ng_initial_guesses.append([0.15, 1., 1000.])
        # ng_initial_guesses.append([0.15, 1., 10000])
        # ng_initial_guesses.append([0.15, 1., 100000])
        # ng_initial_guesses.append([0.15, 1., 1000000])
        # ng_initial_guesses.append([0.15, 5., 0.1])
        # ng_initial_guesses.append([0.15, 5., 1.])
        # ng_initial_guesses.append([0.15, 5., 10.])
        # ng_initial_guesses.append([0.15, 5., 100.])
        # ng_initial_guesses.append([0.15, 5., 1000.])
        # ng_initial_guesses.append([0.15, 5., 10000.])
        # ng_initial_guesses.append([0.15, 5., 100000.])
        # ng_initial_guesses.append([0.15, 5., 1000000.])
        # ng_initial_guesses.append([0.15, 10., 100000.])
        #
        # ng_max_ll = -np.inf
        # ng_best_model = None
        # ng_best_params = None
        # for i in range(25):
        #     sel_params = ng_initial_guesses[i]
        #     lower_bound, upper_bound = [1e-15, 1e-15, 1e-2], [1, 1000, 1000000.]
        #     ng_p0 = dadi.Misc.perturb_params(sel_params, lower_bound=lower_bound,
        #                                   upper_bound=upper_bound)
        #     ng_popt = dadi.Inference.optimize_log_fmin(ng_p0, nonsyn_data,
        #                                           spectra.integrate, pts=None,
        #                                           func_args=[self.neugamma,
        #                                                      theta_nonsyn],
        #                                           lower_bound=lower_bound,
        #                                           upper_bound=upper_bound,
        #                                           verbose=len(sel_params),
        #                                           maxiter=50, multinom=False)
        #     ng_model_sfs = spectra.integrate(
        #         ng_popt, None, self.neugamma, theta_nonsyn, None)
        #     if fold:
        #         if not ng_model_sfs.folded:
        #             ng_model_sfs = ng_model_sfs.fold()
        #     ng_this_ll = dadi.Inference.ll(ng_model_sfs, nonsyn_data)
        #     logger.info('Neutral+gamma DFE start {0}: params = {1}, '
        #                 'll = {2}'.format(i, ng_popt, ng_this_ll))
        #     if ng_this_ll > ng_max_ll:
        #         ng_best_model = ng_model_sfs
        #         ng_max_ll = ng_this_ll
        #         ng_best_params = ng_popt
        #
        # if ng_best_model is None:
        #     logger.error('All neutral+gamma DFE starts returned NaN '
        #                  'likelihoods.')
        #     sys.exit(1)

        logger.info('Computing intermediate summary statistics.')

        mu = 1.5E-8
        Na = theta_syn / (4 * L_syn * mu) # Estimate of N_anc
        #  Theta is given in terms of 4 * L_syn * mu * Na, thus,
        #  we don't need to divide by nu in order to get ancestral popultion
        #

        # Save the best-fit model SFS to disk, in dadi's native Spectrum
        # format. best_model already reflects the requested fold state
        # (folded once per start, above), so no additional fold is applied
        # here -- dadi raises an error if you try to fold a spectrum that
        # is already folded.
        gamma_model_sfs_file = \
            '{0}{1}gamma_dfe_best_model_sfs.fs'.format(
                args['outprefix'], underscore)
        # NEUGAMMA:
        # neugamma_model_sfs_file = \
        #     '{0}{1}neugamma_dfe_best_model_sfs.fs'.format(
        #         args['outprefix'], underscore)

        best_model.to_file(gamma_model_sfs_file)
        # NEUGAMMA:
        # ng_best_model.to_file(neugamma_model_sfs_file)

        logger.info('Saved best-fit gamma-DFE model SFS to {0}\n'.format(
            gamma_model_sfs_file))
        # NEUGAMMA:
        # logger.info(
        #     'Saved best-fit neutral+gamma-DFE model SFS to {0}\n'.format(
        #         neugamma_model_sfs_file))

        # Unscaled (i.e., population-size-scaled) parameter estimates,
        # computed once so they can be reused for both the text and CSV
        # outputs.
        gamma_unscaled_params = np.divide(
            best_params, np.array([1, 2 * Na]))
        # NEUGAMMA:
        # neugamma_unscaled_params = np.divide(
        #     ng_best_params, np.array([1, 1, 2 * Na]))

        logger.info('Outputting results.')

        with open(inferred_DFE, 'w') as f:
            f.write('Assuming a gamma-distributed DFE...\n')
            f.write('The maximum likelihood is: ' +
                    str(max_ll) +  '.\n')
            f.write('The maximum likelihood scaled ' +
                    'gamma-distributed DFE is parameterized ' +
                    'as: ' + str(best_params)  + '.\n')
            # Non-scaled params are divided by
            f.write(
                'The maximum likelihood parameters have an '
                'estimate of: '
                '{0}.\n'.format(gamma_unscaled_params))
            f.write('The empirical nonsynonymous SFS is: ' +
                    str(nonsyn_data) + '.\n')
            f.write('The best fit model SFS is: ' +
                    str(best_model) + '.\n')
            # NEUGAMMA:
            # f.write('Assuming a neu-gamma-distributed DFE...\n')
            # f.write('The maximum likelihood is: ' +
            #         str(ng_max_ll) +  '.\n')
            # f.write('The maximum likelihood scaled ' +
            #         'neu-gamma-distributed DFE is parameterized ' +
            #         'as: ' + str(ng_best_params)  + '.\n')
            # # Non-scaled params are divided by
            # f.write(
            #     'The maximum likelihood parameters have an '
            #     'estimate of: '
            #     '{0}.\n'.format(neugamma_unscaled_params))
            # f.write('The empirical nonsynonymous SFS is: ' +
            #         str(nonsyn_data) + '.\n')
            # f.write('The best fit model SFS is: ' +
            #         str(ng_best_model) + '.\n')

        # Write a CSV summary with one row per DFE model fitted
        # (currently the gamma-distributed DFE only).
        csv_fieldnames = [
            'dfe_model', 'max_log_likelihood', 'scaled_params',
            'unscaled_params', 'theta_syn', 'theta_nonsyn', 'L_syn',
            'ancestral_Na_estimate', 'fold', 'best_model_sfs_file',
            'likelihood', 'ns_s_ratio', 'pts', 'lower_bound', 'upper_bound',
        ]
        with open(dfe_results_csv, 'w', newline='') as csv_handle:
            csv_writer = csv.DictWriter(csv_handle, fieldnames=csv_fieldnames)
            csv_writer.writeheader()
            csv_writer.writerow({
                'dfe_model': 'gamma',
                'max_log_likelihood': max_ll,
                'scaled_params': json.dumps(
                    list(np.atleast_1d(best_params).astype(float))),
                'unscaled_params': json.dumps(
                    list(np.atleast_1d(
                        gamma_unscaled_params).astype(float))),
                'theta_syn': theta_syn,
                'theta_nonsyn': theta_nonsyn,
                'L_syn': L_syn,
                'ancestral_Na_estimate': Na,
                'fold': fold,
                'best_model_sfs_file': gamma_model_sfs_file,
                'likelihood': likelihood,
                'ns_s_ratio': ns_s_ratio,
                'pts': json.dumps(pts_l),
                'lower_bound': json.dumps(lower_bound),
                'upper_bound': json.dumps(upper_bound),
            })
            # NEUGAMMA:
            # csv_writer.writerow({
            #     'dfe_model': 'neutral_plus_gamma',
            #     'max_log_likelihood': ng_max_ll,
            #     'scaled_params': json.dumps(
            #         list(np.atleast_1d(ng_best_params).astype(float))),
            #     'unscaled_params': json.dumps(
            #         list(np.atleast_1d(
            #             neugamma_unscaled_params).astype(float))),
            #     'theta_syn': theta_syn,
            #     'theta_nonsyn': theta_nonsyn,
            #     'L_syn': L_syn,
            #     'ancestral_Na_estimate': Na,
            #     'fold': fold,
            #     'best_model_sfs_file': neugamma_model_sfs_file,
            # })
        logger.info('Wrote CSV summary of DFE fits to {0}\n'.format(
            dfe_results_csv))
        logger.info('Pipeline executed succesfully.')

if __name__ == '__main__':
    DFEInference().main()
