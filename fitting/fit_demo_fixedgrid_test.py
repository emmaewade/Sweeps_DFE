""" Uses dadi to infer the best-fit demographic model for given sfs.

JCM 20231212

Modified to additionally:
    - Write a CSV summary of results across all fitted models.
    - Save the empirical synonymous SFS and the best-fit model SFS
      (both non-scaled and scaled by theta) to disk in dadi's native
      Spectrum format (``.fs``), so they can be reloaded downstream
      (e.g., by ``fitdadi_infer_DFE.py``) without recomputation.
    - Add a --fold/--no-fold flag. When folding, the input SFS is folded
      once at load time, and each candidate model spectrum is folded
      before its likelihood is computed against the (matching) folded
      empirical data.
    - Comment out contraction (nu < 1) starting points for the two_epoch
      and three_epoch models, leaving only stasis/growth starts.
    - Size each model's optimization loop to the length of its own
      initial_guesses list, rather than a hardcoded number of starts.

[202609 fixed-grid test] Copy of 202608/sweeps/fit_demo_jon_output_improved.py
with a --pts option (default 1200 1600 2000) replacing the n*2/n*3/n*4 grid,
and --nu_bounds / --T_bounds (defaults [0.01, open] / [1e-5, 10]) on the
two_epoch and three_epoch optimizations (previously unbounded).
"""


import sys
import os
import csv
import json
import logging
import time
import argparse
import warnings

import numpy
import dadi
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


class DemographicInference():
    """Wrapper class to allow functions to reference each other."""

    def ExistingFile(self, fname):
        """Return *fname* if existing file, otherwise raise ValueError."""
        if os.path.isfile(fname):
            return fname
        else:
            raise ValueError("%s must specify a valid file name" % fname)

    def fitDemographicModelParser(self):
        """Return *argparse.ArgumentParser* for ``fitdadi_infer_DFE.py``."""
        parser = ArgumentParserNoArgHelp(
            description=(
                'Given the number of individuals in population one and two, '
                'this script outputs a `*pops_file.txt` in the format '
                'specified for use by the python package, `easySFS.py`.'),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument(
            'syn_input_sfs', type=self.ExistingFile,
            help=('Synonynomous site-frequency spectrum from which the '
                  'demographic parameters should be inferred.'))
        parser.add_argument(
            '--mask_singletons', dest='mask_singletons',
            help=('Boolean flag for masking singlestons in Spectrum.'),
            action='store_true')
        parser.add_argument(
            '--mask_doubletons', dest='mask_doubletons',
            help=('Boolean flag for masking doublestons in Spectrum.'),
            action='store_true')
        parser.add_argument(
            '--sum_alleles', dest='sum_alleles',
            help=('Boolean flag to determine if genome length should be '
                  'determined by calculated sum of allele count.'),
            action='store_true')
        parser.set_defaults(mask_singletons=False)
        parser.set_defaults(mask_doubletons=False)
        parser.set_defaults(sum_alleles=False)
        parser.add_argument(
            '--L_syn', type=float,
            dest='L_syn',
            help=('Sum of allele counts for synonymous SFS.'),
            default=100000000.0)
        parser.add_argument(
            '--model_type', type=str,
            dest='model_type',
            help=('Type of demographic model.'))
        parser.add_argument(
            '--input_params', type=float,
            dest='input_params', nargs='*',
            help=('Input starting parameters for MLE search.'),
            default=[])
        parser.add_argument(
            '--maxiter', type=int, default=25,
            help=('Maximum L-BFGS-B iterations per start. 25 (the default) '
                  'can stop three_epoch fits short of the optimum; the '
                  'refit of non-converged three_epoch fits '
                  '(refit_three_epoch.sh) uses 100.'))
        parser.add_argument(
            '--fold', dest='fold', action='store_true',
            help=('Fold the input SFS and each candidate model SFS.'))
        parser.add_argument(
            '--no-fold', dest='fold', action='store_false',
            help=('Use unfolded SFS throughout.'))
        parser.set_defaults(fold=True)
        parser.add_argument(
            '--nu_bounds', type=float, nargs=2, default=[1e-2, float('inf')],
            help=('Lower/upper bound on every relative population size '
                  '(nu) in two_epoch / three_epoch. inf = unbounded. The '
                  'upper bound is left open on purpose: if every parameter '
                  'has both bounds, L-BFGS-B treats the problem as "boxed" '
                  'and takes an unscaled first step along the gradient, which '
                  'throws fits to the box edges.'))
        parser.add_argument(
            '--T_bounds', type=float, nargs=2, default=[1e-5, 10.],
            help=('Lower/upper bound on every epoch duration (T, in 2*Na '
                  'generations) in two_epoch / three_epoch. Integration time '
                  'grows with T, and unbounded three_epoch fits ran T off to '
                  '1e4-1e6, taking hours per likelihood evaluation.'))
        parser.add_argument(
            '--pts', type=int, nargs=3, default=[1200, 1600, 2000],
            help=('dadi grid sizes (previously n*2, n*3, n*4, which is only '
                  '20/30/40 points for n = 10).'))
        parser.add_argument(
            'outprefix', type=str,
            help='The file prefix for the output files')
        return parser

    def compute_allele_sum(self, file_path):
        """Return sum of allele counts from input SFS file.
        """
        with open(file_path, 'r') as file:
             lines = file.readlines()
             if len(lines) >= 2:
                 second_line = lines[1].strip()
                 values = list(map(float, second_line.split()))
                 array_sum = sum(values)
                 return array_sum

    def snm(self, notused, ns, pts):
        """Return a standard neutral model.

        ns = (n1, )
            n1: Number of samples in resulting Spectrum object.
        pts: Number of grid points to use in integration.
        """
        xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
        phi = dadi.PhiManip.phi_1D(xx)  # Define initial phi.

        # Construct Spectrum object.
        fs = dadi.Spectrum.from_phi(phi, ns, (xx, ))
        return fs

    def two_epoch(self, params, ns, pts):
        """Define a two-epoch demography, i.e., an instantaneous size change.

        params = (nu, T)
            nu: Ratio of contemporary to ancient population size.
            T: Time in the past at which size change occured,
                in units of 2*N_a.
        ns = (n1, )
            n1: Number of samples in resulting Spectrum object.
        pts: Number of grid points to use in integration.
        """
        nu, T = params  # Define given parameters.
        xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
        phi = dadi.PhiManip.phi_1D(xx)  # Define initial phi.

        phi = dadi.Integration.one_pop(phi, xx, T, nu)  # Integrate.

        # Construct Spectrum object.
        fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
        return fs

    def two_epoch_sel(self, params, ns, pts):
        """Define a two-epoch demography, i.e., an instantaneous size change.

        This method incorporates a gamma parameter.

        params = (nu, T, gamma)
            nu: Ratio of contemporary to ancient population size.
            T: Time in the past at which size change occured,
                in units of 2*N_a.
            gamma: Parameter tuple describing a gamma distribution.
        ns = (n1, )
            n1: Number of samples in resulting Spectrum object.
        pts: Number of grid points to use in integration.
        """
        nu, T, gamma = params  # Define given parameters.
        xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
        phi = dadi.PhiManip.phi_1D(xx, gamma=gamma)  # Define initial phi

        phi = dadi.Integration.one_pop(phi, xx, T, nu, gamma=gamma)

        # Construct Spectrum object.
        fs = dadi.Spectrum.from_phi(phi, ns, (xx, ))
        return fs

    def growth(self, params, ns, pts):
        """Exponential growth beginning some time ago.

        params = (nu, T)
            nu: Ratio of contemporary to ancient population size.
            T: Time in the past at which size change occured,
                in units of 2*N_a.
        ns = (n1, )
            n1: Number of samples in resulting Spectrum object.
        pts: Number of grid points to use in integration.
        """
        nu, T = params  # Define given parameters.
        xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
        phi = dadi.PhiManip.phi_1D(xx)  # Define initial phi.

        def nu_func(t): return numpy.exp(numpy.log(nu) * t / T)  # Exp growth.
        phi = dadi.Integration.one_pop(phi, xx, T, nu_func)  # Integrate.

        # Construct Spectrum object.
        fs = dadi.Spectrum.from_phi(phi, ns, (xx, ))
        return fs

    def bottlegrowth(self, params, ns, pts):
        """Instantaneous size change followed by exponential growth.

        params = (nuB, nuF, T)
            nuB: Ratio of population size after instantaneous change to ancient
                population size.
            nuF: Ratio of contemporary to ancient population size.
            T: Time in the past at which size change occured,
                in units of 2*N_a.
        ns = (n1, )
            n1: Number of samples in resulting Spectrum object.
        pts: Number of grid points to use in integration.
        """
        nuB, nuF, T = params  # Define given parameters.

        xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
        phi = dadi.PhiManip.phi_1D(xx)  # Define initial phi

        # Exponential growth function
        def nu_func(t): return nuB * numpy.exp(numpy.log(nuF / nuB) * t / T)

        phi = dadi.Integration.one_pop(phi, xx, T, nu_func)  # Integrate.

        # Construct Spectrum object.
        fs = dadi.Spectrum.from_phi(phi, ns, (xx, ))
        return fs

    def three_epoch(self, params, ns, pts):
        """Define a three-epoch demography.

        params = (nuB, nuF, TB, TF)
            nuB: Ratio of bottleneck population size to ancient
                population size.
            nuF: Ratio of contemporary to ancient population size.
            TB: Length of bottleneck, in units of 2 * N_a.
            TF: Time since bottleneck recovery, in units of 2 * N_a.
        ns = (n1, )
            n1: Number of samples in resulting Spectrum object.
        pts: Number of grid points to use in integration.
        """
        nuB, nuF, TB, TF = params  # Define given parameters.

        xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
        phi = dadi.PhiManip.phi_1D(xx)  # Define initial phi.

        phi = dadi.Integration.one_pop(phi, xx, TB, nuB)  # Integrate 1 to 2.
        phi = dadi.Integration.one_pop(phi, xx, TF, nuF)  # Integrate 2 to 3.

        fs = dadi.Spectrum.from_phi(phi, ns, (xx, ))
        return fs

    def four_epoch(self, params, ns, pts):
        """Define a four-epoch demography.

        params = (Na, Nb, Nc, Ta, Tb, Tc)
            Na: ratio of population size between epoch 1 and 2.
            Nb: ratio of population size between epoch 2 and 3.
            Nc: ratio of population size between epoch 3 and 4.
            Ta: Bottleneck length between epoch 1 and 2, in units of 2 * N_a.
            Tb: Length of bottleneck between epoch 2 and 3,
                in units of 2 * N_a.
            Tc: Length of bottleneck between epoch 3 and 4,
                in units of 2 * N_a.
        ns = (n1, )
            n1: Number of samples in resulting Spectrum object.
        pts: Number of grid points to use in integration.
        """
        Na, Nb, Nc, Ta, Tb, Tc = params  # Define given parameters.

        xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
        phi = dadi.PhiManip.phi_1D(xx)  # Define initial phi.

        # Integrate epochs.
        phi = dadi.Integration.one_pop(phi, xx, Ta, Na)  # Integrate 1 to 2.
        phi = dadi.Integration.one_pop(phi, xx, Tb, Nb)  # Integrate 2 to 3.
        phi = dadi.Integration.one_pop(phi, xx, Tc, Nc)  # Integrate 3 to 4.

        # Construct spectrum object.
        fs = dadi.Spectrum.from_phi(phi, ns, (xx, ))
        return fs

    def one_epoch(self, params, ns, pts):
        """Define a one-epoch demography.
        params = (nu, T)
            nu: Ratio of contemporary to ancient population size.
            T: Time in the past at which size change occured,
                in units of 2*N_a.
        ns = (n1, )
            n1: Number of samples in resulting Spectrum object.
        pts: Number of grid points to use in integration.
        """
        xx = dadi.Numerics.default_grid(pts)  # Define likelihood surface.
        phi = dadi.PhiManip.phi_1D(xx)  # Define initial phi.

        # Construct spectrum object
        fs = dadi.Spectrum.from_phi(phi, ns, (xx, ))
        return fs


    def one_pop(phi, xx, T, nu=1, gamma=0, h=0.5, theta0=1.0, initial_t=0,
                frozen=False, beta=1):
        """
        Integrate a 1-dimensional phi foward.

        phi: Initial 1-dimensional phi
        xx: Grid upon (0,1) overwhich phi is defined.

        nu, gamma, and theta0 may be functions of time.
        nu: Population size
        gamma: Selection coefficient on *all* segregating alleles
        h: Dominance coefficient. h = 0.5 corresponds to genic selection.
        Heterozygotes have fitness 1+2sh and homozygotes have fitness 1+2s.
        theta0: Propotional to ancestral size. Typically constant.
        beta: Breeding ratio, beta=Nf/Nm.

        T: Time at which to halt integration
        initial_t: Time at which to start integration. (Note that this only
            matters if one of the demographic parameters is a function of time.)

        frozen: If True, population is 'frozen' so that it does not change.
            In the one_pop case, this is equivalent to not running the
            integration at all.
        """
        phi = phi.copy()

        # For a one population integration, freezing means just not integrating.
        if frozen:
            return phi

        if T - initial_t == 0:
            return phi
        elif T - initial_t < 0:
            raise ValueError('Final integration time T (%f) is less than '
                            'intial_time (%f). Integration cannot be run '
                            'backwards.' % (T, initial_t))

        vars_to_check = (nu, gamma, h, theta0, beta)
        if numpy.all([numpy.isscalar(var) for var in vars_to_check]):
            return _one_pop_const_params(phi, xx, T, nu, gamma, h, theta0,
                                     initial_t, beta)

        nu_f = Misc.ensure_1arg_func(nu)
        gamma_f = Misc.ensure_1arg_func(gamma)
        h_f = Misc.ensure_1arg_func(h)
        theta0_f = Misc.ensure_1arg_func(theta0)
        beta_f = Misc.ensure_1arg_func(beta)

        current_t = initial_t
        nu, gamma, h = nu_f(current_t), gamma_f(current_t), h_f(current_t)
        beta = beta_f(current_t)
        dx = numpy.diff(xx)
        while current_t < T:
            dt = _compute_dt(dx,nu,[0],gamma,h)
            this_dt = min(dt, T - current_t)

            # Because this is an implicit method, I need the *next* time's params.
            # So there's a little inconsistency here, in that I'm estimating dt
            # using the last timepoints nu,gamma,h.
            next_t = current_t + this_dt
            nu, gamma, h = nu_f(next_t), gamma_f(next_t), h_f(next_t)
            beta = beta_f(next_t)
            theta0 = theta0_f(next_t)

            if numpy.any(numpy.less([T,nu,theta0], 0)):
                raise ValueError('A time, population size, migration rate, or '
                             'theta0 is < 0. Has the model been mis-specified?')
            if numpy.any(numpy.equal([nu], 0)):
                raise ValueError('A population size is 0. Has the model been '
                             'mis-specified?')

            _inject_mutations_1D(phi, this_dt, xx, theta0)
            # Do each step in C, since it will be faster to compute the a,b,c
            # matrices there.
            phi = int_c.implicit_1Dx(phi, xx, nu, gamma, h, beta, this_dt,
                                 use_delj_trick=use_delj_trick)
            current_t = next_t
        return phi

    def main(self):
        """Execute main function."""
        # Parse command line arguments
        parser = self.fitDemographicModelParser()
        args = vars(parser.parse_args())
        prog = parser.prog

        # Assign arguments
        syn_input_sfs = args['syn_input_sfs']
        outprefix = args['outprefix']
        mask_singletons = args['mask_singletons']
        mask_doubletons = args['mask_doubletons']
        sum_alleles = args['sum_alleles']
        L_syn = args['L_syn']
        input_params = args['input_params']
        model_type = args['model_type']
        fold = args['fold']
        pts_arg = list(args['pts'])
        nu_lo, nu_hi = args['nu_bounds']
        T_lo, T_hi = args['T_bounds']

        # Numpy options
        numpy.set_printoptions(linewidth=numpy.inf)

        # create output directory if needed
        outdir = os.path.dirname(args['outprefix'])
        if outdir:
            if not os.path.isdir(outdir):
                if os.path.isfile(outdir):
                    os.remove(outdir)
                os.mkdir(outdir)

        # Output files: logfile
        # Remove output files if they already exist
        underscore = '' if args['outprefix'][-1] == '/' else '_'
        exponential_growth_demography = \
            '{0}{1}exponential_growth_demography.txt'.format(
                args['outprefix'], underscore)
        two_epoch_demography = \
            '{0}{1}two_epoch_demography.txt'.format(
                args['outprefix'], underscore)
        bottleneck_growth_demography = \
            '{0}{1}bottleneck_growth_demography.txt'.format(
                args['outprefix'], underscore)
        three_epoch_demography = \
            '{0}{1}three_epoch_demography.txt'.format(
                args['outprefix'], underscore)
        one_epoch_demography = \
           '{0}{1}one_epoch_demography.txt'.format(
                args['outprefix'], underscore)
        expected_sfs_comparison = \
           '{0}{1}expected_sfs_comparison.txt'.format(
                args['outprefix'], underscore)
        logfile = '{0}{1}log.log'.format(args['outprefix'], underscore)
        # CSV summary output, collecting one row per fitted model.
        demography_results_csv = \
            '{0}{1}demography_results.csv'.format(
                args['outprefix'], underscore)
        # Directory/prefix for saved SFS (.fs) files.
        empirical_syn_sfs_file = \
            '{0}{1}empirical_syn_sfs.fs'.format(
                args['outprefix'], underscore)
        to_remove = [logfile, demography_results_csv]
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

        # Construct initial Spectrum object from input synonymous sfs.
        syn_data = dadi.Spectrum.from_file(syn_input_sfs)
        if fold:
            if not syn_data.folded:
                syn_data = syn_data.fold()
        if mask_singletons:
            syn_data.mask[1] = True
        if mask_doubletons:
            syn_data.mask[2] = True
        if sum_alleles:
            L_syn = self.compute_allele_sum(syn_input_sfs)
            logger.info('Calculating allele sum based on number of variants')

        # Save a copy of the (possibly masked) empirical synonymous SFS
        # that was actually used for inference, in dadi's native format.
        syn_data.to_file(empirical_syn_sfs_file)
        logger.info('Saved empirical synonymous SFS to {0}\n'.format(
            empirical_syn_sfs_file))

        syn_ns = syn_data.sample_sizes  # Number of samples.
        pts_l = pts_arg
        if min(pts_l) <= syn_ns[0]:
            logger.error('Grid pts {0} must exceed the sample size {1}.'.format(
                pts_l, syn_ns[0]))
            sys.exit(1)
        logger.info('dadi grid pts: {0}\n'.format(pts_l))

        if model_type:
            model_list = [model_type]
        else:
            model_list = ['one_epoch', 'two_epoch', 'three_epoch']

        # CSV summary file: header, one row appended per fitted model below.
        csv_fieldnames = [
            'model', 'best_params', 'multinomial_log_likelihood',
            'poisson_log_likelihood', 'theta_syn', 'theta_nonsyn',
            'L_syn', 'ancestral_Ne_estimate', 'extra_estimates', 'fold',
            'non_scaled_sfs_file', 'scaled_sfs_file',
        ]
        csv_file_exists = os.path.isfile(demography_results_csv)
        csv_handle = open(demography_results_csv, 'a', newline='')
        csv_writer = csv.DictWriter(csv_handle, fieldnames=csv_fieldnames)
        if not csv_file_exists:
            csv_writer.writeheader()

        # Fit different epoch models and compute likelihood
        for model in model_list:
            if model == 'exponential_growth':
                initial_guesses = []
                initial_guesses.append([0.1, 1])
                initial_guesses.append([0.2, 1])
                initial_guesses.append([0.3, 1])
                initial_guesses.append([0.4, 1])
                initial_guesses.append([1, 1])
                initial_guesses.append([6, 1])
                initial_guesses.append([7, 1])
                initial_guesses.append([8, 1])
                initial_guesses.append([9, 1])
                initial_guesses.append([10, 1])
                initial_guesses.append([12.5, 1])
                initial_guesses.append([15, 1])
                initial_guesses.append([17.5, 1])
                initial_guesses.append([20, 1])
                initial_guesses.append([22.5, 1])
                initial_guesses.append([25, 1])
                initial_guesses.append([20.75, 1])
                initial_guesses.append([30, 1])
                initial_guesses.append([33.3, 1])
                initial_guesses.append([36.6, 1])
                initial_guesses.append([40, 1])
                initial_guesses.append([50, 1])
                initial_guesses.append([60, 1])
                initial_guesses.append([70, 1])
                initial_guesses.append([80, 1])
                file = exponential_growth_demography
                func_ex = dadi.Numerics.make_extrap_log_func(self.growth)
                logger.info('Beginning demographic inference for exponential '
                            'growth demographic model.')
            elif model == 'two_epoch':
                initial_guesses = []
                # Contraction starting points (nu < 1) commented out --
                # only stasis/growth starts remain.
                # initial_guesses.append([0.05, 0.1])
                # initial_guesses.append([0.075, 0.1])
                # initial_guesses.append([0.1, 0.1])
                # initial_guesses.append([0.25, 0.1])
                # initial_guesses.append([0.5, 0.1])
                # initial_guesses.append([0.75, 0.1])
                initial_guesses.append([1, 0.1])
                initial_guesses.append([1.25, 0.1])
                initial_guesses.append([1.5, 0.1])
                initial_guesses.append([1.75, 0.1])
                initial_guesses.append([2, 0.1])
                initial_guesses.append([2.25, 0.1])
                initial_guesses.append([2.5, 0.1])
                initial_guesses.append([2.75, 0.1])
                initial_guesses.append([3, 0.1])
                initial_guesses.append([3.25, 0.1])
                initial_guesses.append([3.5, 0.1])
                initial_guesses.append([3.75, 0.1])
                initial_guesses.append([4, 0.1])
                initial_guesses.append([4.25, 0.1])
                initial_guesses.append([4.5, 0.1])
                initial_guesses.append([4.75, 0.1])
                initial_guesses.append([5, 0.1])
                # More extreme growth starting points.
                initial_guesses.append([6, 0.1])
                initial_guesses.append([7, 0.1])
                initial_guesses.append([8, 0.1])
                initial_guesses.append([9, 0.1])
                initial_guesses.append([10, 0.1])
                initial_guesses.append([12.5, 0.1])
                initial_guesses.append([15, 0.1])
                initial_guesses.append([17.5, 0.1])
                initial_guesses.append([20, 0.1])
                initial_guesses.append([25, 0.1])
                initial_guesses.append([30, 0.1])
                initial_guesses.append([40, 0.1])
                initial_guesses.append([50, 0.1])
                initial_guesses.append([60, 0.1])
                initial_guesses.append([70, 0.1])
                initial_guesses.append([80, 0.1])
                file = two_epoch_demography
                func_ex = dadi.Numerics.make_extrap_log_func(self.two_epoch)
                logger.info('Beginning demographic inference for two-epoch '
                            'demographic model.')
            elif model == 'bottleneck_growth':
                initial_guess = [0.1, 0.1, 0.1]
                initial_guesses = []
                initial_guesses.append([0.1, 0.1, 0.1])
                initial_guesses.append([0.2, 0.2, 0.1])
                initial_guesses.append([0.3, 0.3, 0.1])
                initial_guesses.append([0.4, 0.4, 0.1])
                initial_guesses.append([1, 1, 0.1])
                initial_guesses.append([6, 6, 0.1])
                initial_guesses.append([7, 7, 0.1])
                initial_guesses.append([8, 8, 0.1])
                initial_guesses.append([9, 9, 0.1])
                initial_guesses.append([10, 10, 0.1])
                initial_guesses.append([12.5, 12.5, 0.1])
                initial_guesses.append([15, 15, 0.1])
                initial_guesses.append([17.5, 17.5, 0.1])
                initial_guesses.append([20, 20, 0.1])
                initial_guesses.append([22.5, 22.5, 0.1])
                initial_guesses.append([25, 25, 0.1])
                initial_guesses.append([27.5, 27.5, 0.1])
                initial_guesses.append([30, 30, 0.1])
                initial_guesses.append([33.3, 33.3, 0.1])
                initial_guesses.append([36.6, 36.6, 0.1])
                initial_guesses.append([40, 40, 0.1])
                initial_guesses.append([50, 50, 0.1])
                initial_guesses.append([60, 60, 0.1])
                initial_guesses.append([70, 70, 0.1])
                initial_guesses.append([80, 80, 0.1])
                file = bottleneck_growth_demography
                func_ex = dadi.Numerics.make_extrap_log_func(self.bottlegrowth)
                logger.info('Beginning demographic inference for bottleneck + '
                            'growth demographic model.')
            elif model == 'three_epoch':
                initial_guesses = []
                # Contraction starting points (nu < 1) commented out --
                # only stasis/growth starts remain.
                # initial_guesses.append([0.05, 0.05, 0.1, 0.1])
                # initial_guesses.append([0.075, 0.075, 0.1, 0.1])
                # initial_guesses.append([0.1, 0.1, 0.1, 0.1])
                # initial_guesses.append([0.25, 0.25, 0.1, 0.1])
                # initial_guesses.append([0.5, 0.5, 0.51, 0.1])
                # initial_guesses.append([0.75, 0.75, 0.1, 0.1])
                initial_guesses.append([1, 1, 0.1, 0.1])
                initial_guesses.append([1.25, 1.25, 0.1, 0.1])
                initial_guesses.append([1.5, 1.5, 0.1, 0.1])
                initial_guesses.append([1.75, 1.75, 0.1, 0.1])
                initial_guesses.append([2, 2, 0.1, 0.1])
                initial_guesses.append([2.25, 2.25, 0.1, 0.1])
                initial_guesses.append([2.5, 2.5, 0.1, 0.1])
                initial_guesses.append([2.75, 2.75, 0.1, 0.1])
                initial_guesses.append([3, 3, 0.1, 0.1])
                initial_guesses.append([3.25, 3.25, 0.1, 0.1])
                initial_guesses.append([3.5, 3.5, 0.1, 0.1])
                initial_guesses.append([3.75, 3.75, 0.1, 0.1])
                initial_guesses.append([4, 4, 0.1, 0.1])
                initial_guesses.append([4.25, 4.25, 0.1, 0.1])
                initial_guesses.append([4.5, 4.5, 0.1, 0.1])
                initial_guesses.append([4.75, 4.75, 0.1, 0.1])
                initial_guesses.append([5, 5, 0.1, 0.1])
                # More extreme growth starting points.
                initial_guesses.append([6, 6, 0.1, 0.1])
                initial_guesses.append([7, 7, 0.1, 0.1])
                initial_guesses.append([8, 8, 0.1, 0.1])
                initial_guesses.append([9, 9, 0.1, 0.1])
                initial_guesses.append([10, 10, 0.1, 0.1])
                initial_guesses.append([12.5, 12.5, 0.1, 0.1])
                initial_guesses.append([15, 15, 0.1, 0.1])
                initial_guesses.append([17.5, 17.5, 0.1, 0.1])
                initial_guesses.append([20, 20, 0.1, 0.1])
                initial_guesses.append([25, 25, 0.1, 0.1])
                initial_guesses.append([30, 30, 0.1, 0.1])
                initial_guesses.append([40, 40, 0.1, 0.1])
                initial_guesses.append([50, 50, 0.1, 0.1])
                initial_guesses.append([60, 60, 0.1, 0.1])
                initial_guesses.append([70, 70, 0.1, 0.1])
                initial_guesses.append([80, 80, 0.1, 0.1])
                file = three_epoch_demography
                func_ex = dadi.Numerics.make_extrap_log_func(self.three_epoch)
                logger.info('Beginning demographic inference for three-epoch '
                            'demographic model.')
            elif model == 'one_epoch':
                initial_guesses = []
                initial_guesses.append([0.1])
                initial_guesses.append([0.2])
                initial_guesses.append([0.3])
                initial_guesses.append([0.4])
                initial_guesses.append([0.5])
                initial_guesses.append([0.6])
                initial_guesses.append([0.7])
                initial_guesses.append([0.8])
                initial_guesses.append([0.9])
                initial_guesses.append([1])
                initial_guesses.append([1.25])
                initial_guesses.append([1.5])
                initial_guesses.append([1.75])
                initial_guesses.append([2])
                initial_guesses.append([2.25])
                initial_guesses.append([2.5])
                initial_guesses.append([2.75])
                initial_guesses.append([3])
                initial_guesses.append([3.33])
                initial_guesses.append([3.66])
                initial_guesses.append([4])
                initial_guesses.append([5])
                initial_guesses.append([6])
                initial_guesses.append([7])
                initial_guesses.append([8])
                file = one_epoch_demography
                func_ex = dadi.Numerics.make_extrap_log_func(self.snm)
                logger.info('Beginning demographic inference for one-epoch '
                            'demographic model.')
            # Parameter bounds. nu >= 0.01 (dadi's time step shrinks as
            # nu -> 0) and T in [1e-5, 10] (integration time grows with T).
            # nu keeps an open upper bound (inf) so the problem is not fully
            # "boxed" -- see --nu_bounds.
            if model == 'two_epoch':            # [nu, T]
                lower_bound = [nu_lo, T_lo]
                upper_bound = [nu_hi, T_hi]
            elif model == 'three_epoch':        # [nuB, nuF, TB, TF]
                lower_bound = [nu_lo, nu_lo, T_lo, T_lo]
                upper_bound = [nu_hi, nu_hi, T_hi, T_hi]
            else:
                lower_bound = None
                upper_bound = None
            logger.info('Parameter bounds for {0}: lower {1}, upper {2}'.format(
                model, lower_bound, upper_bound))
            with open(file, 'w') as f:
                max_likelihood = -1e25
                for i in range(len(initial_guesses)):
                    # Start at initial guess
                    if len(input_params) > 0:
                        p0 = input_params
                    else:
                        p0 = initial_guesses[i]
                    # Randomly perturb parameters before optimization,
                    # except the first start when --input_params is given:
                    # that start is used as is (e.g. the two_epoch MLE
                    # embedded in three_epoch), so the fit cannot end below
                    # the nested model.
                    if not (len(input_params) > 0 and i == 0):
                        p0 = dadi.Misc.perturb_params(
                            p0, fold=1, upper_bound=upper_bound,
                            lower_bound=lower_bound)
                    logger.info(
                        'Beginning optimization with guess, {0}.'.format(p0))
                    popt = dadi.Inference.optimize_log_lbfgsb(
                        p0=p0, data=syn_data, model_func=func_ex, pts=pts_l,
                        lower_bound=lower_bound,
                        upper_bound=upper_bound,
                        verbose=len(p0), maxiter=args['maxiter'])
                    logger.info(
                        'Finished optimization with guess, ' + str(p0) + '.')
                    logger.info('Best fit parameters: {0}.'.format(popt))
                    # Calculate the best-fit model allele-frequency spectrum.
                    # Note, this spectrum needs to be multiplied by "theta".
                    non_scaled_spectrum = func_ex(popt, syn_ns, pts_l)
                    # Fold the model spectrum here, before scaling/likelihood,
                    # conditional on `fold` -- matches the fold state of
                    # syn_data set above.
                    if fold:
                        if not non_scaled_spectrum.folded:
                            non_scaled_spectrum = non_scaled_spectrum.fold()
                    # Likelihood of the data given the model AFS.
                    multinomial_ll_non_scaled_spectrum = \
                        dadi.Inference.ll_multinom(
                            model=non_scaled_spectrum, data=syn_data)
                    logger.info(
                        'Maximum log composite likelihood: {0}.'.format(
                            multinomial_ll_non_scaled_spectrum))
                    theta = dadi.Inference.optimal_sfs_scaling(
                        non_scaled_spectrum, syn_data)
                    logger.info(
                        'Optimal value of theta: {0}.'.format(theta))
                    if multinomial_ll_non_scaled_spectrum > max_likelihood:
                        best_params = popt
                        best_non_scaled_spectrum = non_scaled_spectrum
                        max_likelihood = multinomial_ll_non_scaled_spectrum
                        theta_syn = theta
                best_scaled_spectrum = theta_syn * best_non_scaled_spectrum
                theta_nonsyn = theta_syn * 2.14
                poisson_ll = dadi.Inference.ll(
                    model=best_scaled_spectrum, data=syn_data)
                f.write('Best fit parameters: {0}.\n'.format(best_params))
                f.write(
                    'Maximum multinomial log composite '
                    'likelihood: {0}.\n'.format(
                        max_likelihood))
                f.write(
                    'Maximum poisson log composite likelihood: {0}.\n'.format(
                        poisson_ll))
                f.write('Non-scaled best-fit model spectrum: {0}.\n'.format(
                    best_non_scaled_spectrum))
                f.write('Optimal value of theta_syn: {0}.\n'.format(theta_syn))
                f.write('Optimal value of theta_nonsyn: {0}.\n'.format(
                    theta_nonsyn))
                f.write('Empirical synonymous sfs: {0}.\n'.format(
                    syn_data))
                f.write('Scaled best-fit model spectrum: {0}.\n'.format(
                    best_scaled_spectrum))
                logger.info('Converting and interpreting Tau estimates.')
                simulation_mu = 1.5E-8
                # Dict to hold any model-specific derived estimates
                # (e.g., generation times), later serialized into the CSV.
                extra_estimates = {}
                if model == 'one_epoch':
                    ancestral_Ne = theta_syn / (4 * L_syn * simulation_mu)
                    f.write(
                        'Estimate for ancestral population size is ' +
                        str(ancestral_Ne) + '.\n')
                elif model == 'two_epoch':
                    best_nu = best_params[0]
                    best_tau = best_params[1]
                    generations = 2 * best_tau * theta_syn / (4 * simulation_mu * L_syn)
                    ancestral_Ne = theta_syn / (4 * L_syn * simulation_mu)
                    f.write(
                        'Estimate for time is ' + str(generations) + 
                        ' generations.\n')
                    f.write(
                        'Estimated ancestral population size is ' +
                        str(ancestral_Ne) + '.\n')
                    extra_estimates['generations'] = generations
                else:
                    best_nu_b = best_params[0]
                    best_nu_f = best_params[1]
                    best_tau_b = best_params[2]
                    best_tau_f = best_params[3]
                    # Correct the algebra
                    generations_b = 2 * best_tau_b * theta_syn / (4 * simulation_mu * L_syn)
                    generations_f = 2 * best_tau_f * theta_syn / (4 * simulation_mu * L_syn)
                    generations_total = generations_b + generations_f
                    ancestral_Ne = theta_syn / (4 * L_syn * simulation_mu)
                    f.write(
                        'Estimate for bottleneck length in generations is ' +
                        str(generations_b) + ' generations.\n')
                    f.write(
                        'Estimate for time since bottleneck in generations is ' +
                        str(generations_f) + ' generations.\n')
                    f.write(
                        'Estimate for total time in generations is ' +
                        str(generations_total) + '.\n')
                    f.write(
                        'Estimate for ancestral population size is ' +
                        str(ancestral_Ne) + '.\n')
                    extra_estimates['generations_bottleneck'] = generations_b
                    extra_estimates['generations_since_bottleneck'] = \
                        generations_f
                    extra_estimates['generations_total'] = generations_total

            # Save the best-fit model SFS (non-scaled and theta-scaled)
            # for this model to disk in dadi's native Spectrum format.
            non_scaled_sfs_file = \
                '{0}{1}{2}_best_model_sfs_non_scaled.fs'.format(
                    args['outprefix'], underscore, model)
            scaled_sfs_file = \
                '{0}{1}{2}_best_model_sfs_scaled.fs'.format(
                    args['outprefix'], underscore, model)
            best_non_scaled_spectrum.to_file(non_scaled_sfs_file)
            best_scaled_spectrum.to_file(scaled_sfs_file)
            
            logger.info(
                'Saved best-fit model SFS for {0} model to {1} '
                'and {2}.\n'.format(model, non_scaled_sfs_file,
                                     scaled_sfs_file))

            # Append a row summarizing this model's fit to the CSV.
            csv_writer.writerow({
                'model': model,
                'best_params': json.dumps(list(numpy.atleast_1d(
                    best_params).astype(float))),
                'multinomial_log_likelihood': max_likelihood,
                'poisson_log_likelihood': poisson_ll,
                'theta_syn': theta_syn,
                'theta_nonsyn': theta_nonsyn,
                'L_syn': L_syn,
                'ancestral_Ne_estimate': ancestral_Ne,
                'extra_estimates': json.dumps(extra_estimates),
                'fold': fold,
                'non_scaled_sfs_file': non_scaled_sfs_file,
                'scaled_sfs_file': scaled_sfs_file,
            })
            csv_handle.flush()

        csv_handle.close()
        logger.info('Wrote CSV summary of all fitted models to {0}\n'.format(
            demography_results_csv))
        logger.info('Finished demographic inference.')
        logger.info('Pipeline executed succesfully.')


if __name__ == '__main__':
    DemographicInference().main()
