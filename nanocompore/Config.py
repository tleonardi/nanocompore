from schema import Schema, And, Or, Optional
from nanocompore.common import is_valid_fasta


CONFIG_SCHEMA = Schema({
    'data': And(
        {
            str: { # condition
                str: { # sample (replicate)
                    'bam': lambda f: open(f, 'r'),
                    'pod5': lambda f: open(f, 'r')
                }
            }
        },
        And(lambda d: len(d) == 2, error='Only two conditions allowed')
    ),
    'fasta': And(is_valid_fasta, error='Invalid fasta file'),
    Optional('kit'): Or('RNA002', 'RNA004'),
    Optional('bed'): And(lambda f: open(f, 'r'), error='Invalid bed file'),
    Optional('nthreads'): And(lambda n: n >= 3, error='nthreads must be >= 3'),
    Optional('min_coverage'): And(int, lambda n: n >= 0, error='min_coverage must be >= 0'),
    Optional('downsample_high_coverage'): And(int, lambda n: n >= 0, error='downsample_high_coverage must be >= 0'),
    Optional('min_ref_length'): And(int, lambda n: n >= 0, error='min_ref_length must be >= 0'),
    Optional('comparison_methods'): And(['GMM',
                                         'KS',
                                         'TT',
                                         'MW',
                                         'GAUSSIAN_MIXTURE_MODEL',
                                         'KOLMOGOROV_SMIRNOV',
                                         'T_TEST',
                                         'MANN_WHITNEY'],
                                        lambda l: len(l) > 0),
    Optional('sequence_context'): And(int, lambda n: n >= 0 and n <= 4, error='sequence_context must be >= 0 and <= 4'),
    Optional('sequence_context_weights'): Or('uniform', 'harmonic'),
    Optional('pvalue_threshold'): And(float, lambda n: n >= 0 and n <= 1, error='pvalue_threshold must be >= 0 and <= 1'),
    Optional('logit'): bool,
    Optional('anova'): bool,
    Optional('bool'): bool,
    Optional('allow_warnings'): bool,
    Optional('outpath'): str,
    Optional('outprefix'): str,
    Optional('overwrite'): bool,
    Optional('log_level'): Or('warning', 'info', 'debug'),
    Optional('progress'): bool,
    Optional('correction_method'): 'fdr_bh',
})


DEFAULT_KIT = 'RNA002'
DEFAULT_NTHEARDS = 3
DEFAULT_MIN_COVERAGE = 30
DEFAULT_MAX_READS = 5000
DEFAULT_DOWNSAMPLE_HIGH_COVERAGE = 5000
DEFAULT_MIN_REF_LENGTH = 100
DEFAULT_COMPARISON_METHODS = ['GMM', 'KS']
DEFAULT_SEQUENCE_CONTEXT = 0
DEFAULT_SEQUENCE_CONTEXT_WEIGHTS = 'uniform'
DEFAULT_PVALUE_THRESHOLD = 0.05
DEFAULT_LOGIT = True
DEFAULT_ANOVA = False
DEFAULT_ALLOW_WARNINGS = False
DEFAULT_OUTPATH = 'nanocompore_output'
DEFAULT_OUTPREFIX = 'out_'
DEFAULT_OVERWRITE = False
DEFAULT_LOG_LEVEL = 'info'
DEFAULT_PROGRESS = False
DEFAULT_CORRECTION_METHOD = 'fdr_bh'


class Config:
    """
    Represents the input configuration for a Nanocompore experiment.
    An instance contains information about the input data files and
    any additoinal parameters for Nanocompore.
    """

    def __init__(self, config_file):
        self._config = CONFIG_SCHEMA.validate(config_file)


    def get_data(self):
        return self._config['data']


    def get_kit(self):
        return self._config.get('kit', DEFAULT_KIT)


    def get_nthreads(self):
        return self._config.get('nthreads', DEFAULT_NTHEARDS)


    def get_fasta_ref(self):
        """
        Reference fasta file that was used for mapping the data.
        """
        return self._config['fasta']


    def get_bed(self):
        """
        BED file with annotation of transcriptome used for mapping.
        """
        return self._config.get('bed')


    def get_min_coverage(self):
        """
        Minimum coverage required in each condition to do the comparison.
        """
        return self._config.get('min_coverage', DEFAULT_MIN_COVERAGE)


    def get_downsample_high_coverage(self):
        """
        Transcripts with high coverage will be downsampled.
        """
        return self._config.get('downsample_high_coverage', DEFAULT_DOWNSAMPLE_HIGH_COVERAGE)


    def get_min_ref_length(self):
        """
        Minimum length of a reference transcript to include it in the analysis.
        """
        return self._config.get('min_ref_length', DEFAULT_MIN_REF_LENGTH)


    def get_comparison_methods(self):
        """
        List of comparison methods. Valid methods are: GMM,KS,TT,MW.
        """
        remappings = {
            'GAUSSIAN_MIXTURE_MODEL': 'GMM',
            'KOLMOGOROV_SMIRNOV': 'KS',
            'T_TEST': 'TT',
            'MANN_WHITNEY': 'MW',
        }
        return [remappings.get(m, m)
                for m in self._config.get('comparison_methods', DEFAULT_COMPARISON_METHODS)]


    def get_sequence_context(self):
        """
        Sequence context for combining p-values.
        """
        return self._config.get('sequence_context', DEFAULT_SEQUENCE_CONTEXT)


    def get_sequence_context_weights(self):
        """
        Type of weights to use for combining p-values.
        """
        return self._config.get('sequence_context_weights', DEFAULT_SEQUENCE_CONTEXT_WEIGHTS)


    def get_pvalue_threshold(self):
        """
        Adjusted p-value threshold for reporting significant sites.
        """
        return self._config.get('pvalue_thr', DEFAULT_PVALUE_THRESHOLD)


    def get_logit(self):
        """
        Use logistic regression testing downstream of GMM method.
        """
        return self._config.get('logit', DEFAULT_LOGIT)


    def get_anova(self):
        """
        Use Anova test downstream of GMM method.
        """
        return self._config.get('anova', DEFAULT_ANOVA)


    def get_allow_warnings(self):
        """
        If True runtime warnings during the ANOVA tests don't raise an error.
        """
        return self._config.get('allow_warnings', DEFAULT_ALLOW_WARNINGS)


    def get_outpath(self):
        """
        Path to the output folder.
        """
        return self._config.get('outpath', DEFAULT_OUTPATH)


    def get_outprefix(self):
        """
        Text outprefix for all the files generated.
        """
        prefix = self._config.get('outprefix', DEFAULT_OUTPREFIX)
        if prefix and prefix.endswith('_'):
            return prefix
        return prefix + '_'


    def get_overwrite(self):
        """
        Use <outpath> even if it exists already.
        """
        return self._config.get('overwrite', DEFAULT_OVERWRITE)


    def get_log_level(self):
        """
        Get the log level.
        """
        return self._config.get('log_level', DEFAULT_LOG_LEVEL)


    def get_progress(self):
        """
        Display a progress bar during execution.
        """
        return self._config.get('progress', DEFAULT_PROGRESS)


    def get_correction_method(self):
        """
        Get the multiple test correction method.
        """
        return self._config.get('correction_method', DEFAULT_CORRECTION_METHOD)
