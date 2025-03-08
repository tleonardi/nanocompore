import re
import os

from schema import Schema, And, Or, Optional

from nanocompore.common import is_valid_fasta
from nanocompore.common import Kit


def valid_device(value):
    """
    Check if a device is valid.

    Parameters
    ----------
    value : str
        Device name

    Returns
    -------
        bool
    """
    return value in ['cpu', 'cuda'] or re.compile(r"cuda:\d+").match(value)


def validate_device(value):
    """
    Validate the device configuration.

    Parameters
    ----------
    value : str | dict
        Device configuration.

    Returns
    -------
        bool
    """
    if isinstance(value, str):
        return valid_device(value)

    if isinstance(value, dict):
        for k, v in value.items():
            if not valid_device(k) or not isinstance(v, int) or v < 0:
                return False
        return True
    raise ValueError("The value for 'devices' in the configuration "
                     f"is of unexpected type: {type(value)}")


def validate_eventalign_data(config):
    """
    Validate if the input data for eventalign contains
    all necessary fields.

    Parameters
    ----------
    config : dict
        The configuration object.

    Returns
    -------
        bool
    """
    if config['resquiggler'] != 'eventalign':
        return True
    for condition, cond_data in config['data'].items():
        for sample, sample_data in cond_data.items():
            if ('eventalign_tsv' not in sample_data and
                'eventalign_db' not in sample_data):
                return False
    return True


def validate_uncalled4_data(config):
    """
    Validate if the input data for Uncalled4 contains
    all necessary fields.

    Parameters
    ----------
    config : dict
        The configuration object.

    Returns
    -------
        bool
    """
    if config['resquiggler'] != 'uncalled4':
        return True
    for condition, cond_data in config['data'].items():
        for sample, sample_data in cond_data.items():
            if 'bam' not in sample_data:
                return False
    return True


def depleted_condition_exists(config):
    return config['depleted_condition'] in config['data']


CONFIG_SCHEMA = Schema(And({
    'data': And(
        {
            str: { # condition
                str: { # sample (replicate)
                    Optional('bam'): lambda f: open(f, 'r'),
                    Optional('pod5'): str,
                    Optional('eventalign_tsv'): lambda f: open(f, 'r'),
                    Optional('eventalign_db'): lambda f: open(f, 'r')
                }
            }
        },
        And(lambda d: len(d) == 2, error='Only two conditions allowed')
    ),
    'depleted_condition': str,
    'fasta': And(is_valid_fasta, error='Invalid fasta file'),
    'resquiggler': Or('uncalled4', 'eventalign', 'remora'),
    'kit': Or(*[v.name for v in Kit]),
    Optional('devices'): validate_device,
    Optional('preprocessing_db'): str,
    Optional('bed'): And(lambda f: open(f, 'r'), error='Invalid bed file'),
    Optional('nthreads'): And(lambda n: n >= 2, error='nthreads must be >= 2'),
    Optional('min_coverage'): And(int, lambda n: n >= 0, error='min_coverage must be >= 0'),
    Optional('downsample_high_coverage'): And(int, lambda n: n >= 0, error='downsample_high_coverage must be >= 0'),
    Optional('min_ref_length'): And(int, lambda n: n >= 0, error='min_ref_length must be >= 0'),
    Optional('max_invalid_kmers_freq'): And(float, lambda n: n >= 0, error='max_invalid_kmers_freq must be in the [0, 1] range'),
    Optional('comparison_methods'): And(['GMM',
                                         'GOF',
                                         'KS',
                                         'TT',
                                         'MW',
                                         'auto',
                                         'AUTO',
                                         'GAUSSIAN_MIXTURE_MODEL',
                                         'GOODNESS_OF_FIT',
                                         'KOLMOGOROV_SMIRNOV',
                                         'T_TEST',
                                         'MANN_WHITNEY'],
                                        lambda l: len(l) > 0),
    Optional('motor_dwell_offset'): And(int, lambda n: n >= 0, error='motor_dwell_offset must be >= 0'),
    Optional('sequence_context'): And(int, lambda n: n >= 0 and n <= 4, error='sequence_context must be >= 0 and <= 4'),
    Optional('sequence_context_weights'): Or('uniform', 'harmonic'),
    Optional('pvalue_threshold'): And(float, lambda n: n >= 0 and n <= 1, error='pvalue_threshold must be >= 0 and <= 1'),
    Optional('logit'): bool,
    Optional('anova'): bool,
    Optional('bool'): bool,
    Optional('allow_warnings'): bool,
    Optional('outpath'): str,
    Optional('outprefix'): str,
    Optional('result_exists_strategy'): Or("stop", "continue", "overwrite"),
    Optional('log_level'): Or('warning', 'info', 'debug'),
    Optional('progress'): bool,
    Optional('correction_method'): 'fdr_bh'},
    # Additional validation of the full configuration
    And(validate_eventalign_data,
        error='When using the "eventalign" resquiggler ' +
              'each sample in the configuration must contain ' +
              'either the field "eventalign_tsv" with a path ' +
              'to the eventalign tsv file or "eventalign_db" ' +
              'with the alreday collapsed eventalign data.'),
    And(validate_uncalled4_data,
        error='When using the "uncalled4" resquiggler ' +
              'each sample in the configuration must contain ' +
              'the field "bam" with a path to the aligned bam file ' +
              'that contains the resquiggling tags produced ' +
              'by Uncalled4.'),
    And(depleted_condition_exists,
        error="The condition set in 'depleted_condition' is not " + \
              "defined in 'data'.")))


DEFAULT_KIT = 'RNA002'
DEFAULT_DEVICES = 'cpu'
DEFAULT_PREPROCESSING_DB = 'preprocessing_db.sqlite'
DEFAULT_NTHEARDS = 2
DEFAULT_MIN_COVERAGE = 30
DEFAULT_MAX_READS = 5000
DEFAULT_DOWNSAMPLE_HIGH_COVERAGE = 5000
DEFAULT_MIN_REF_LENGTH = 100
DEFAULT_MAX_INVALID_KMERS_FREQ = 0.1
DEFAULT_COMPARISON_METHODS = ['GMM', 'KS']
DEFAULT_MOTOR_DWELL_OFFSET = 0
DEFAULT_SEQUENCE_CONTEXT = 0
DEFAULT_SEQUENCE_CONTEXT_WEIGHTS = 'uniform'
DEFAULT_PVALUE_THRESHOLD = 0.05
DEFAULT_LOGIT = True
DEFAULT_ANOVA = False
DEFAULT_ALLOW_WARNINGS = False
DEFAULT_OUTPATH = 'nanocompore_output'
DEFAULT_OUTPREFIX = 'out_'
DEFAULT_RESULT_EXISTS_STRATEGY = 'stop'
DEFAULT_LOG_LEVEL = 'info'
DEFAULT_PROGRESS = False
DEFAULT_CORRECTION_METHOD = 'fdr_bh'
DEFAULT_READ_LEVEL_DATA = False
DEFAULT_READ_LEVEL_DATA_TRANSCRIPTS = []


class Config:
    """
    Represents the input configuration for a Nanocompore experiment.
    An instance contains information about the input data files and
    any additoinal parameters for Nanocompore.
    """

    def __init__(self, config_file):
        self._config = CONFIG_SCHEMA.validate(config_file)

        self._test_condition = [cond
                                for cond in self.get_condition_labels()
                                if cond != self.get_depleted_condition()][0]


    def get_data(self):
        return self._config['data']


    def get_preprocessing_db(self):
        """
        Returns the path where the kmer data should
        be stored after resquiggling/preprocessing.
        """
        path = self._config.get('preprocessing_db', DEFAULT_PREPROCESSING_DB)
        if os.path.isabs(path):
            return path
        else:
            return os.path.join(self.get_outpath(), path)


    def get_resquiggler(self):
        """
        Returns the resquiggler software specified
        by the user.
        """
        return self._config['resquiggler']


    def get_kit(self):
        """
        Returns an instance of the Kit enum representing
        the sequencing kit, e.g. RNA002.
        """
        kit_name = self._config.get('kit', DEFAULT_KIT)
        return Kit[kit_name]


    def get_devices(self):
        """
        Returns the device to be used for computations.
        E.g. cpu or cuda.
        """
        return self._config.get('devices', DEFAULT_DEVICES)


    def get_nthreads(self):
        """
        Number of threads to use.
        """
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


    def get_max_invalid_kmers_freq(self):
        """
        Maximum allowed ratio of invalid kmers in the read.
        """
        return self._config.get('max_invalid_kmers_freq', DEFAULT_MAX_INVALID_KMERS_FREQ)


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


    def get_motor_dwell_offset(self):
        """
        Interactions between a modification and the motor protein
        can affect the translocation speed. If present, this
        should affect the dwell time of the downstream position.
        This returns the expected offset for this effect.
        """
        return self._config.get('motor_dwell_offset', DEFAULT_MOTOR_DWELL_OFFSET)


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


    def get_result_exists_strategy(self):
        """
        What to do if <outpath> already exists.
        Options are: stop, continue, and overwrite.
        Default: stop
        """
        return self._config.get('result_exists_strategy', DEFAULT_RESULT_EXISTS_STRATEGY)


    def get_log_level(self):
        """
        Get the log level.
        """
        return self._config.get('log_level', DEFAULT_LOG_LEVEL).upper()


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


    def get_read_level_data(self):
        """
        Save read level data to database.
        (That means the intensity/dwell-time measurements
        for every position of every read will be stored.)
        Warning: This can take up a lot of space.
        """
        return self._config.get('read_level_data', DEFAULT_READ_LEVEL_DATA)


    def get_read_level_data_transcripts(self):
        """
        If read_level_data is True, only save read level data for these transcripts.
        The transcript ids must match the RNAME field from the bam file.
        Note: empty list means all transcripts
        """
        return self._config.get('read_level_data_transcripts', DEFAULT_READ_LEVEL_DATA_TRANSCRIPTS)


    def get_sample_labels(self):
        return [sample
                for samples in self.get_data().values()
                for sample in samples]


    def get_condition_labels(self):
        return list(self.get_data().keys())


    def get_depleted_condition(self):
        return self._config['depleted_condition']


    def get_test_condition(self):
        return self._test_condition


    def is_multi_replicate(self):
        return all(len(reps) > 1
                   for reps in self.get_data().values())


    def get_sample_ids(self):
        labels = self.get_sample_labels()
        return dict(zip(labels, range(len(labels))))


    def get_condition_ids(self):
        labels = self.get_condition_labels()
        return dict(zip(labels, range(len(labels))))


    def sample_to_condition(self):
        return {sample: cond
                for cond, sample_defs in self.get_data().items()
                for sample in sample_defs.keys()}


    def get_condition_samples(self, condition_label):
        return list(self.get_data()[condition_label].keys())


    def get_conditions_to_samples(self):
        return {cond: list(sample_defs.keys())
                for cond, sample_defs in self.get_data().items()}


    def get_sample_condition_bam_data(self):
        return [(sample, condition, samp_def['bam'])
                for condition, samples in self.get_data().items()
                for sample, samp_def in samples.items()]


    def get_sample_pod5_bam_data(self):
        return [(sample, samp_def['pod5'], samp_def['bam'])
                for _, samples in self.get_data().items()
                for sample, samp_def in samples.items()]

