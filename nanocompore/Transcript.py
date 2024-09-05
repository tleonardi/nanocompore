from loguru import logger


class Transcript():
    def __init__(self,
                 ref_id='',
                 experiment='',
                 config=None,
                 ref_seq=0):

        logger.trace(f"Creating transcript object for {ref_id}")
        self._experiment = experiment
        self._name = ref_id
        self._ref_seq = ref_seq
        self._length = len(self._ref_seq)


    ########## Public ##########

    @property
    def name(self):
        return self._name


    @property
    def length(self):
        return self._length


    @property
    def seq(self):
        return self._ref_seq


    @property
    def condition_labels(self):
        return list(self._experiment.get_condition_labels())


    @property
    def sample_labels(self):
        return self._experiment.get_sample_labels()


    @property
    def sample_2_condition(self, sample_label):
        return self._experiment.sample_to_condition(sample_label)

