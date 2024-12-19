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


    @property
    def name(self):
        return self._name


    @property
    def length(self):
        return self._length


    @property
    def seq(self):
        return self._ref_seq

