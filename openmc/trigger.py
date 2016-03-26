from numbers import Real
from xml.etree import ElementTree as ET
import sys
import warnings
from collections import Iterable

import openmc.checkvalue as cv

if sys.version_info[0] >= 3:
    basestring = str


class Trigger(object):
    """A criterion for when to finish a simulation based on tally uncertainties.

    Parameters
    ----------
    trigger_type : {'variance', 'std_dev', 'rel_err'}
        Determine whether to trigger on the variance, standard deviation, or
        relative error of scores.
    threshold : float
        The threshold for the trigger type.

    Attributes
    ----------
    trigger_type : {'variance', 'std_dev', 'rel_err'}
        Determine whether to trigger on the variance, standard deviation, or
        relative error of scores.
    threshold : float
        The threshold for the trigger type.
    scores : list of str
        Scores which should be checked against the trigger

    """

    def __init__(self, trigger_type, threshold):
        # Initialize Mesh class attributes
        self.trigger_type = trigger_type
        self.threshold = threshold
        self._scores = []

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is first time we have tried to copy this object, create a copy
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._trigger_type = self._trigger_type
            clone._threshold = self._threshold

            clone.scores = self.scores

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    def __eq__(self, other):
        if str(self) == str(other):
            return True
        else:
            return False

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        string = 'Trigger\n'
        string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._trigger_type)
        string += '{0: <16}{1}{2}\n'.format('\tThreshold', '=\t', self._threshold)
        string += '{0: <16}{1}{2}\n'.format('\tScores', '=\t', self._scores)
        return string

    @property
    def trigger_type(self):
        return self._trigger_type

    @property
    def threshold(self):
        return self._threshold

    @property
    def scores(self):
        return self._scores

    @trigger_type.setter
    def trigger_type(self, trigger_type):
        cv.check_value('tally trigger type', trigger_type,
                    ['variance', 'std_dev', 'rel_err'])
        self._trigger_type = trigger_type

    @threshold.setter
    def threshold(self, threshold):
        cv.check_type('tally trigger threshold', threshold, Real)
        self._threshold = threshold

    @scores.setter
    def scores(self, scores):
        cv.check_type('trigger scores', scores, Iterable, basestring)

        # Set scores making sure not to have duplicates
        self._scores = []
        for score in scores:
            if score not in self._scores:
                self._scores.append(score)


    def add_score(self, score):
        """Add a score to the list of scores to be checked against the trigger.

        Parameters
        ----------
        score : str
            Score to append

        """

        warnings.warn('Trigger.add_score(...) has been deprecated and may be '
                      'removed in a future version. Tally trigger scores should '
                      'be defined using the scores property directly.',
                      DeprecationWarning)
        self.scores.append(score)

    def get_trigger_xml(self, element):
        """Return XML representation of the trigger

        Returns
        -------
        element : xml.etree.ElementTree.Element
            XML element containing trigger data

        """

        subelement = ET.SubElement(element, "trigger")
        subelement.set("type", self._trigger_type)
        subelement.set("threshold", str(self._threshold))
        if len(self._scores) != 0:
            subelement.set("scores", ' '.join(map(str, self._scores)))
