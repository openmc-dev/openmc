from numbers import Real
from xml.etree import ElementTree as ET
import sys
import warnings
from collections.abc import Iterable

import openmc.checkvalue as cv


class Trigger:
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
        self.trigger_type = trigger_type
        self.threshold = threshold
        self._scores = []

    def __eq__(self, other):
        return str(self) == str(other)

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
        cv.check_type('trigger scores', scores, Iterable, str)

        # Set scores making sure not to have duplicates
        self._scores = []
        for score in scores:
            if score not in self._scores:
                self._scores.append(score)

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
