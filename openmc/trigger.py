from collections.abc import Iterable
from numbers import Real

import lxml.etree as ET

import openmc.checkvalue as cv
from .mixin import EqualityMixin


class Trigger(EqualityMixin):
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

    def __init__(self, trigger_type: str, threshold: float):
        self.trigger_type = trigger_type
        self.threshold = threshold
        self._scores = []

    def __repr__(self):
        string = 'Trigger\n'
        string += '{: <16}=\t{}\n'.format('\tType', self._trigger_type)
        string += '{: <16}=\t{}\n'.format('\tThreshold', self._threshold)
        string += '{: <16}=\t{}\n'.format('\tScores', self._scores)
        return string

    @property
    def trigger_type(self):
        return self._trigger_type

    @trigger_type.setter
    def trigger_type(self, trigger_type):
        cv.check_value('tally trigger type', trigger_type,
                       ['variance', 'std_dev', 'rel_err'])
        self._trigger_type = trigger_type

    @property
    def threshold(self):
        return self._threshold

    @threshold.setter
    def threshold(self, threshold):
        cv.check_type('tally trigger threshold', threshold, Real)
        self._threshold = threshold

    @property
    def scores(self):
        return self._scores

    @scores.setter
    def scores(self, scores):
        cv.check_type('trigger scores', scores, Iterable, str)

        # Set scores making sure not to have duplicates
        self._scores = []
        for score in scores:
            if score not in self._scores:
                self._scores.append(score)

    def to_xml_element(self):
        """Return XML representation of the trigger

        Returns
        -------
        element : lxml.etree._Element
            XML element containing trigger data

        """

        element = ET.Element("trigger")
        element.set("type", self._trigger_type)
        element.set("threshold", str(self._threshold))
        if len(self._scores) != 0:
            element.set("scores", ' '.join(self._scores))
        return element

    @classmethod
    def from_xml_element(cls, elem: ET.Element):
        """Generate trigger object from an XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element

        Returns
        -------
        openmc.Trigger
            Trigger object

        """
        # Generate trigger object
        trigger_type = elem.get("type")
        threshold = float(elem.get("threshold"))
        trigger = cls(trigger_type, threshold)

        # Add scores if present
        scores = elem.get("scores")
        if scores is not None:
            trigger.scores = scores.split()

        return trigger
