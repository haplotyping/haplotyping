import logging,requests,re
from haplotyping.graph.graph import Graph
import haplotyping

class Sections():
    """constructing linear ordered sections from candidates in the De Bruijn graph"""
    
    def __init__(self, graph: Graph):
        """
        construct sections
        """
        assert isinstance(graph, Graph)
        #logger
        self._logger = logging.getLogger(__name__)