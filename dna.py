"""
Author: Yixin Zhang
Last updated: Dec. 11, 2023

This file contains DNA definitions.
"""

class DNA:
    def __init__(self, name, sequence, oligo_type):
        self.name = name
        self.sequence = str(sequence)
        self.ext5 = None
        self.ext3 = None
        self.is_double_stranded = True
        self.is_circular = False
        self.mod_ext5 = None
        self.mod_ext3 = None

        if oligo_type == 'primer':
            self.mod_ext5 = 'hydroxyl'
            self.is_double_stranded = False
        if oligo_type == 'plasmid':
            self.is_circular = True
    
    def get_json_dict(self):
        """
        Generates a dictionary of DNA information for json interpretation.

        Returns:
        dict: DNA information.
        """
        return {"sequence": self.sequence,
                "ext5": self.ext5,
                "ext3": self.ext3,
                "is_double_stranded": self.is_double_stranded,
                "is_circular": self.is_circular,
                "mod_ext5": self.mod_ext5,
                "mod_ext3": self.mod_ext3
                }
