#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Exceptions.

COPYRIGHT

@author John D. Chodera <jchodera@gmail.com>

"""

#=============================================================================================
# EXCEPTIONS
#=============================================================================================

class UnitsException(Exception):
    """
    Exception denoting that an argument has the incorrect units.

    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class ValueException(Exception):
    """
    Exception denoting that an argument has the incorrect value.

    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

