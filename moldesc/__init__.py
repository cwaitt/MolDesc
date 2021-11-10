"""ASE integrated modules to compute physical quantites of molecules (volume, elongation, etc.)"""

# Add imports here
from .moldesc import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
