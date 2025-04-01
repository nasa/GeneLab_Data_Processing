from pathlib import Path

# Import for access at the module level
from . import checks
from . import protocol
from . import schemas

# Set config path
config = Path(__file__).parent / "config.yaml"