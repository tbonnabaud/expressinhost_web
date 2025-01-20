import logging
from logging.handlers import RotatingFileHandler

logger = logging.getLogger("uvicorn.error")

# Set the logging level
logger.setLevel(logging.WARNING)

# Create a rotating file handler
handler = RotatingFileHandler(
    "expressinhost.log",  # Log file name
    maxBytes=5 * 1024 * 1024,  # Maximum file size (5 MB)
    backupCount=2,  # Number of backup files to keep
)

# Create a formatter and set it for the handler
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
# formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)

# Add the handler to the logger
logger.addHandler(handler)
