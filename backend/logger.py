import logging
import sys
from logging.handlers import RotatingFileHandler

uvicorn_logger = logging.getLogger("uvicorn.error")
logger = logging.getLogger("expressinhost")
scraping_logger = logging.getLogger("expressinhost.scraping")

# Set the logging level
logger.setLevel(logging.INFO)

# Create a rotating file handler
file_handler = RotatingFileHandler(
    "logs/expressinhost.log",  # Log file name
    maxBytes=5 * 1024 * 1024,  # Maximum file size (5 MB)
    backupCount=2,  # Number of backup files to keep
)

# Create a formatter and set it for the handler
# formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
file_handler.setFormatter(formatter)
stdout_handler = logging.StreamHandler(sys.stdout)
stdout_handler.setFormatter(formatter)

# Add the file handler to the logger
logger.addHandler(file_handler)
# Add the stdout handler to the logger
logger.addHandler(stdout_handler)
