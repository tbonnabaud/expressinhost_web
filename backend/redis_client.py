from redis import Redis

from .settings import settings

# Create a Redis connection for token storage
redis_client = Redis(
    host=settings.REDIS_HOST,
    port=settings.REDIS_PORT,
    decode_responses=True,  # Automatically decode responses to strings
)
