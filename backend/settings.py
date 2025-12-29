from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    POSTGRES_USER: str
    POSTGRES_PASSWORD: str
    POSTGRES_DB: str
    DB_HOST: str = "localhost"
    DB_PORT: int = 5432
    # openssl rand -hex 32
    JWT_SECRET_KEY: str
    JWT_ALGORITHM: str = "HS256"
    # Cookie security: set to False in development, True in production
    COOKIE_SECURE: bool = True
    # No reply e-mail to send password reset link
    MAIL_SERVER: str | None = None
    MAIL_PORT: int = 465
    MAIL_USERNAME: str | None = None
    MAIL_PASSWORD: str | None = None
    MAIL_ADDRESS: str | None = None
    REDIS_HOST: str = "localhost"
    REDIS_PORT: int = 6379

    model_config = SettingsConfigDict(env_file=".env", extra="ignore")


settings = Settings()
