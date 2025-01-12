from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    POSTGRES_USER: str
    POSTGRES_PASSWORD: str
    POSTGRES_DB: str
    DB_HOST: str = "localhost"
    # openssl rand -hex 32
    JWT_SECRET_KEY: str = (
        "26192e0eba1df7e20106f2a2d55a65f1dde175361103721dd4b747cfbb5f2ab5"
    )
    JWT_ALGORITHM: str = "HS256"
    # No reply e-mail to send password reset link
    MAIL_USERNAME: str | None = None
    MAIL_PASSWORD: str | None = None
    MAIL_FROM: str | None = None
    MAIL_PORT: int | None = None
    MAIL_SERVER: str | None = None
    MAIL_FROM_NAME: str | None = None

    model_config = SettingsConfigDict(env_file=".env", extra="ignore")


settings = Settings()
