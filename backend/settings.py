from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    DATABASE_URL: str
    # openssl rand -hex 32
    JWT_SECRET_KEY: str = (
        "26192e0eba1df7e20106f2a2d55a65f1dde175361103721dd4b747cfbb5f2ab5"
    )
    JWT_ALGORITHM: str = "HS256"

    model_config = SettingsConfigDict(env_file=".env")


settings = Settings()
