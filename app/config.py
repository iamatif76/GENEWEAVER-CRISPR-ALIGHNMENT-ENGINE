"""
Application configuration for GeneWeaver.
Loads environment variables securely using python-dotenv.
"""

from pathlib import Path
from pydantic_settings import BaseSettings, SettingsConfigDict


BASE_DIR = Path(__file__).resolve().parent.parent


class Settings(BaseSettings):
    # Security: Read secrets exclusively from environment variables
    API_KEY: str = ""

    # Server settings
    HOST: str = "0.0.0.0"
    PORT: int = 8000
    DEBUG: bool = False
    LOG_LEVEL: str = "INFO"

    # Database settings
    DATABASE_URL: str = "sqlite:///./geneweaver.db"

    # File storage
    UPLOAD_DIR: Path = BASE_DIR / "uploads"
    DATA_DIR: Path = BASE_DIR / "data"

    # CRISPR defaults
    DEFAULT_PAM: str = "NGG"
    DEFAULT_MAX_MISMATCHES: int = 3
    DEFAULT_GUIDE_LENGTH: int = 20
    DEFAULT_SEED_LENGTH: int = 10  # PAM-proximal seed region (3' end of protospacer)

    # Hardware acceleration mode: "auto", "cpu", "gpu"
    GPU_ENABLED: str = "auto"

    model_config = SettingsConfigDict(
        env_file=str(BASE_DIR / ".env"),
        env_file_encoding="utf-8",
        extra="ignore"
    )


settings = Settings()

# Ensure required directories exist
settings.UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
settings.DATA_DIR.mkdir(parents=True, exist_ok=True)
