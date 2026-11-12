"""
Security and authentication dependencies for GeneWeaver.
Validates API keys securely using constant-time comparison.
Never logs or exposes secret credentials.
"""

import hmac
import logging
from typing import Optional
from fastapi import Header, HTTPException, Security, status
from fastapi.security import APIKeyHeader, HTTPBearer, HTTPAuthorizationCredentials

from app.config import settings

logger = logging.getLogger(__name__)

api_key_header_scheme = APIKeyHeader(name="X-API-Key", auto_error=False)
bearer_scheme = HTTPBearer(auto_error=False)


def verify_api_key(
    x_api_key: Optional[str] = Security(api_key_header_scheme),
    auth_credentials: Optional[HTTPAuthorizationCredentials] = Security(bearer_scheme),
) -> bool:
    """
    Authenticate protected endpoints using API key provided via X-API-Key or Bearer token.
    Reads expected key from environment configuration.
    """
    configured_key = settings.API_KEY.strip()

    # If no API key is configured in environment, permit open access
    if not configured_key:
        return True

    provided_key = None
    if x_api_key:
        provided_key = x_api_key.strip()
    elif auth_credentials and auth_credentials.credentials:
        provided_key = auth_credentials.credentials.strip()

    if not provided_key:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Missing API key. Provide via 'X-API-Key' header or 'Authorization: Bearer <key>'.",
            headers={"WWW-Authenticate": "Bearer"},
        )

    # Use constant-time comparison to prevent timing side-channel attacks
    if not hmac.compare_digest(provided_key, configured_key):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid API key provided.",
            headers={"WWW-Authenticate": "Bearer"},
        )

    return True

