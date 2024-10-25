import bcrypt


def hash_password(password: str) -> str:
    """
    The function `hash_password` takes a password as input, hashes it using bcrypt, and returns the
    hashed password as a string.
    """
    return bcrypt.hashpw(password.encode(), bcrypt.gensalt()).decode()


def check_password(password: str, hashed_password: str) -> bool:
    """
    The function `check_password` compares a plain text password with a hashed password using bcrypt to
    determine if they match.
    """
    return bcrypt.checkpw(password.encode(), hashed_password.encode())
