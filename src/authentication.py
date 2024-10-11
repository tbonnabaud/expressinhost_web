import bcrypt


def hash_password(password: str) -> str:
    """
    The function `hash_password` takes a password as input, hashes it using bcrypt, and returns the
    hashed password as a string.

    Args:
      password (str): The `hash_password` function takes a password string as input and returns the
    hashed version of the password using the bcrypt hashing algorithm. The `password` parameter is the
    string that you want to hash.

    Returns:
      A hashed version of the input password is being returned.
    """
    return bcrypt.hashpw(password.encode(), bcrypt.gensalt()).decode()


def check_password(password: str, hashed_password: str) -> bool:
    """
    The function `check_password` compares a plain text password with a hashed password using bcrypt to
    determine if they match.

    Args:
      password (str): The `password` parameter is a string that represents the plain text password that
    a user enters.
      hashed_password (str): The `hashed_password` parameter is the hashed version of the original
    password. It is the result of applying a cryptographic hash function to the password to securely
    store it or compare it with other hashed passwords for authentication purposes.

    Returns:
      A boolean value indicating whether the provided password matches the hashed password.
    """
    return bcrypt.checkpw(password.encode(), hashed_password.encode())
