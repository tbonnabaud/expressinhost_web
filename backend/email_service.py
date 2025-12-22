import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

from .logger import logger
from .settings import settings


def send_email(email_to: str, subject: str, content: str):
    message = MIMEMultipart()
    message["From"] = settings.MAIL_ADDRESS
    message["To"] = email_to
    message["Subject"] = subject
    message.attach(MIMEText(content, "plain"))

    try:
        with smtplib.SMTP_SSL(settings.MAIL_SERVER, settings.MAIL_PORT) as server:
            server.login(settings.MAIL_USERNAME, settings.MAIL_PASSWORD)
            server.sendmail(settings.MAIL_ADDRESS, email_to, message.as_string())
            logger.info(f"E-mail sent to {email_to} with success")

    except Exception as e:
        logger.error(e)
