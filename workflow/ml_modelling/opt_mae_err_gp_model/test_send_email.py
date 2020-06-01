# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# +
# Import smtplib for the actual sending function
import smtplib

# Import the email modules we'll need
from email.mime.text import MIMEText

# Open a plain text file for reading.  For this example, assume that
# the text file contains only ASCII characters.
# textfile = "textfile"
# with open(textfile, 'rb') as fp:
#     # Create a text/plain message
#     msg = MIMEText(fp.read())

msg = MIMEText(
"""
Hello Bello
"""
)

# msg = """
# Hello bello.
# """

# me == the sender's email address
# you == the recipient's email address
# msg['Subject'] = 'The contents of %s' % textfile
msg['Subject'] = 'The contents of TEMP'
# msg['From'] = "raulf2012@gmail.com"
msg['From'] = "raulf2012@hotmail.com"
msg['To'] = "raulf2012@hotmail.com"

# Send the message via our own SMTP server, but don't include the
# envelope header.

s = smtplib.SMTP('localhost')
s.sendmail(me, [you], msg.as_string())
s.quit()
